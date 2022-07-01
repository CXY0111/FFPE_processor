import vcf
import pysam
import re
import os
import subprocess
import argparse


def forward_backward(vcf_file, bam_file, output_file, argref, argalt, validate, count_orphans, min_BQ):
    #  header should be:
    #       chrom, pos, ref, alt, total_read, supported_forward, supported_reversed
    with open(output_file, 'w') as f:
        f.write("chrom\tpos\tref\talt\ttotal_read\tsupported_forward\tsupported_reversed\n")

    # open the vcf file and bam file
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    samfile = pysam.AlignmentFile(bam_file, "rb")
    for record in vcf_reader:
        chrom = record.CHROM
        pos = record.POS
        ref = record.REF
        alt = record.ALT

        # choose only SNP
        if len(ref) != 1 or len(alt) != 1 or len(alt[0]) != 1:
            continue

        # choose variant only C to T
        if ref != argref or alt[0] != argalt:
            continue

        # print out the record
        print('for', record, ':')

        forward_support, reverse_support = 0, 0
        total = 0

        # get the aligned base in a given position
        for pileupcolumn in samfile.pileup(chrom, pos - 1, pos, min_base_quality=min_BQ, ignore_orphans=not count_orphans):
            if pileupcolumn.pos == pos - 1:
                for pileupread in pileupcolumn.pileups:
                    # query position is None if is_del or is_refskip is set
                    if pileupread.query_position is None:
                        # print(pileupread)
                        total += 1
                        continue
                    base = pileupread.alignment.query_sequence[pileupread.query_position]
                    if base == alt[0]:
                        if pileupread.alignment.flag & 16:
                            reverse_support += 1
                        else:
                            forward_support += 1

                    total += 1
        print('backward_suport=', reverse_support, 'forward_suport=', forward_support, 'total=', total)
        with open(output_file, 'a') as f:
            f.write(chrom + "\t" + str(pos) + "\t" + ref + "\t" + str(alt[0]) + "\t" + str(total) + "\t" + str(
                forward_support) + "\t" + str(reverse_support) + "\n")

        if validate:
            # this is the part to validate the res by samtools
            if not count_orphans:
                cmd = 'samtools mpileup -f /diskmnt/Datasets_public/Reference/GRCh38.d1.vd1/GRCh38.d1.vd1.fa -s ' + \
                      bam_file + ' -r chr1:' + str(pos) + '-' + str(pos) + ' -Q '+ str(min_BQ) +' | cut -f 1,2,3,4,5 '
            else:
                cmd = 'samtools mpileup -f /diskmnt/Datasets_public/Reference/GRCh38.d1.vd1/GRCh38.d1.vd1.fa -s ' + \
                      bam_file + ' -r chr1:' + str(pos) + '-' + str(pos) + ' -Q '+ str(min_BQ) +' -A | cut -f 1,2,3,4,5 '
            ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            output = ps.communicate()[0]
            res = output.decode('ascii').split('\n')[1]
            result = res.split('\t')[4]

            forward, backward, other = 0, 0, 0

            # get rid of inserts and deletes
            new_res = ''
            i = 0
            while i < len(result):
                if result[i] == '+' or result[i] == '-':
                    if 48 <= ord(result[i + 1]) <= 57:
                        i += int(result[i + 1]) + 1
                else:
                    new_res += result[i]
                i += 1
            # get rid of begining and ending character
            new_res = re.sub("\^.|\$", '', new_res)

            for s in new_res:

                if s == str(alt[0]):
                    forward += 1
                elif s == str(alt[0]).lower():
                    backward += 1
                elif s == ',' or 97 <= ord(s) <= 122 or s == '.' or 65 <= ord(s) <= 90 or s == '*':
                    other += 1
            print('samtools_mpileup:')
            if backward == reverse_support and forward == forward_support and total == forward + backward + other:
                print(res, 'matches!')
            else:
                print(res, 'not match')

        print('')


if __name__ == '__main__':
    # to read command line arguments
    parser = argparse.ArgumentParser(description="Calculate the Reverse Strand proportion from a vcf file")

    parser.add_argument("-i", "--input", type=str, dest="vcf_file", help="input vcf file name")
    parser.add_argument("-b", "--bam", type=str, dest="bam_file", help="input bam file name")
    parser.add_argument("-o", "--output", type=str, dest="output_file", help="output file name")
    parser.add_argument("-r", "--ref", type=str, default='C', dest="ref", choices=['A', 'T', 'C', 'G'],
                        help="reference of the variant")
    parser.add_argument("-a", "--alt", type=str, default='T', dest="alt", choices=['A', 'T', 'C', 'G'],
                        help="alternative of the variant")
    parser.add_argument("-v", "--validate", dest="validate", default=False, action='store_true',
                        help="validate my result by samtools")
    parser.add_argument("-A", "--count-orphans", dest="count_orphans", default=False, action='store_true',
                        help='Do not skip anomalous read pairs in variant calling.')
    parser.add_argument("-Q",  "--min-BQ", dest="min_BQ", type=int, default=0,
                        help='Minimum base quality for a base to be considered' )

    args = parser.parse_args()

    abs_path_vcf = os.path.abspath(args.vcf_file)
    abs_path_bam = os.path.abspath(args.bam_file)
    abs_path_output = os.path.abspath(args.output_file)

    print("Input VCF: {}".format(abs_path_vcf))
    print("Input BAM: {}".format(abs_path_bam))
    print("Output VCF: {}".format(abs_path_output))
    print("Reference: {}".format(args.ref))
    print("Alternative: {}".format(args.alt))
    print("Validate by samtools or not: {}".format(args.validate))
    print("count-orphans: {}".format(args.count_orphans))
    print("min-BQ: {}".format(args.min_BQ))

    # vcf_file = '/diskmnt/Projects/Users/chen.xiangyu/dash/0f45d954-d951-4927-a2ba-476e319a6a88/call' \
    #            '-snp_indel_proximity_filter/execution/output/ProximityFiltered.vcf'
    # bam_file = '/diskmnt/Projects/Users/chen.xiangyu/ffpe_analysis/bam_file/CTSP-AD3X.WGS.F.hg38/104214d9-0138-4eea-8330-6df0cfca32c4_wgs_gdc_realn.bam'
    forward_backward(abs_path_vcf, abs_path_bam, abs_path_output, args.ref, args.alt, args.validate, args.count_orphans,
                     args.min_BQ)
