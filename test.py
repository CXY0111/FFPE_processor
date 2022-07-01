import vcf
import pysam
import re
import os
import subprocess
import argparse


def forward_backward(vcf_file, bam_file, output, argref, argalt, validate):
    # two cases:
    # 1. if no validation by samtools, header should be:
    #    chrom, pos, ref, alt, total_read, supported_forward, supported_reversed
    # 2. if validate by samtools, header should be:
    #    chrom, pos, ref, alt, total_read, supported_forward, supported_reversed, sam_result
    with open(output, 'w') as f:
        if validate:
            f.write("chrom\tpos\tref\talt\ttotal_read\tsupported_forward\tsupported_reversed\tsam_result\n")
        else:
            f.write("chrom\tpos\tref\talt\ttotal_read\tsupported_forward\tsupported_reversed\n")

    i = 0

    # open the vcf file and bam file, open the output file
    vcf_reader = vcf.Reader(open(vcf_file, 'r'))
    samfile = pysam.AlignmentFile(bam_file, "rb")
    w = open(output, 'a')

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

        print('for', record, ':')

        forward_support, reverse_support = 0, 0
        total = 0

        # get the aligned base in a given position
        for pileupcolumn in samfile.pileup(chrom, pos - 1, pos, min_base_quality=0, ignore_orphans=False):
            if pileupcolumn.pos == pos - 1:
                for pileupread in pileupcolumn.pileups:
                    # query position is None if is_del or is_refskip is set
                    if pileupread.query_position is None:
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
        w.write(chrom + "\t" + str(pos) + "\t" + ref + "\t" + str(alt[0]) + "\t" + str(total) + "\t" + str(
            forward_support) + "\t" + str(reverse_support) + "\n")

        # this is the part to validate the res by samtools
        cmd = 'samtools mpileup -f /diskmnt/Datasets_public/Reference/GRCh38.d1.vd1/GRCh38.d1.vd1.fa -s ' + bam_file + \
              ' -r chr1:' + str(pos) + '-' + str(pos) + ' -Q 0 -A|cut -f 1,2,3,4,5 '
        ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output = ps.communicate()[0]
        res = output.decode('ascii').split('\n')[1]

        forward, backward, other = 0, 0, 0

        # get rid of inserts and deletes
        new_res = re.sub("[+|-][0-9]+[ACGTNacgtn*#]+", '', res.split('\t')[4])
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
        #
        # # print 10 records
        # i += 1
        # if i > 10:
        #     break


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

    # vcf_file = '/diskmnt/Projects/Users/chen.xiangyu/dash/0f45d954-d951-4927-a2ba-476e319a6a88/call' \
    #            '-snp_indel_proximity_filter/execution/output/ProximityFiltered.vcf'
    # bam_file = '/diskmnt/Projects/Users/chen.xiangyu/ffpe_analysis/bam_file/CTSP-AD3X.WGS.F.hg38/104214d9-0138-4eea-8330-6df0cfca32c4_wgs_gdc_realn.bam'
    forward_backward(abs_path_vcf, abs_path_bam, abs_path_output, args.ref, args.alt, args.validate)
