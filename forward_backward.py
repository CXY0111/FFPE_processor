import vcf
import pysam
import re
import os
import subprocess
import argparse
from utils import get_sam_results, exact_region


def forward_backward(vcf_file, bam_file, output_file, argref, argalt, validate, fasta_ref, count_orphans, min_BQ, region):
    # generate the index file
    pysam.tabix_index(vcf_file,preset="vcf", force=True, keep_original=True)

    #  header should be:
    #       chrom, pos, ref, alt, total_read, total_forward, total_reversed, supported_forward, supported_reversed
    with open(output_file, 'w') as f:
        f.write("chrom\tpos\tref\talt\ttotal_read\ttotal_forward\ttotal_reversed\tsupported_forward\tsupported_reversed\n")

    # open the vcf file and bam file
    vcf_reader = vcf.Reader(filename=vcf_file+'.gz')
    samfile = pysam.AlignmentFile(bam_file, "rb")

    if region == 'all':                # iterate all the vcf file
        for record in vcf_reader:
            chrom = record.CHROM
            pos = record.POS
            ref = record.REF
            alt = record.ALT

            # choose only SNP
            if len(ref) != 1 or len(alt) != 1 or len(alt[0]) != 1:
                continue

            # choose variant only ref to alt
            if ref != argref or alt[0] != argalt:
                continue

            # print out the record
            print('for', record, ':')

            forward_support, reverse_support = 0, 0
            reverse_all, forward_all = 0, 0
            total = 0

            # get the aligned base in a given position
            for pileupcolumn in samfile.pileup(chrom, pos - 1, pos, min_base_quality=min_BQ,
                                               ignore_orphans=not count_orphans):
                if pileupcolumn.pos == pos - 1:
                    for pileupread in pileupcolumn.pileups:
                        # query position is None if is_del or is_refskip is set
                        if pileupread.query_position is not None:
                            base = pileupread.alignment.query_sequence[pileupread.query_position]
                            if base == alt[0]:
                                if pileupread.alignment.flag & 16:
                                    reverse_support += 1
                                else:
                                    forward_support += 1
                        if pileupread.alignment.flag & 16:
                            reverse_all += 1
                        else:
                            forward_all += 1

                        total += 1
            print('backward_suport=', reverse_support, 'forward_suport=', forward_support, 'reverse_all=', reverse_all,
                  'forward_all=', forward_all, 'total=', total)
            with open(output_file, 'a') as f:
                f.write(chrom + "\t" + str(pos) + "\t" + ref + "\t" + str(alt[0]) + "\t" + str(total) + "\t" +
                        str(forward_all) + "\t" + str(reverse_all) + "\t" + str(forward_support) + "\t" +
                        str(reverse_support) + "\n")

            # this is the part to validate the res by samtools
            if validate:
                res, sam_forward_support, sam_reverse_support, sam_forward, sam_reverse, sam_total = get_sam_results(
                    fasta_ref, bam_file, pos, alt, count_orphans, min_BQ)
                print('samtools_mpileup:')
                if sam_reverse_support == reverse_support and sam_forward_support == forward_support and total == sam_total and reverse_all == sam_reverse and forward_all == sam_forward and total == sam_total:
                    print(res, 'matches!')
                else:
                    print(res, 'not match')

            print('')
    else:                              # iterate only the selected region
        chrom, start, end = exact_region(region)
        for record in vcf_reader.fetch(chrom, start, end):
            chrom = record.CHROM
            pos = record.POS
            ref = record.REF
            alt = record.ALT

            # choose only SNP
            if len(ref) != 1 or len(alt) != 1 or len(alt[0]) != 1:
                continue

            # choose variant only ref to alt
            if ref != argref or alt[0] != argalt:
                continue

            # print out the record
            print('for', record, ':')

            forward_support, reverse_support = 0, 0
            reverse_all, forward_all = 0, 0
            total = 0

            # get the aligned base in a given position
            for pileupcolumn in samfile.pileup(chrom, pos - 1, pos, min_base_quality=min_BQ,
                                               ignore_orphans=not count_orphans):
                if pileupcolumn.pos == pos - 1:
                    for pileupread in pileupcolumn.pileups:
                        # query position is None if is_del or is_refskip is set
                        if pileupread.query_position is not None:
                            base = pileupread.alignment.query_sequence[pileupread.query_position]
                            if base == alt[0]:
                                if pileupread.alignment.flag & 16:
                                    reverse_support += 1
                                else:
                                    forward_support += 1
                        if pileupread.alignment.flag & 16:
                            reverse_all += 1
                        else:
                            forward_all += 1

                        total += 1
            print('backward_suport=', reverse_support, 'forward_suport=', forward_support, 'reverse_all=', reverse_all,
                  'forward_all=', forward_all, 'total=', total)
            with open(output_file, 'a') as f:
                f.write(chrom + "\t" + str(pos) + "\t" + ref + "\t" + str(alt[0]) + "\t" + str(total) + "\t" +
                        str(forward_all) + "\t" + str(reverse_all) + "\t" + str(forward_support) + "\t" +
                        str(reverse_support) + "\n")

            # this is the part to validate the res by samtools
            if validate:
                res, sam_forward_support, sam_reverse_support, sam_forward, sam_reverse, sam_total = get_sam_results(
                    fasta_ref, bam_file, chrom, pos, alt, count_orphans, min_BQ)
                print('samtools_mpileup:')
                if sam_reverse_support == reverse_support and sam_forward_support == forward_support and total == sam_total and reverse_all == sam_reverse and forward_all == sam_forward and total == sam_total:
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
    parser.add_argument("-f", "--fasta-ref", dest="fasta_ref", default=None,
                        help="The faidx-indexed reference file in the FASTA format.Only needed when use -v")
    parser.add_argument("-A", "--count-orphans", dest="count_orphans", default=False, action='store_true',
                        help='Do not skip anomalous read pairs in variant calling.')
    parser.add_argument("-Q", "--min-BQ", dest="min_BQ", type=int, default=0,
                        help='Minimum base quality for a base to be considered')
    parser.add_argument("-R", "--region", dest="region", type=str, default='all',
                        help='Only generate results in region. If no specify, process all the vcf file.')

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
    print("region: {}".format(args.region))

    # vcf_file = '/diskmnt/Projects/Users/chen.xiangyu/dash/0f45d954-d951-4927-a2ba-476e319a6a88/call' \
    #            '-snp_indel_proximity_filter/execution/output/ProximityFiltered.vcf'
    # bam_file = '/diskmnt/Projects/Users/chen.xiangyu/ffpe_analysis/bam_file/CTSP-AD3X.WGS.F.hg38/104214d9-0138-4eea-8330-6df0cfca32c4_wgs_gdc_realn.bam'
    # fasta_ref = '/diskmnt/Datasets_public/Reference/GRCh38.d1.vd1/GRCh38.d1.vd1.fa'
    forward_backward(abs_path_vcf, abs_path_bam, abs_path_output, args.ref, args.alt, args.validate, args.fasta_ref,
                     args.count_orphans, args.min_BQ,args.region)
