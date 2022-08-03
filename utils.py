import subprocess
import re


def get_sam_results(fasta_ref, bam_file, chrom, pos, alt, count_orphans, min_BQ, min_MQ, max_depth,ignore_overlaps):
    """
    this is the part to validate the res by samtools

    :param fasta_ref: The faidx-indexed reference file in the FASTA format.
    :param bam_file: Path to bam file
    :param pos: position to be validated
    :param alt: alternatives
    :param count_orphans: whether to skip anomalous read pairs in variant calling
    :param min_BQ: Minimum base quality for a base to be considered
    :return: res, forward_support, reverse_support, forward, reverse, total
    :rtype: multiple
    """
    # if not count_orphans:
    #     cmd = 'samtools mpileup -f ' + fasta_ref + ' -s ' + bam_file + ' -r ' + chrom + ':' + str(pos) + '-' + \
    #           str(pos) + ' -Q ' + str(min_BQ) + ' --reverse-del -B | cut -f 1,2,3,4,5 '
    # else:
    #     cmd = 'samtools mpileup -f ' + fasta_ref + ' -s ' + bam_file + ' -r ' + chrom + ':' + str(pos) + '-' + \
    #           str(pos) + ' -Q ' + str(min_BQ) + ' -A --reverse-del -B | cut -f 1,2,3,4,5 '
    # # print(cmd)
    cmd = '/usr/src/samtools/bin/samtools mpileup -f ' + fasta_ref + ' -s ' + bam_file + ' -r ' + chrom + ':' + str(pos) + '-' + str(pos) + ' -Q ' + str(min_BQ) + ' --reverse-del -B -d ' + str(max_depth) + ' -q ' + str(min_MQ)
    if count_orphans:
        cmd += ' -A'
    if not ignore_overlaps:
        cmd += ' -x'
    cmd += ' | cut -f 1,2,3,4,5 '
    # print(cmd)
    ps = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output = ps.communicate()[0]
    res = output.decode('ascii').split('\n')[1]
    result = res.split('\t')[4]

    forward_support, reverse_support = 0, 0
    forward, reverse = 0, 0
    total = 0

    # get rid of begining and ending character
    result = re.sub("\^.|\$", '', result)

    # get rid of inserts and deletes
    new_res = ''
    i = 0
    while i < len(result):
        if result[i] == '+' or result[i] == '-':
            num = re.match('[0-9]+', result[i + 1:]).group(0)
            i += int(num) + len(num)
        else:
            new_res += result[i]
        i += 1
    # # get rid of begining and ending character
    # new_res = re.sub("\^.|\$", '', new_res)

    for s in new_res:

        if s == str(alt[0]):
            forward_support += 1
        elif s == str(alt[0]).lower():
            reverse_support += 1

        if s == ',' or 97 <= ord(s) <= 122 or s == '#':
            reverse += 1
        elif s == '.' or 65 <= ord(s) <= 90 or s == '*':
            forward += 1
        # elif s == ',' or 97 <= ord(s) <= 122 or s == '.' or 65 <= ord(s) <= 90 or s == '*':
        #     other += 1
        total += 1
    # print('backward_suport=', reverse_support, 'forward_suport=', forward_support, 'reverse_all=', reverse,
    #       'forward_all=', forward, 'total=', total)
    return res, forward_support, reverse_support, forward, reverse, total


def exact_region(region):
    """
    convert str region to chrom, start, end
    :param region: STR, usr input
    :return: chrom, start, end
    :rtype: multiple
    """
    chrom = re.findall('.+:', region)[0][:-1]
    start = int(re.findall(':.+-', region)[0][1:-1])
    end = int(re.findall('-.+', region)[0][1:])
    return chrom, start, end
