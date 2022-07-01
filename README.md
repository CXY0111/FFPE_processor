The code is to iterate over a VCF file. For every variant in it, search it in the bam file and calcualte the reverse strand porpotion. In order to use samtools, add environment variable first.
1.  export PATH=/diskmnt/Projects/Users/chen.xiangyu/samtools/bin:$PATH
2.  source activate ffpe
3.  python forward_backward.py -i /diskmnt/Projects/Users/chen.xiangyu/dash/0f45d954-d951-4927-a2ba-476e319a6a88/call-snp_indel_proximity_filter/execution/output/ProximityFiltered.vcf -b /diskmnt/Projects/Users/chen.xiangyu/ffpe_analysis/bam_file/CTSP-AD3X.WGS.F.hg38/104214d9-0138-4eea-8330-6df0cfca32c4_wgs_gdc_realn.bam -o output.txt
