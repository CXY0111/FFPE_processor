# FFPE_processor

The code is to iterate over a VCF file. For every variant in it, search it in the bam file and calcualte the reverse 
strand porpotion.

### directly
install the libraries in requirements.txt
```
python forward_backward.py -i <location to vcf> -b <location to bam> -o <output file name>
```

### on kaimai
1. In order to use samtools, add environment variable first.
```
export PATH=/diskmnt/Projects/Users/chen.xiangyu/samtools/bin:$PATH
```

2. Activate the conda eneironment
```
source activate ffpe
```
3. Run
```
python forward_backward.py -i /diskmnt/Projects/Users/chen.xiangyu/dash/0f45d954-d951-4927-a2ba-476e319a6a88/call-snp_indel_proximity_filter/execution/output/ProximityFiltered.vcf -b /diskmnt/Projects/Users/chen.xiangyu/ffpe_analysis/bam_file/CTSP-AD3X.WGS.F.hg38/104214d9-0138-4eea-8330-6df0cfca32c4_wgs_gdc_realn.bam -o output.txt
```

### on computer1
1. git the code from github
```
git clone https://github.com/Brave-banana/FFPE_processor.git && cd FFPE_processor
```
2. Set the `LSF_DOCKER_VOLUMES` environment variable to mount the location of  files
```
export LSF_DOCKER_VOLUMES="/storage1/fs1/m.wyczalkowski/:/storage1/fs1/m.wyczalkowski/ /storage1/fs1/dinglab/Active/Projects/:/storage1/fs1/dinglab/Active/Projects/"
```
3. submit through `bsub` (Notice that this step should under directory `FFPE_processor`)
```
bsub -Is -q general-interactive -a 'docker(bravebanana/ffpe_processor)' /bin/bash
```
4. Run FFPE_processor
```
python forward_backward.py -i /storage1/fs1/m.wyczalkowski/Active/cromwell-data/cromwell-workdir/cromwell-executions/tindaisy2.ffpe.cwl/0f45d954-d951-4927-a2ba-476e319a6a88/call-snp_indel_proximity_filter/execution/output/ProximityFiltered.vcf -b /storage1/fs1/dinglab/Active/Projects/CTSP-DLBCL/primary/AWG/WGS/534d8f65-e46e-44ec-82ed-f15192b1517d/104214d9-0138-4eea-8330-6df0cfca32c4_wgs_gdc_realn.bam -o output.txt
```

### Options

***-i,--input*** STR input vcf file name    
***-b,--bam*** STR input bam file name  
***-o,--output*** STR output file name  
***-r,--ref*** STR reference of the variant, default to be 'C', choose from ['A', 'T', 'C', 'G']    
***-a,--alt*** STR alternative of the variant, default to be 'T', choose from ['A', 'T', 'C', 'G']  
***-v,--validate*** validate my result by samtools,default to be False  
***-f,--fasta-ref***STR The faidx-indexed reference file in the FASTA format.Only needed when use -v    
***-A,--count-orphans*** skip anomalous read pairs in variant calling,default to be False    
***-Q,--min-BQ*** INT Minimum base quality for a base to be considered,default to be 0  
***-R,--region*** STR Only generate results in region, default to be 'all'. 

