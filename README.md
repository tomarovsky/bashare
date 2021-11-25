```assembly_stats.sh -i genome.fasta -n <num_of_chr>```
Creates len, syn, whitelist, renamelist stats files, taking into account number of chromosomes.

```bwa_index.sh -i genome.fasta```.
Creates a folder, softlinks to fasta in it, bwa index and goes back.

```bowtie2_index.sh -i genome.fasta```
Creates a folder and in it a soft reference to fasta, bowtie2-build and goes back.

```draw_coverage.sh```.
Draws genomecov coverage. 

```draw_densities_of_hetero.sh```.
Draws snps/indels coverage. 

```Krater.sh -i <fastq_prefix>```.
Fastq file prefix (file name without .final_[12].fastq). The rest of the default parameters are.

````pseudosomal_region.sh -i <file *_10000_windows_stats name> -s <name of sex scaffold>```
Creates folder pseudosomal_region, moves transferred file with stats into it, creates new one with stats only by sex chromosome, starts Biocrutch/scripts/genomecov/pseudoautosomal_region.py  

````same_filtration_but_using_bcftools.sh -i <file.vcf>```
Run in the folder prepare_region_list. Will create the same_filtration_but_using_bcftools folder and do everything up to the homo/hetero files. Will run draw_densities_of_hetero.sh.

```template.sh''.
Template for creating new scripts.
