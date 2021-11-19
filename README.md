```assembly_stats.sh -i genome.fasta -n <num_of_chr>```
Создает файлы статистик len, syn, whitelist, renamelist с учетом кол-ва хромосом.

```bwa_index.sh -i genome.fasta```
Создает папку, в ней мягкую ссылку на fasta, bwa index и возвращается обратно.

```bowtie2_index.sh -i genome.fasta```
Создает папку, в ней мягкую ссылку на fasta, bowtie2-build и возвращается обратно.

```draw_coverage.sh```
Отрисует genomecov coverage. 

```draw_densities_of_hetero.sh```
Отрисует snps/indels coverage. 

```krater.sh -i <fastq_prefix>```
Fastq file prefix (file name without .final_[12].fastq). Остальные параметры по умолчанию.

```pseudosomal_region.sh -i <file *_10000_windows_stats name> -s <name of sex scaffold>```
Создаст папку pseudosomal_region, переместит в нее переданный файл со статистиками, создаст новый со статистикамии только по половой хромосоме, запустит Biocrutch/scripts/genomecov/pseudoautosomal_region.py.  

```same_filtration_but_using_bcftools.sh -i <file.vcf>```
Запускать в папке prepare_region_list. Создаст папку same_filtration_but_using_bcftools и сделает все до файлов homo/hetero. Запустит draw_densities_of_hetero.sh.

```template.sh```
Шаблон для создания новых скриптов.
