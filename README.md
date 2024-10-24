# Pathogen detection for DNA-based-mNGS intraocular fluid data

This is a simple framework which can be used to analysis DNA-based-mNGS intraocular fluid data. It based on mainstream analysis tools, including quality control, mapping, remove host sequences and fast multiple alignment. We designed several scripts to execute analysis.

## Turorial

### Step1: First you need to check whether following softwares has been installed on your computer

    Fastp(0.20.0)
	bwa(0.7.17)
	minirmd(V1)
	Kraken2(2.1.3)
	samtools
	seqkit
	R(4.4.1)
	Rscript

### Step2: Clone this repo

```shell
git clone https://github.com/ShuhuanYu/mNGS_pathogen_detection.git
```

### Step3: Download "nt" database

You can download indexed "nt" database by click [https://benlangmead.github.io/aws-indexes/k2](https://benlangmead.github.io/aws-indexes/k2)

![image](https://github.com/ShuhuanYu/mNGS_pathogen_detection/blob/main/images/kraken2-nt.png)

### Step4: Download "hg38" genome and index it

```shell
cd PATH_TO_YOUR_DIR/ref
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
bwa index -p hg38 GCF_000001405.40_GRCh38.p14_genomic.fna.gz
```

### Step5: Modify script parameters

For estimate_attention_index.r,

```R
concerned_pathogens_V1 <- read.csv("PATH_TO_YOUR_PATHOGEN_FILE/known_pathogen_databse_V1.0.csv", header = T, fill = TRUE)[,1:3]
```

Pathogen database is very important for mngs detection and annotation. Different organizations and different directions correspond to different reference databases. Reference databases need to be regularly updated and maintained to ensure the accuracy of pathogen detection and annotation.

For run.sh,

```shell
kraken2 --db PATH_TO_YOU_NT_DATABASE \
                --threads 8 \
                --output $sample_id.classified.nt.output \
                -report $sample_id.classified.nt.report \
                --memory-mapping \
                --use-names \
                $sample_id.unmapped.rmDup.sorted.fastq
```

### Step6:  Execute analysis

```shell
bash run.sh SAMPLE-ID NC-ID
```

## Output files

There are 3 output files in this framework.

### 1-$sample_id.classified.nt.report

This is a standard report of kraken2, including 6 columns.

```shell
90.68  446288  446288  U       0       unclassified
  9.32  45867   838     R       1       root
  9.13  44936   277     R1      131567    cellular organisms
  6.20  30492   364     D       2759        Eukaryota
  5.96  29321   7       D1      33154         Opisthokonta
  5.46  26865   0       K       33208           Metazoa
  5.46  26863   1       K1      6072              Eumetazoa
  5.46  26855   135     K2      33213               Bilateria
  5.33  26251   0       K3      33511                 Deuterostomia
  5.33  26249   0       P       7711                    Chordata
```

* column 1: reads number/sample reads, this is a relative index
* column2: number of reads aligned to this taxonomy
* column3: number of unique reads aligned to this taxonomy
* column4: domain
* column5: taxonomy id
* column6: scientific name

### 2-$sample_id.pathogen_with_genius_and_domain.txt

This file contains all "species" level pathogen detected by kraken2. Besides, species normalization and relative aboundance wa caculated by R.

```shell
Rank    tax_id  Domain  Genius  Species all_reads       CPM     percent aboundance_based_domain aboundance_based_species
1       1747    Bacteria        Cutibacterium   Cutibacterium acnes     7981    16216.44        1.62    56.44   17.4
2       5207    Eukaryota       Cryptococcus    Cryptococcus neoformans 1573    3196.15 0.32    5.16    3.43
3       180562  Eukaryota       Morina  Morina longifolia       298     605.5   0.06    0.98    0.65
4       2589401 Eukaryota       Diabelia        Diabelia stenophylla    289     587.21  0.06    0.95    0.63
5       1630339 Eukaryota       Diabelia        Diabelia stenophylla var. tetrasepala   289     587.21  0.06    0.95    0.63
6       76775   Eukaryota       Malassezia      Malassezia restricta    173     351.52  0.04    0.57    0.38
7       33011   Bacteria        Cutibacterium   Cutibacterium granulosum        162     329.16  0.03    1.15    0.35
8       40324   Bacteria        Stenotrophomonas        Stenotrophomonas maltophilia    106     215.38  0.02    0.75    0.23
9       38304   Bacteria        Corynebacterium Corynebacterium tuberculostearicum      90      182.87  0.02    0.64    0.2
```

### 3-$sample_id.pathogen_report.txt

Based on $sample_id.pathogen_with_genius_and_domain.txt and know experience, attention level was infered to assist the interpreter in making judgments. Only pathogens in ref database was included.

```shell
attention_level speciescn       tax_id  all_reads       CPM     Rank    Domain  Genius  Species percent aboundance_based_domain aboundance_based_species Reads_num_in_NC label_experience
背景微生物      痤疮丙酸杆菌    1747    7981    16216.44        1       Bacteria        Cutibacterium   Cutibacterium acnes     1.62%    56.44%  17.4%   11      Both
        新型隐球菌      5207    1573    3196.15 2       Eukaryota       Cryptococcus    Cryptococcus neoformans 0.32%   0%      3.43%    NA      Pathogenicity
        限制性马拉色菌  76775   173     351.52  6       Eukaryota       Malassezia      Malassezia restricta    0.04%   0.57%   0.38%    6       Both
        颗粒状皮肤杆菌  33011   162     329.16  7       Bacteria        Cutibacterium   Cutibacterium granulosum        0.03%   1.15%    0.35%   2       Background
        嗜麦芽窄食单胞菌        40324   106     215.38  8       Bacteria        Stenotrophomonas        Stenotrophomonas maltophilia     0.02%   0%      0.23%   NA      Both
        结核硬脂酸棒状杆菌      38304   90      182.87  9       Bacteria        Corynebacterium Corynebacterium tuberculostearicum       0.02%   0.64%   0.2%    3       Background
        球形马拉色菌    76773   42      85.34   15      Eukaryota       Malassezia      Malassezia globosa      0.01%   0%      0.09%    NA      Background
        奥斯陆莫拉菌    34062   42      85.34   16      Bacteria        Moraxella       Moraxella osloensis     0.01%   0.3%    0.09%    2       Both
        琼氏不动杆菌    40215   25      50.8    28      Bacteria        Acinetobacter   Acinetobacter junii     0.01%   0.18%   0.05%    2       Both
```
## Coverage output
Besides above 3 output files, this framework can also calculate coverage fro each pathogen, and provide personalized visualization results, such as genome coverage across reference. You can see more details from "coverage_output_demo" folder.
![image](https://github.com/yushuhuan/mNGS_pathogen_detection/blob/main/images/genome_coverage_across_reference.png)
## Contribute

Thanks to Prof. TaoYong. Thanks to Xue yong, Dr. Qian zhuyun and Wang lu. Thanks for their supporting and helping.
