# Pathogen detection for DNA-based-mNGS intraocular fluid data

This is a simple framework which can be used to analysis DNA-based-mNGS intraocular fluid data. It based on mainstream analysis tools, including qualiity control, mapping, remove host sequences and fast multiple alignment. We designed several scripts to execute analysis.

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

### Step3:Download "nt" database

You can download indexed "nt" database by click [https://benlangmead.github.io/aws-indexes/k2](https://benlangmead.github.io/aws-indexes/k2)


Download "hg38" genome and index it

Execute analysis

* Fastp(0.23.4) was used for quality control, the parameter --adapter_sequence AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA was used to remove adapters.
* Bowtie2(2.3.5.1) was used to remove human sequences to improve microbiome enrichment level. The reference genome of human is hg38.
* Kraken2(2.1.3) was used to get taxonomic information by mapping to "nt" database.
* samtools
* seqkit
* minirmd
* R
* Rscript

All taxonomic abundance results then be combined to get overall results by Krakentools(1.2).

To get eyes-related pathogens, a filter step was designed by using filter_bracken.out.py in Krakentools(1.2), setting the parameters --include and --exclude could get interested taxonomic results.

Filtered abundance tables were combined by combine_kreports.py in Krakentools(1.2). Then use ImportTaxonomy.pl in Krona(2.8.1) tranfer bracken results to krona format(html). It can be viewed visually in the browser.

Docker record:

Softwares version:
  fastp 0.20.1
  bowtie2 2.4.4
  kraken2 2.1.2
  bracken 2.9
  krona 2.8.1
  krakentools 1.2

########################################
V2.1:

1. Fastp(0.23.4) was used for quality control, the parameter --adapter_sequence AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA was used to remove adapters.
2. bwa(2.3.5.1) was used to remove human sequences to improve microbiome enrichment level. The reference genome of human is hg38.
3. Kraken2(2.1.3) was used to get taxonomic information by mapping to "blast/nt" database.
4. R scripts and shell scripts were used to get more specific information.
