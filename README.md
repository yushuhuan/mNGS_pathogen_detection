# DNA-based_mNGS_pathogen_detection
For DNA-based mNGS analysis pipeline(Aliyun):

V2.0
A general analysis module was designed to obtain microbiome abundance in the sample.
Module A:
1. Fastp(0.23.4) was used for quality control, the parameter --adapter_sequence AGTCGGAGGCCAAGCGGTCTTAGGAAGACAA was used to remove adapters.
2. Bowtie2(2.3.5.1) was used to remove human sequences to improve microbiome enrichment level. The reference genome of human is hg38.
3. Kraken2(2.1.3) was used to get taxonomic information by mapping to "standard", "pluspf" and  "eupathdb" databases respectively.
4. Bracken(2.9) was used to estimate microbiome abundance on the species level.

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
