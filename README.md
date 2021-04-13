# NIPT by cfDNA screening - NUS SoC Summer Workshop 2021
Non-Invasive Prenatal Testing by cell-free DNA screening

National University of Singapore - School of Computing

Summer Workshop 2021 


# Introduction

The genome can be seen as a very long string composed from an alphabet of 4 "letters": {A,C,G,T,(sometimes N for unknown)}. The letters are called **nucleotides or bases**.
**Next-generation sequencing (NGS)** technologies such as [Illumina](https://www.youtube.com/watch?v=fCd6B5HRaZ8) are able to recognize the nucleotides on short "substrings" or fragments called **reads** of length 100-300 nucleotides. Sequencing a sample generates millions of such substrings/reads so that genome is tiled to about 0.1 to 1 times in the dataset. We say it is an **ultra-low pass whole genome sequencing (ULP-WGS)** with a 0.1x to 1x depth of coverage.
For example, the picture below is a part of genome for a deeper depth of coverage (~30x). You can see its sequence at the bottom. The grey bars in the middle are reads.

![68747470733a2f2f692e696d6775722e636f6d2f393559733372452e706e67](https://user-images.githubusercontent.com/30068394/114520681-456f6600-9c74-11eb-8898-f4893cd96b93.png)


# Part 0 - Setup

The following packages are needed.
The following presumes $HOME/bin is in your $PATH. (Q: do I need to specify how to do this?)

- [bwa](https://sourceforge.net/projects/bio-bwa/files/): for aligning reads to the reference genome
```
git clone https://github.com/lh3/bwa.git
make -C bwa
cp bwa/bwa $HOME/bin
```
- [samtools](https://sourceforge.net/projects/samtools/files/): for aligned file manipulation
```
git clone https://github.com/samtools/samtools.git
make -C samtools
cp samtools/samtools $HOME/bin
```
- [ivg](http://software.broadinstitute.org/software/igv/download): for visualisation
It requires X11 to show a window on your host, if you run it on a server without monitor. For windows user, it can be easily solved by using MobaXterm as terminal.

- [wget](): to download files from the Internet using the command line
If wget is not yet installed, you may want to run the following command in Linux.
```
sudo apt-get install wget
```

# Part 1 - Upstream analysis

## The reference human genome

The human genome is composed of 3 trillions of nucleotides or base pairs and organised into 22+XY chromosomes. Each chromosome is an independent long string.
Here is the karyotype (picture overview) of a human male genome.

![DNA_human_male_chromosomes](https://user-images.githubusercontent.com/30068394/114515152-a136f080-9c6e-11eb-9cfa-0260f121cf67.png)

The human reference genome file has first been elaborated from 2000 to 2003 by the [Human Genome Project](https://www.genome.gov/human-genome-project), a long-term international research program. 

### 1. Download the human genome on your computer.
```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
```
Unzip the downloaded file.
```
gzip -d hg38.fa.gz
```

### 2. Reference indexes

To allow tools to quickly access certain regions in the reference, we need to generate index files. We will generate here both indexes from samtools and BWA.

#### A. FASTA file index from samtools

The .fai FASTA index format records the lengths of the various sequences in the reference and their offsets from the beginning of the file. Generate this index by running the following:
```
samtools faidx hg38.fa
```
You can now quickly obtain any part of the human reference genome. For instance, we can get the 200bp in chromosome 1 from position 1000000 to 1000200 by using a special format to describe the target region [chr]:[start]-[end].
```
samtools faidx hg38.fa chr1:1000000-1000200
```
You should get the following output describing the region:
```
>chr1:1000000-1000200
GGTGGAGCGCGCCGCCACGGACCACGGGCGGGCTGGCGGGCGAGCGGCGAGCGCGCGGCG
ATCCGAGCCCCTAGGGCGGATCCCGGCTCCAGGCCCGCGCGCGCCTCAGGCCGTTTCCCT
ATTTAAGGCCTCGCCGCCGCGGGGTGTGTGAACCCGGCTCCGCATTCTTTCCCACACTCG
CCCCAGCCAATCGACGGCCGC
```

#### B. BWA's FM index

BWA uses the FM-index, which a compressed full-text substring index based around the Burrows-Wheeler transform [XXX]. Generate this index by running the following:
``` 
bwa index hg38.fa
```
You should see bwa generate some information about the build process:
```
[bwa_index] Pack FASTA... 18.88 sec
[bwa_index] Construct BWT for the packed sequence...
[BWTIncCreate] textLength=6434693834, availableWord=464768632
[BWTIncConstructFromPacked] 10 iterations done. 99999994 characters processed.
[BWTIncConstructFromPacked] 20 iterations done. 199999994 characters processed.
[BWTIncConstructFromPacked] 30 iterations done. 299999994 characters processed.
...
[BWTIncConstructFromPacked] 700 iterations done. 6411146922 characters processed.
[BWTIncConstructFromPacked] 710 iterations done. 6432978554 characters processed.
[bwt_gen] Finished constructing BWT in 711 iterations.
[bwa_index] 2257.81 seconds elapse.
[bwa_index] Update BWT... 26.14 sec
[bwa_index] Pack forward-only FASTA... 11.69 sec
[bwa_index] Construct SA from BWT and Occ... 1002.41 sec
[main] Version: 0.7.17-r1188
[main] CMD: bwa index hg38.fa
[main] Real time: 3324.306 sec; CPU: 3316.926 sec
```

Finally, let's quickly check that the new index files have been created using the FASTA file name as prefix:
```
ls -lh hg38.fa*
```
You should get something like the following:
```
3.0G hg38.fa
20K  hg38.fa.amb
438K hg38.fa.ann
3.0G hg38.fa.bwt
151K hg38.fa.fai
767M hg38.fa.pac
1.5G hg38.fa.sa
```

## A liquid biposy sample

Let's now download a sequenced liquid biopsy sample from a pregnant mother for NIPT.
(Q: Do we take a sample from the igene company? If yes, where do we store the FASTQ file? On GitHub? -> privacy issue + large file Or do we take a sample from NCBI?)
The raw output of the Illumina NGS sequencing is stored in a FASTQ [XXX] file.

We can have a look at the first read of the file:
```
head -4 Ion.fastq
```
You can obtain the following. In a FASTQ, each read is stored into 4 lines with:
1. the read ID
2. the raw sequence of letters/nucleotides
3. 2. a line with the '+' symbol
4. the base quality score [XXX] for each base in line 2 in ASCII.

```
@R1FMJ:04530:11901
CCTAACCTAACCCTAACCCTAACCTACCTA
+
;48808/58/55*5:/55/5:6:/57/)//
```

### Read mapping

#### 1. Align the fastq file against the reference human genome.
Replace $threads by the number of CPUs you want to use.
This generates the aligned BAM file. XXX
```
bwa mem -t $threads hg38.fa IonXpress_001_18B0246378_1_2019-02-21.fastq
```

#### 2. Sort the BAM file according to the genome coordinates.
```
samtools sort -@ $threads -o IonXpress_001_18B0246378_1_2019-02-21.sorted.bam IonXpress_001_18B0246378_1_2019-02-21.bam
```

#### 3. Index the sorted BAM file.
Many analysis tools will require the BAM file to be indexed to be able to quickly access reads mapped to specific genome regions.
```
samtools index IonXpress_001_18B0246378_1_2019-02-21.sorted.bam
```
Finally, you can quickly look at the header of you aligned sample. You can check you file is sorted by looking at the first line. Then, the header gives you the list of chromosome or reference names. 
```
samtools view -H IonXpress_001_18B0246378_1_2019-02-21.sorted.bam
```
You can also print the 5 first reads.
```
samtools view IonXpress_001_18B0246378_1_2019-02-21.sorted.bam | head -5
```
### Summary statistics

You can use samtools to quickly get some statistics about your sample. (Q: Is this interesting?)

### Visualisation

IGV is a popular API to visualise genetic files like FASTA or BAM.
After IGV starts, you can see a window.

Then,
1. Load genome fasta file

![68747470733a2f2f692e696d6775722e636f6d2f7a746d34534f422e706e67](https://user-images.githubusercontent.com/30068394/114517489-f8d65b80-9c70-11eb-8ee6-95f793c35da3.png)

2. Load the bam

![68747470733a2f2f692e696d6775722e636f6d2f394f39525a737a2e706e67](https://user-images.githubusercontent.com/30068394/114517537-0390f080-9c71-11eb-9b0f-315d35a5c851.png)

3. Play with it

![68747470733a2f2f692e696d6775722e636f6d2f636651796c36572e706e67](https://user-images.githubusercontent.com/30068394/114520743-50c29180-9c74-11eb-871c-96d9e09444b0.png)



# Part 2 - Downstream analysis

To extract read information from BAM files in coding manner, an easy way is to learn Python and use the [Pysam](https://pysam.readthedocs.io/en/latest/) package.
(Q: what do I put there? Shall I stop the tutorial there?)
