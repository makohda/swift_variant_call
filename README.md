# swift_variant_call

## Running example shell scripts (single sample analysis using bwa, gatk4, lofreq, and annovar)
**swift_amplicon.SRA5AF_sample.04.sh** is for single sample processing.
The script is designed to work with low allele frequency variants (e.g. more than 5% allele frequencies).
Before running it , YOU will specify your sample name in the script.
If you had Accel-Amplicon_TP53_acrometrix_1.fastq.gz and Accel-Amplicon_TP53_acrometrix_2.fastq.gz, you will write your sample name like this.
```
id=Accel-Amplicon_TP53_acrometrix
```
You can find this description at line 5 in swift_amplicon.SRA5AF_sample.04.sh.

Run this script, you just type
```
./swift_amplicon.SRA5AF_sample.04.sh
```

----

# Preparations, notes
Assuming to use recent macOS
```
Version information of programs are:
FastQC: v0.11.7
Trimmomatic: 0.38
bwa: 0.7.17-r1188
samtools 1.9
Using htslib 1.9
picard: 2.18.23-SNAPSHOT
The Genome Analysis Toolkit (GATK): v4.1.0.0
lofreq: 2.1.3.1
Annovar: $Date: 2018-04-16 00:47:49 -0400 (Mon, 16 Apr 2018) $
```

## Install programs using homebrew
```
$ brew install bwa
$ brew install fastqc
$ brew install picard-tools
$ brew install samtools
$ brew install lofreq
$ brew install coreutils
$ brew install grep
```

## Install Trimmomatic, IGV, GATK4, PrimerClip, Annovar
- Trimmomatic: A flexible read trimming tool for Illumina NGS data http://www.usadellab.org/cms/?page=trimmomatic

- Integrative Genomics Viewer https://software.broadinstitute.org/software/igv/download
Download IGV to run on Linux / MacOS command line
Mynote is also useful https://github.com/makohda/IDRC_first_step_course_2018#install-igv--3-min

- swiftbiosciences/primerclip: Swift Accel-Amplicon primer trimming tool for fast alignment-based primer trimming https://github.com/swiftbiosciences/primerclip

- GATK | Home https://software.broadinstitute.org/gatk/

- Annovar Download ANNOVAR - ANNOVAR Documentation http://annovar.openbioinformatics.org/en/latest/user-guide/download/


## Installing primercip into macOS is a bit tricky
Firstly, install haskell stack to build primerclip.
```
$ curl -sSL https://get.haskellstack.org/ | sh
```
Then, clone primerclip from github
```
$ git clone https://github.com/swiftbiosciences/primerclip.git
$ cd primerclip
```
Read README.md in primerclip direcotry, and run two commands for build.
```
$ stack build
$ stack install
```
But, I met following error.
=> ld: library not found for -lcrt0.o

If I understood this problem correctly, this error was caused by macOS manner.
In macOS, we can't compile softwares with -static option.
Avoiding this compile error, modify the build rule file.
```
$ emacs primerclip.cabal
```
Remove -static option from the file.
```
-O2 -static -optl-static -optl-pthread
=>
-O2 -optl-pthread
```

Now, you can compile primerclip.
```
$ stack build
$ stack install
```

## Setup ANNOVAR
I use several annotation databases in my example shell script.
- refGeneWithVer
- genomicSuperDups
- clinvar_20180603
- exac03
- gnomad_genome
- avsnp150
- ljb26_all
- cosmic87_coding
- cosmic87_noncoding
- tommo-3.5kjpnv2-20180625-af_snvall.MAF.genericdb
Of these, cosmic87 will be described in the next paragraph.
tommo-3.5kjpnv2 was generated using ToMMoVcf2Annovar https://github.com/makohda/ToMMoVcf2Annovar
Other annotations were installed using Annovar native function. See ANNOVAR User Guide http://annovar.openbioinformatics.org/en/latest/user-guide/startup/

## Download COSMIC data files
Download files from here https://cancer.sanger.ac.uk/cosmic/download?genome=37
Version was 87 at 2019/02/09.
I choosed GRCh37, because Swift masterfile is based on GRCh37.

Following four files were downloaded. These will be integrated into Annovar for annotating variants.
- VCF/CosmicCodingMuts.vcf.gz
- VCF/CosmicNonCodingVariants.vcf.gz
- CosmicMutantExport.tsv.gz
- CosmicNCV.tsv.gz

File size are relatively large, so "Scripted downloading" is prefermable (login required).
```
$ echo "registred e-mail address:***password***" | base64
=> xxxxxx

$ curl -H "Authorization: Basic xxxxxx" https://cancer.sanger.ac.uk/cosmic/file_download/GRCh37/cosmic/v87/VCF/CosmicCodingMuts.vcf.gz
=> {"url":"https://cog.sanger.ac.uk/cosmic/GRCh37/cosmic/v87/VCF/CosmicCodingMuts.vcf.gz?AWSAccessKeyId=yyyy"}

$ wget -c "https://cog.sanger.ac.uk/cosmic/GRCh37/cosmic/v87/VCF/CosmicCodingMuts.vcf.gz?AWSAccessKeyId=yyyy"" -O CosmicCodingMuts.vcf.gz
```

After completing downloads, convert these four files to annovar data format.
```
$ prepare_annovar_user.pl -dbtype cosmic CosmicMutantExport.tsv -vcf CosmicCodingMuts.vcf > hg19_cosmic87_coding.txt
$ prepare_annovar_user.pl -dbtype cosmic CosmicNCV.tsv -vcf CosmicNonCodingVariants.vcf > hg19_cosmic87_noncoding.txt
Sort by chromosomes and positions.
$ cat hg19_cosmic87_coding.txt | gsort -V -k1 > hg19_cosmic87_coding_sorted.txt
$ cat hg19_cosmic87_noncoding.txt | gsort -V -k1 > hg19_cosmic87_noncoding_sorted.txt
```
hg19_cosmic87_coding_sorted.txt and hg19_cosmic87_noncoding_sorted.txt will be located in annovar/humandb/ directory.

You can get more detail information here.
http://annovar.openbioinformatics.org/en/latest/user-guide/filter/#cosmic-annotations

## Setup GATK4
Now, we will use GATK4. But, these notes are probably useful. I wrote this at last year 2018.
https://github.com/makohda/IDRC_first_step_course_2018#install-gatk--3-min
https://github.com/makohda/IDRC_first_step_course_2018#install-gatk-bundle-resource-10-min

## Download masterfile and merged_targets from Swift web site.
Log in required.
https://swiftbiosci.com/protected-content/protected-content_amplicon-bed-files/

Accel-Amplicon Comprehensive TP53 Panel BED file and masterfile.
- Accel-Amplicon_TP53_masterfile_170228.txt
- Accel-Amplicon_TP53_merged_targets
- Accel-Amplicon_TP53_nonmerged_targets

Demonstration Data.
- Accel-Amplicon_TP53_acrometrix_IndelRealigned_BQSR.bam
- Accel-Amplicon_TP53_acrometrix_IndelRealigned_BQSR.bai
- Accel-Amplicon_TP53_acrometrix_gatkHC.vcf
- Accel-Amplicon_TP53_acrometrix_lofreq.vcf

**Accel-Amplicon_TP53_masterfile_170228.txt will be required for primerclip step.**

**Accel-Amplicon_TP53_merged_targets will be required for GATK processing. It will be used for specifying genomic regions for calcuration.**

Convert .bam to .fastq format, because it could be used as test data to test our pipeline sensitivity.
```
$ picard SamToFastq  INPUT=Accel-Amplicon_TP53_acrometrix_IndelRealigned_BQSR.bam F=Accel-Amplicon_TP53_acrometrix_1.fastq F2=Accel-Amplicon_TP53_acrometrix_2.fastq
$ gzip Accel-Amplicon_TP53_acrometrix_1.fastq
$ gzip Accel-Amplicon_TP53_acrometrix_2.fastq
```

----

## note: Reference human genome version
https://github.com/makohda/IDRC_first_step_course_2018#tips-1kg-b37-decoy

----
## Primerclip is needed to trim (softclip) Swift-generated primer sequneces from your mapped sequene files (.sam).
Sample code is like this.

Make a .sam files from fastq reads.
```
$ bwa mem -t ${thread} -M \
            -R "@RG\tID:FLOWCELLID\tSM:${id}\tPL:illumina\tLB:${id}_library_1" \
            human_g1k_v37_decoy.fasta \
            ${id}_1.fastq.gz ${id}_2.fastq.gz > ${id}.aligned_reads.sam
```
Run primerclip using masterfile. Masterfile could be downloaded from Swift web site (require registration).
```
$ primerclip Accel-Amplicon_TP53_masterfile_170228.txt ${id}.aligned_reads.sam ${id}.aligned_reads_clipped.sam
$ samtools view -@ ${thread} -1 ${id}.aligned_reads_clipped.sam | samtools sort -@4 - -o - > ${id}.aligned_reads_clipped_sorted.bam
$ samtools index -@ ${thread} ${id}.aligned_reads_clipped_sorted.bam
```
At the same time, generate .bam file from the .sam file without primerclip treatment.
```
$ samtools view -@ ${thread} -1 ${id}.aligned_reads.sam | samtools sort -@4 - -o - > ${id}.aligned_reads_sorted.bam
$ samtools index -@ ${thread} ${id}.aligned_reads_sorted.bam
```
Compare two .bam files (with, without primerclip processing).
```
$ sh IGV_2.4.18/igv.sh ${id}.aligned_reads_sorted.bam ${id}.aligned_reads_clipped_sorted.bam
```
![IGV_primerclip](https://raw.githubusercontent.com/makohda/swift_variant_call/master/images/primerclip.png)

----

# Deprecated
## Running example shell scripts#1 (single sample)
swift_amplicon.single_sample.sh is for single sample processing.
You will specify your sample name in the script.
If you had Accel-Amplicon_TP53_acrometrix_1.fastq.gz and Accel-Amplicon_TP53_acrometrix_2.fastq.gz, you will write your sample name like this.
```
id=Accel-Amplicon_TP53_acrometrix
```
You can find this description at line 13 in swift_amplicon.single_sample.sh.

Run this script, you just type
```
./swift_amplicon.single_sample.sh
```

## Running example shell scripts#2 (multiple samples)
swift_amplicon.multi_sample.sh is for multiple sample processing.
You **don't** specify your sample names in the script.

This script automatically collect your *_1.fastq.gz and generate each samle name.

So, you always set your file name as **_1.fastq.gz**

Run this script, you just type
```
./swift_amplicon.multi_sample.sh
```
**DO NOT USE** _**R**1.fastq.gz or _1.**fq**.gz in your file names.
