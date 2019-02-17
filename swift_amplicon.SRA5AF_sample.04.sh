#!/bin/bash

set -eux # bash スクリプトを書く際のオプション。慣れてきてからググれば良い

id=Accel-Amplicon_TP53_acrometrix
thread=4 # 並列に計算出来るときは４並列で計算するためのオプション。プログラムによっては並列で走る、並列化できないものもある

# 最初にソフトウェアのバージョンを出力しておく
echo "version information"
fastqc --version | perl -pe 's/FastQC/FastQC:/'
echo "Trimmomatic: `java -Xmx4g -jar Trimmomatic-0.38/trimmomatic-0.38.jar -version`"
bwa 2>&1 | grep Version | perl -pe 's/Version/bwa/'
samtools --version | grep -v Copyright
picard MergeVcfs 2>&1 | grep Version | perl -pe 's/Version/picard/'
./gatk-4.1.0.0/gatk --version 2>&1 | grep Toolkit | perl -pe 's/ v/: v/'
lofreq version 2>&1 | grep version | perl -pe 's/version/lofreq/'
./annovar/table_annovar.pl -h 2>&1 | grep Version | perl -pe 's/ +Version/Annovar/'


date

if [ 1 -eq 0 ]; then
  :
fi

if [ ! -e ./fastqc_original ]; then
  mkdir fastqc_original
fi
fastqc -t ${thread} ${id}_1.fastq.gz ${id}_2.fastq.gz  -o fastqc_original

date

java -Xmx4g -jar Trimmomatic-0.38/trimmomatic-0.38.jar PE \
            -threads ${thread} -phred33 -trimlog ${id}.trimlog \
            ${id}_1.fastq.gz \
            ${id}_2.fastq.gz \
            ${id}_1.paired.fastq.gz ${id}_1.unpaired.fastq.gz \
            ${id}_2.paired.fastq.gz ${id}_2.unpaired.fastq.gz \
            ILLUMINACLIP:./Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10 \
            MINLEN:30
date

# fastqc_after_trimming ディレクトリが無い時だけ作る
if [ ! -e ./fastqc_after_trimming ]; then
  mkdir fastqc_after_trimming
fi
fastqc -t ${thread} ${id}_1.paired.fastq.gz ${id}_2.paired.fastq.gz  -o fastqc_after_trimming

date

bwa mem -t ${thread} -M \
            -R "@RG\tID:FLOWCELLID\tSM:${id}\tPL:illumina\tLB:${id}_library_1" \
            human_g1k_v37_decoy.fasta \
            ${id}_1.fastq.gz ${id}_2.fastq.gz > ${id}.aligned_reads.sam

samtools view -@ ${thread} -1 ${id}.aligned_reads.sam | samtools sort -@4 - -o - > ${id}.aligned_reads_sorted.bam
# indexをはる
samtools index -@ ${thread} ${id}.aligned_reads_sorted.bam

# primerclip をかける
primerclip Accel-Amplicon_TP53_masterfile_170228.txt ${id}.aligned_reads.sam ${id}.aligned_reads_clipped.sam
# sortして、samからbamに変換
samtools view -@ ${thread} -1 ${id}.aligned_reads_clipped.sam | samtools sort -@4 - -o - > ${id}.aligned_reads_clipped_sorted.bam
# indexをはる
samtools index -@ ${thread} ${id}.aligned_reads_clipped_sorted.bam

# primerclip前の.samも比較用に.bamにしておく
samtools view -@ ${thread} -1 ${id}.aligned_reads.sam | samtools sort -@4 - -o - > ${id}.aligned_reads_sorted.bam
samtools index -@ ${thread} ${id}.aligned_reads_sorted.bam


date

./gatk-4.1.0.0/gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=8" BaseRecalibrator \
                    -R human_g1k_v37_decoy.fasta \
                    --known-sites dbsnp_138.b37.vcf \
                    --known-sites Mills_and_1000G_gold_standard.indels.b37.vcf \
                    -I ${id}.aligned_reads_clipped_sorted.bam \
                    -O ${id}_recal.table

date

./gatk-4.1.0.0/gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=8" ApplyBQSR \
                    -R human_g1k_v37_decoy.fasta \
                    -I ${id}.aligned_reads_clipped_sorted.bam \
                    -bqsr ${id}_recal.table \
                    -O ${id}.aligned_reads_clipped_recal_sorted.bam

samtools index ${id}.aligned_reads_clipped_recal_sorted.bam

date

# ここから lofreq + Annovar

lofreq indelqual --dindel -f human_g1k_v37_decoy.fasta ${id}.aligned_reads_clipped_recal_sorted.bam -o ${id}.aligned_reads_clipped_recal_qindel_sorted.bam
samtools index ${id}.aligned_reads_clipped_recal_qindel_sorted.bam
# TP53のみ
#lofreq call --force-overwrite -f human_g1k_v37_decoy.fasta -o lofreq_vars.vcf -r 17:7570000-7580000 --call-indels ${id}.aligned_reads_clipped_recal_qindel_sorted.bam
# 全ゲノム
#lofreq call --force-overwrite -f human_g1k_v37_decoy.fasta -o lofreq_vars.${id}.vcf --call-indels ${id}.aligned_reads_clipped_recal_qindel_sorted.bam
##source=lofreq call -q 20 -Q 20 -m 30 -C 50 --call-indels -f /seq/refGenomes/Homo_sapiens_assembly19broad.fasta -l merged_targets_5col.bed -o tp53-acro-1_S1.realigned.bam_bqsr.lf.vcf tp53-acro-1_S1.realign
lofreq call --force-overwrite -q 20 -Q 20 -m 30 -C 50 --call-indels -f human_g1k_v37_decoy.fasta -l tp53_170228_merged_targets.bed -o lofreq_vars.${id}.vcf ${id}.aligned_reads_clipped_recal_qindel_sorted.bam

./gatk-4.1.0.0/gatk --java-options "-Xmx4G" VariantFiltration \
                    -R human_g1k_v37_decoy.fasta \
                    -V lofreq_vars.${id}.vcf \
                    --filter-expression 'AF < 0.01' --filter-name 'ExtremeLowAlleleFrequency' \
                    -O lofreq_vars.${id}.filtered.vcf

./annovar/convert2annovar.pl -format vcf4 \
                             --includeinfo \
                             --withzyg \
                             --allsample \
                             lofreq_vars.${id}.filtered.vcf \
                             --outfile lofreq_vars.${id}.filtered.avinput

date

./annovar/table_annovar.pl lofreq_vars.${id}.filtered.avinput annovar/humandb/ \
                           -buildver hg19 \
                           -protocol refGeneWithVer,genomicSuperDups,cosmic87_coding_sorted,cosmic87_noncoding_sorted,clinvar_20180603,generic,exac03,gnomad_genome,avsnp150,ljb26_all \
                           -genericdb tommo-3.5kjpnv2-20180625-af_snvall.MAF.genericdb \
                           -operation g,r,f,f,f,f,f,f,f,f \
                           -nastring NA \
                           --otherinfo \
                           --argument '--time --chromosome 17 --hgvs --exonicsplicing --splicing_threshold 2',,,,,,,,, \
                           --remove \
                           --polish \
                           --thread ${thread} \
                           -out lofreq_vars.${id}.filtered.avoutput

date

cat lofreq_vars.${id}.filtered.avoutput.hg19_multianno.txt | grep -v "AF=0\.00" > lofreq_vars.${id}.filtered.avoutput.hg19_multianno.rm_extreme_low.txt
cat lofreq_vars.${id}.filtered.avoutput.hg19_multianno.rm_extreme_low.txt | grep -e Chr -e exonic -e splcing > lofreq_vars.${id}.filtered.avoutput.hg19_multianno.rm_extreme_low.exonic_splicing.txt
cat lofreq_vars.${id}.filtered.avoutput.hg19_multianno.rm_extreme_low.exonic_splicing.txt | csvlook -t | less -S


# file move/remove
if [ ! -e ./fastq_paired_unpaired ]; then
  mkdir fastq_paired_unpaired
fi
mv *paired.fastq.gz fastq_paired_unpaired/

if [ ! -e ./bam_dir ]; then
  mkdir bam_dir
fi
mv *.bam *.bai bam_dir/
rm *.sam

if [ ! -e ./trimlog_dir ]; then
  mkdir trimlog_dir
fi
mv *.trimlog trimlog_dir/

if [ ! -e ./gvcf_dir ]; then
  mkdir gvcf_dir
fi
mv *.g.vcf *.g.vcf.idx gvcf_dir/

if [ ! -e ./vcf_dir ]; then
  mkdir vcf_dir
fi
mv combined_genotyped*.vcf combined_genotyped*.vcf.idx vcf_dir/

if [ ! -e ./log_dir ]; then
  mkdir log_dir
fi
mv *primerclip_runstats.log *.trimlog *_recal.table masterparsefails.log log_dir/
echo "Finish!"

exit
