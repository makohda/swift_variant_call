#!/bin/bash

## Usage
## give a permission
# chmod +x test_run.gatk4.basic.sh
## run with message logging
# ./test_run.gatk4.basic.sh 2>&1 | tee test_run.gatk4.basic.log
## When you want to stop running, Ctr+C

# 使わないところは # でコメントアウトできて、プログラムとしては無視される
set -eux # bash スクリプトを書く際のオプション。慣れてきてからググれば良い
thread=4 # 並列に計算出来るときは４並列で計算するためのオプション。プログラムによっては並列で走る、並列化できないものもある

# 最初にソフトウェアのバージョンを出力しておく
echo "version information"
fastqc --version | perl -pe 's/FastQC/FastQC:/'
echo "Trimmomatic: `java -Xmx4g -jar Trimmomatic-0.38/trimmomatic-0.38.jar -version`"
bwa 2>&1 | grep Version | perl -pe 's/Version/bwa/'
samtools --version | grep -v Copyright
picard MergeVcfs 2>&1 | grep Version | perl -pe 's/Version/picard/'
./gatk-4.1.0.0/gatk --version 2>&1 | grep Toolkit | perl -pe 's/ v/: v/'
./annovar/table_annovar.pl -h 2>&1 | grep Version | perl -pe 's/ +Version/Annovar/'
##

date

if [ 1 -eq 0 ]; then
  :
fi

for id in *_1.fastq.gz; do
  id=`echo ${id} | perl -pe 's/_1.fastq.gz//'`
  echo $id
  # Sequence base quality 見たい時はここで fastq ファイルに fastQC かける
  # fastqc_original ディレクトリが無い時だけ作る
  if [ ! -e ./fastqc_original ]; then
    mkdir fastqc_original
  fi
  fastqc -t ${thread} ${id}_1.fastq.gz ${id}_2.fastq.gz  -o fastqc_original

  date

  # quality trimming (exome や whole genome だと極端な話しなくても良いが、それ以外のアプリケーションでは必ずすべき)
  java -Xmx4g -jar Trimmomatic-0.38/trimmomatic-0.38.jar PE \
              -threads ${thread} -phred33 -trimlog ${id}.trimlog \
              ${id}_1.fastq.gz \
              ${id}_2.fastq.gz \
              ${id}_1.paired.fastq.gz ${id}_1.unpaired.fastq.gz \
              ${id}_2.paired.fastq.gz ${id}_2.unpaired.fastq.gz \
              ILLUMINACLIP:./Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10 \
              MINLEN:30
  #           TRAILING:20   \
  #           trimmer=["ILLUMINACLIP:" + config["illuminaclip_file"] + ":2:30:10", "MINLEN:30"]
  # 設定の参考はここから https://github.com/clinical-genomics-uppsala/accel_amplicon_trimming/blob/master/rules/accel_amplicon.smk
  # 2:30:10は本家のデフォルト設定っぽいが、もう少し詳細はこちらから Trimmomaticのilluminaclipについて - kuroの覚え書き http://k-kuro.hatenadiary.jp/entry/20170829/p1
  date

  # fastqc_after_trimming ディレクトリが無い時だけ作る
  if [ ! -e ./fastqc_after_trimming ]; then
    mkdir fastqc_after_trimming
  fi
  fastqc -t ${thread} ${id}_1.paired.fastq.gz ${id}_2.paired.fastq.gz  -o fastqc_after_trimming

  date

  # ヒトゲノムへのマッピング
  # 参照にするヒトゲノム配列は human_g1k_v37_decoy
  # これは GRCh37をベースにして 1000 genome project の中でエクソーム用にカスタムされたもの
  # 1-22, X, Y はそのまま GRCh37 なので、このまま使っても、GRCh37や38に切り替えても問題はないです
  # ID:FLOWCELLID は本来なら実験時のフローセル番号を入れると良いが、手作業実験の場合は大変なので入れていない
  bwa mem -t ${thread} -M \
              -R "@RG\tID:FLOWCELLID\tSM:${id}\tPL:illumina\tLB:${id}_library_1" \
              human_g1k_v37_decoy.fasta \
              ${id}_1.fastq.gz ${id}_2.fastq.gz > ${id}.aligned_reads.sam

  # primerclip をかける
  primerclip Accel-Amplicon_TP53_masterfile_170228.txt ${id}.aligned_reads.sam ${id}.aligned_reads_clipped.sam
  # sortして、samからbamに変換
  samtools view -@4 -1 ${id}.aligned_reads_clipped.sam | samtools sort -@4 - -o - > ${id}.aligned_reads_clipped_sorted.bam
  # indexをはる
  samtools index -@ ${thread} ${id}.aligned_reads_clipped_sorted.bam

  # primerclip前の.samも比較用に.bamにしておく
  samtools view -@4 -1 ${id}.aligned_reads.sam | samtools sort -@4 - -o - > ${id}.aligned_reads_sorted.bam
  samtools index -@ ${thread} ${id}.aligned_reads_sorted.bam

  sleep 10
  echo "*******************************************************************************************"
  echo "sh IGV_2.4.18/igv.sh ${id}.aligned_reads_sorted.bam ${id}.aligned_reads_clipped_sorted.bam"
  echo "上のコマンドで、IGV が primerclip の前後の .bam ファイルを開いてくれます"
  echo "*******************************************************************************************"
  sleep 10

  date

  ## dedup は amplicon-sequence なのでしない

  # GATK のメモ、たぶんググった方が早いので、特に時間が余ったときだけ眺めれば良い
  # GATK4 User Guide — GATK-Forum https://gatkforums.broadinstitute.org/gatk/categories/gatk4-user-guide
  # Best Practices Workflows — GATK-Forum https://gatkforums.broadinstitute.org/gatk/categories/best-practices-workflows
  ## ここからそれぞれの best practice に辿れるし、github へも辿れる
  ## Germline short variant discovery (SNPs + Indels) — GATK-Forum https://gatkforums.broadinstitute.org/gatk/discussion/11145/germline-short-variant-discovery-snps-indels
  # Tutorials — GATK-Forum https://gatkforums.broadinstitute.org/gatk/categories/gatk4-tutorials
  ## GATK | Doc #11813 | (How to) Consolidate GVCFs for joint calling with GenotypeGVCFs https://software.broadinstitute.org/gatk/documentation/article?id=11813
  ## GATK | Doc #11090 | (How to) Run GATK in a Docker container https://software.broadinstitute.org/gatk/documentation/article?id=11090
  # Dictionary — GATK-Forum https://gatkforums.broadinstitute.org/gatk/categories/gatk4-dictionary
  ## Read groups — GATK-Forum https://gatkforums.broadinstitute.org/gatk/discussion/11015/read-groups
  ## GenomicsDB — GATK-Forum https://gatkforums.broadinstitute.org/gatk/discussion/11091/genomicsdb
  # Frequently Asked Questions — GATK-Forum https://gatkforums.broadinstitute.org/gatk/categories/gatk4-faqs
  # Solutions to Problems — GATK-Forum https://gatkforums.broadinstitute.org/gatk/categories/gatk4-solutions

  #./gatk-4.1.0.0/gatk BaseRecalibrator \
  ./gatk-4.1.0.0/gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=8" BaseRecalibrator \
                      -R human_g1k_v37_decoy.fasta \
                      --known-sites dbsnp_138.b37.vcf \
                      --known-sites Mills_and_1000G_gold_standard.indels.b37.vcf \
                      -I ${id}.aligned_reads_clipped_sorted.bam \
                      -O ${id}_recal.table

  # UseParallelGCについては MarkDuplicatesSparkが2-4, BaseRecalibratorが~24とか主張している
  # https://www.biorxiv.org/content/biorxiv/early/2018/06/18/348565.full.pdf
  date

  ./gatk-4.1.0.0/gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=8" ApplyBQSR \
                      -R human_g1k_v37_decoy.fasta \
                      -I ${id}.aligned_reads_clipped_sorted.bam \
                      -bqsr ${id}_recal.table \
                      -O ${id}.aligned_reads_clipped_recal_sorted.bam

  samtools index ${id}.aligned_reads_clipped_sorted.bam

  date

  # ここらへんから -L tp53_170228_merged_targets.bed というオプションを使っている
  # これは **.bed にかかれたゲノム領域だけで計算してくれ、という範囲指定オプション
  # これをしないと全ゲノムで計算されるので時間がかかるため、計算時間節約のために使っている
  # 用いた実験キットのパネルの種類により、与えるファイルが変わることに注意
  # -L をそもそも削除すると全ゲノム、-L 17 とかくと17番染色体のみになる
  ./gatk-4.1.0.0/gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=8" HaplotypeCaller \
                      -R human_g1k_v37_decoy.fasta \
                      -I ${id}.aligned_reads_clipped_recal_sorted.bam \
                      --dbsnp dbsnp_138.b37.vcf \
                      -ERC GVCF \
                      -A AlleleFraction \
                      -A DepthPerAlleleBySample \
                      -A DepthPerSampleHC \
                      -G AS_StandardAnnotation \
                      -G StandardAnnotation \
                      -G StandardHCAnnotation \
                      -L tp53_170228_merged_targets.bed \
                      -O ${id}_raw_variants.g.vcf
  # HaplotypeCaller https://software.broadinstitute.org/gatk/documentation/tooldocs/4.1.0.0/org_broadinstitute_hellbender_tools_walkers_haplotypecaller_HaplotypeCaller.php

done

# GenomicDBImport のための .sample_map ファイルをまず作ります
ls -1 *_raw_variants.g.vcf | perl -ne 's/(.*)_raw_variants.g.vcf/$1/; print "$1\t$1_raw_variants.g.vcf\n"' > gatk.sample_map
cat gatk.sample_map

# GenomicDBImport start
./gatk-4.1.0.0/gatk --java-options "-Xmx4G -Xms4G" GenomicsDBImport \
                    --batch-size 50 \
                    -L tp53_170228_merged_targets.bed \
                    --sample-name-map gatk.sample_map \
                    --reader-threads 5 \
                    --genomicsdb-workspace-path genomicsdb
# GenomicsDBImport https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.7.0/org_broadinstitute_hellbender_tools_genomicsdb_GenomicsDBImport.php
# GATK | Doc #11009 | Intervals and interval lists https://software.broadinstitute.org/gatk/documentation/article?id=11009

./gatk-4.1.0.0/gatk --java-options "-Xmx4G" GenotypeGVCFs \
                    -R human_g1k_v37_decoy.fasta \
                    -V gendb://genomicsdb \
                    -O combined_genotyped.vcf
# GenotypeGVCFs https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.7.0/org_broadinstitute_hellbender_tools_walkers_GenotypeGVCFs.php

date

# ここで出てくる combined_genotyped_raw_snps.vcf が一塩基変化で何もフィルターをかけていない結果になる
./gatk-4.1.0.0/gatk --java-options "-Xmx4G" SelectVariants \
                    -R human_g1k_v37_decoy.fasta \
                    -V combined_genotyped.vcf \
                    --select-type-to-include SNP \
                    -O combined_genotyped_raw_snps.vcf
# SelectVariants https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.7.0/org_broadinstitute_hellbender_tools_walkers_variantutils_SelectVariants.php

date

./gatk-4.1.0.0/gatk --java-options "-Xmx4G" VariantFiltration \
                    -R human_g1k_v37_decoy.fasta \
                    -V combined_genotyped_raw_snps.vcf \
                    --cluster-size 3 --cluster-window-size 10 \
                    --filter-expression 'QD < 2.0'                               --filter-name 'LowQD' \
                    --filter-expression 'FS > 60.0'                              --filter-name 'HighFisherStrand' \
                    --filter-expression 'HaplotypeScore > 13.0'                  --filter-name 'HighHaplotypeScore' \
                    --filter-expression 'MQ < 40.0'                              --filter-name 'lowRMSMappingQuality' \
                    --filter-expression 'MQRankSum < -12.5'                      --filter-name 'LowMQRankSum' \
                    --filter-expression 'ReadPosRankSum < -8.0'                  --filter-name 'LowReadPosRankSum' \
                    --filter-expression 'MQ0 >= 4 && ((MQ0 / (1.0 * DP)) > 0.1)' --filter-name 'HARD_TO_VALI' \
                    --filter-expression 'QUAL < 30.0'                            --filter-name 'VeryLowQual' \
                    --filter-expression 'QUAL >= 30.0 && QUAL < 50.0'            --filter-name 'LowQual' \
                    --genotype-filter-expression 'DP < 10'                        --genotype-filter-name 'LowDP' \
                    -O combined_genotyped_filtered_snps.vcf

date

./gatk-4.1.0.0/gatk --java-options "-Xmx4G" SelectVariants \
                    -R human_g1k_v37_decoy.fasta \
                    -V combined_genotyped.vcf \
                    --select-type-to-include INDEL \
                    -O combined_genotyped_raw_indels.vcf

# ここにMIXEDが必要
# INDEL, SNP, MIXED, MNP, SYMBOLIC, NO_VARIATION
# MIXED と MNP はテストしておくほうが良いな
date

./gatk-4.1.0.0/gatk --java-options "-Xmx4G" VariantFiltration \
                    -R human_g1k_v37_decoy.fasta \
                    -V combined_genotyped_raw_indels.vcf \
                    --filter-expression 'QD < 2.0'               --filter-name 'LowQD' \
                    --filter-expression 'FS > 200.0'             --filter-name 'HighFisherStrand' \
                    --filter-expression 'ReadPosRankSum < -20.0' --filter-name 'LowReadPosRankSum' \
                    --genotype-filter-expression 'DP < 10'        --genotype-filter-name 'LowDP' \
                    -O combined_genotyped_filtered_indels.vcf

date

picard MergeVcfs \
       I=combined_genotyped_filtered_snps.vcf \
       I=combined_genotyped_filtered_indels.vcf \
       O=combined_genotyped_filtered_snps_indels_mixed.vcf
# MergeVcfs https://software.broadinstitute.org/gatk/documentation/tooldocs/4.0.7.0/picard_vcf_MergeVcfs.php

date

./gatk-4.1.0.0/gatk --java-options "-Xmx4G" SelectVariants \
                    -R human_g1k_v37_decoy.fasta \
                    -V combined_genotyped_filtered_snps_indels_mixed.vcf \
                    --exclude-filtered --exclude-non-variants \
                    -O combined_genotyped_filtered_snps_indels_mixed.PASS.vcf

date
# 以上で germline variant call が完了です!

# combined_genotyped_filtered_snps_indels_mixed.PASS.vcf は多サンプルのデータが集合した .vcf なので、ここで個別サンプルごとの .avinput (Annovarの入力ファイル) に変換する
./annovar/convert2annovar.pl -format vcf4 \
                             --includeinfo \
                             --withzyg \
                             --allsample \
                             combined_genotyped_filtered_snps_indels_mixed.PASS.vcf \
                             --outfile combined_genotyped_filtered_snps_indels_mixed.PASS

date

for avinput in combined_genotyped_filtered_snps_indels_mixed.PASS.*.avinput; do
  id=`echo ${avinput} | cut -d. -f3`
  echo "Process: ${id} using Annovar"
  ./annovar/table_annovar.pl combined_genotyped_filtered_snps_indels_mixed.PASS.${id}.avinput annovar/humandb/ \
                            -buildver hg19 \
                            -protocol refGeneWithVer,genomicSuperDups,cosmic87_coding,cosmic87_noncoding,clinvar_20180603,generic,exac03,gnomad_genome,avsnp150,ljb26_all \
                            -genericdb tommo-3.5kjpnv2-20180625-af_snvall.MAF.genericdb \
                            -operation g,r,f,f,f,f,f,f,f,f \
                            -nastring NA \
                            --otherinfo \
                            --argument '--hgvs --exonicsplicing --splicing_threshold 2',,,,,,,,, \
                            --remove \
                            -out ${id}.avoutput2
done

date

grep -wF -e Func.refGeneWithVer -e exonic -e splicing ${id}.avoutput2.hg19_multianno.txt | grep -vwF -e "synonymous SNV" > ${id}.avoutput2.hg19_multianno.exonic.txt

#cat ${id}.avoutput2.hg19_multianno.exonic.txt | perl -F"\t" -lane 'print $_ if $F[11] <= 0.0001 || $. == 1' > ${id}.avoutput2.hg19_multianno.exonic.filtered_1.txt
#cat ${id}.avoutput2.hg19_multianno.exonic.txt | perl -F"\t" -lane 'print $_ if $. == 1 || ($F[11] <= 0.0001 && $F[14] <= 0.0001)' | grep -wF -e Chr -e hom | grep -vwF LowDP > ${id}.avoutput2.hg19_multianno.exonic.filtered_2.txt

echo "Finish!"