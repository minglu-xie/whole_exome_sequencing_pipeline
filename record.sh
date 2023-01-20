#################################################
#  File Name:record.sh
#  Author: Minglu.Xie
#  Mail: minglu.xie.1023@student.uu.se,mingluxie@gmail.com
#  Created Time: Fri 06 Jan 2023 01:07:36 PM CET
#################################################

# create the environment by conda/mamba.
mamba create -n wes
mamba activate wes
mamba install bcftools snpEFF multiqc qualimap fastqc sra-tools

#Download the 'Accession List' from the https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA643767. Here we only want to download the whole exome-sequencing data.
# $cat SRR_Acc_List.txt
# SRR12135911
# SRR12135912

#before use the prefetch command, please install sratools from NCBI.
cat SRR_Acc_List.txt|while read id;
do nohup prefetch ${id} &
done

#the downloaded data will be save in the location where you installed your sratools. Move data to the folder that you created for your project.
mkdir /disk1/minglu/wges/project/GSE153707
cd /disk1/minglu/wges/project/GSE153707
path1=`which prefetch`
path2=`echo ${path1%bin*}`
ls $path2/down/sra
cat SRR_Acc_List.txt|while read id;
do
mv $path2/down/sra/${id}* ./
done

mkdir 1.fastq
#transfer from sra to fastq
for id in `ls *sra*`;
do
sample=${id%%.*}
nohup fastq-dump --gzip --split-3 -O ./1.fastq ${id} 1>${sample}.log 2>&1 &
done

#quality control, use xargs command
mkdir 3.mapping
ls 1.fastq/*|xargs fastqc -t 10 -o 3.mapping
multiqc ./3.mapping -o ./3.mapping


# ##sometimes, the data quality is not good. fg:adapter content. We may need to trim the adapter
# dir=/disk1/minglu/wges/project/sra/00.trim_data/
# cat config |while read id
# do
# arr=($id)
# fq1=${arr[1]}
# fq2=${arr[2]}
# sample=${arr[0]}
# nohup trim_galore -q 25 --phred33 --gzip --length 36 -e 0.1 --stringency 3 --paired -o $dir $fq1 $fq2 2>${sample}_trim_err 1>&2 &
# done
# #only can use 1 core, otherwise will cause error. The reason may be due to the python version. Python -2.7

#configure sample information, generate the config file which include the information for following analysis
ls 1.fastq/*1.fastq.gz >fq1
ls 1.fastq/*2.fastq.gz >fq2
paste list.txt fq1 fq2 >sample.txt
rm fq1 fq2

# $cat sample.txt 
# SRR12135911	1.fastq/SRR12135911_1.fastq.gz	1.fastq/SRR12135911_2.fastq.gz
# SRR12135912	1.fastq/SRR12135912_1.fastq.gz	1.fastq/SRR12135912_2.fastq.gz

#bwa mapping
#index must be build up by bwa index
INDEX=/disk1/pengweixing/database/hg38/index/hg38.fa
cat sample.txt | while read id
do
arr=(${id})
fq1=${arr[1]}
fq2=${arr[2]}
sample=${arr[0]}
bwa mem -t 20 -R '@RG\tID:$sample\tSM:$sample\tLB:WGS\tPL:Illumina' $INDEX $fq1 $fq2 | samtools sort -@ 10 -O BAM -o ./3.mapping/$sample.sort.bam &
done

# #for plot-bamstats, you may need to install gnuplot
# mkdir 4.align
# INDEX=/disk1/pengweixing/database/hg38/index/hg38.fa
# cat sample.txt|while read id;
# do
# arr=($id)
# sample=${arr[0]}
# bam=./3.mapping/${sample}_sort.bam
# samtools stats -@ 15 --reference $INDEX ${bam} > ./4.align/stats/${sample}.stat
# plot-bamstats -p ./4.align/stats/${sample} ./4.align/stats/${sample}.stat
# done


# #qualimap to check coverage & depth
# cat sample.txt | while read id
# do
# arr=($id)
# sample=${arr[0]}
# bam=./3.mapping/${sample}_sort.bam
# qualimap bamqc --java-mem-size=10G -gff /disk1/minglu/wges/biosoft/GATK_resources_hg38/hg38.exon.bed -nr 100000 -nw 500 -nt 16 -bam ${bam} -outdir ./4.align/qualimap/${sample}
# done
# -nr <arg>     Number of reads analyzed in a chunk
#               (default is 1000)
# -nt <arg>     Number of threads (default is 8)
# -nw <arg>     Number of windows (default is 400)


##Before use GATK, we need to download lots of data from their website.
##https://gatk.broadinstitute.org/hc/en-us
##https://console.cloud.google.com/storage/browser/genomics-public-data/resources/broad/hg38/v0;tab=objects?prefix=&forceOnObjectsSortingFiltering=false
mkdir 5.gatk

# Markduplicates
cat SRR_Acc_List.txt|while read sample
do
bam=./3.mapping/${sample}.sort.bam
echo "start MarkDuplicates for ${sample}" `date`
/disk1/minglu/wges/biosoft/gatk4/gatk-4.3.0.0/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" MarkDuplicates -I ${bam} -O ./5.gatk/${sample}_marked.bam 	-M ./5.gatk/${sample}.metrics 	1>./5.gatk/${sample}_log.mark 2>&1 &
done


# FixMateInformation
# 37.18 minutes for 7.7G BAM
cd ./5.gatk
ls *_marked.bam|while read sample;
do
/disk1/minglu/wges/biosoft/gatk4/gatk-4.3.0.0/gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" FixMateInformation -I ${sample} -O ${sample%%_*}_marked_fixed.bam -SO coordinate 1>${sample%%_*}_fixed.log 2>&1 &
done

#build index for *_fixed.bam
ls *fixed.bam|xargs -i samtools index -@ 10 {} &

#BaseRecalibrator:Generates recalibration table for Base Quality Score Recalibration (BQSR)
ref=/disk1/minglu/wges/biosoft/GATK_resources_hg38/Homo_sapiens_assembly38.fasta
snp=/disk1/minglu/wges/biosoft/GATK_resources_hg38/Homo_sapiens_assembly38.dbsnp138.vcf
indel=/disk1/minglu/wges/biosoft/GATK_resources_hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
bed=/disk1/minglu/wges/biosoft/GATK_resources_hg38/hg38.exon.bed
germine_vcf=/disk1/minglu/wges/biosoft/GATK_resources_hg38/af-only-gnomad.hg38.vcf.gz
GATK=/disk1/minglu/wges/biosoft/gatk4/gatk-4.3.0.0/gatk

ls *_marked_fixed.bam|while read sample
do 
$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" \
         BaseRecalibrator \
         -R $ref \
         -I $sample \
         --known-sites $snp \
         --known-sites $indel \
         -O ${sample%%_*}_recal.table \
         1>${sample%%_*}_log.recal 2>&1 &
done

##ApplyBQSR
## Recalibrate the base qualities of the input reads based on the recalibration table produced by the BaseRecalibrator tool, and outputs a recalibrated BAM file.
ls *_marked_fixed.bam|while read sample
do 
$GATK --java-options "-Xmx20G -Djava.io.tmpdir=./"   ApplyBQSR \
         -R $ref  \
         -I $sample  \
         -bqsr ${sample%%_*}_recal.table \
         -O ${sample%%_*}_bqsr.bam \
         1>${sample%%_*}_log.ApplyBQSR  2>&1 &
done

#check the results
for sample in  `ls *bqsr.bam`
do
id=${sample%%_*}
samtools flagstat -@ 10 $sample > ${id}_bqsr.stat &
done


# ##HaplotypeCaller: Call germline SNPs and indels via local re-assembly of haplotypes
# ls *bqsr.bam|while read sample
# do
# $GATK --java-options "-Xmx20G -Djava.io.tmpdir=./" HaplotypeCaller \
#          -R $ref  \
#          -I $sample \
#          --dbsnp $snp \
#          -O ./vcf/${sample%%_*}_raw.vcf \
#          1>./vcf/${sample%%_*}_log.HC 2>&1
# done


## Mutect2
### Call somatic SNVs and indels via local assembly of haplotypes
mkdir Mutect2
$GATK --java-options "-Xmx20G" Mutect2 \
-R $ref -I SRR12135912_bqsr.bam \
-tumor SRR12135912 \
-I SRR12135911_bqsr.bam \
-normal SRR12135911 \
--germline-resource $germine_vcf \
--af-of-alleles-not-in-resource 0.0000025 \
--bam-output HQ461-untreated-Mutect2.bam \
-O ./Mutect2/HQ461-untreated.Mutect2.vcf

$GATK SelectVariants \
-R $ref -V $germine_vcf \
--select-type-to-include SNP \
--restrict-alleles-to BIALLELIC \
-O /disk1/minglu/wges/biosoft/GATK_resources_hg38/af-only-gnomad.hg38.SNP_biallelic.vcf.gz

germine_biallelic_vcf=/disk1/minglu/wges/biosoft/GATK_resources_hg38/af-only-gnomad.hg38.SNP_biallelic.vcf.gz
genomic_intervals=/disk1/minglu/wges/biosoft/GATK_resources_hg38/wgs_calling_regions.hg38.interval_list

$GATK GetPileupSummaries \
-I SRR12135911_bqsr.bam \
-L $genomic_intervals \
-V $germine_biallelic_vcf \
-O ./Mutect2/SRR12135911.pileups.table

$GATK GetPileupSummaries -I SRR12135912_bqsr.bam -L $genomic_intervals -V $germine_biallelic_vcf -O ./Mutect2/SRR12135912.pileups.table

$GATK GetPileupSummaries -I SRR12135912_bqsr.bam -L $genomic_intervals -V $germine_biallelic_vcf -O ./Mutect2/SRR12135912.pileups.table
wait

$GATK CalculateContamination \
-I ./Mutect2/SRR12135912.pileups.table \
-matched ./Mutect2/SRR12135911.pileups.table \
-O ./Mutect2/HQ461-untreated.calculatecontamination.table \

$GATK FilterMutectCalls \
-R $ref \
-V ./Mutect2/HQ461-untreated.Mutect2.vcf \
--contamination-table ./Mutect2/HQ461-untreated.calculatecontamination.table \
-O ./Mutect2/HQ461-untreated.filtered.Mutect2.vcf

##install snpEFF
#download the database
nohup java -jar /disk1/minglu/software/miniconda3/envs/wes/share/snpeff-5.1-2/snpEff.jar download GRCh38.99 &
nohup java -jar /disk1/minglu/software/miniconda3/envs/wes/share/snpeff-5.1-2/snpEff.jar download hg38 &

####snpEFF annotation
java -Xmx20g -jar /disk1/minglu/software/miniconda3/envs/wes/share/snpeff-5.1-2/snpEff.jar -v GRCh38.99 ./HQ461-untreated.filtered.Mutect2.vcf > HQ461-untreated.filtered.ann.vcf

#### Keep those mutations that have a "pass" label
## Keep "HIGH" and "Moderate" record under IMPACT of the ANN field.
cat HQ461-untreated.filtered.ann.vcf | java -jar /disk1/minglu/software/miniconda3/envs/wes/share/snpsift-5.1-0/SnpSift.jar filter "( na FILTER ) | (FILTER = 'PASS')" > HQ461-untreated.filtered.ann.snpsift.vcf

cat HQ461-untreated.filtered.ann.snpsift.vcf | java -jar /disk1/minglu/software/miniconda3/envs/wes/share/snpsift-5.1-0/SnpSift.jar filter "(ANN[0].IMPACT has 'HIGH') | (ANN[0].IMPACT has 'MODERATE')" > HQ461-untreated.filtered.ann.snpsift.second.vcf

#check the results with the original publications
grep "CDK12" HQ461-untreated.filtered.ann.snpsift.second.vcf

# $grep "CDK12" HQ461-untreated.filtered.ann.snpsift.second.vcf
# chr17	39492834	.	G	A	.	PASS	AS_FilterStatus=SITE;AS_SB_TABLE=82,51|11,10;DP=157;ECNT=1;GERMQ=93;MBQ=20,20;MFRL=231,224;MMQ=60,60;MPOS=40;NALOD=1.81;NLOD=18.63;POPAF=5.60;TLOD=40.07;ANN=A|missense_variant|MODERATE|CDK12|ENSG00000167258|transcript|ENST00000447079.5|protein_coding|4/14|c.2192G>A|p.Gly731Glu|2778/8287|2192/4473|731/1490||,A|missense_variant|MODERATE|CDK12|ENSG00000167258|transcript|ENST00000430627.6|protein_coding|4/14|c.2192G>A|p.Gly731Glu|2503/5918|2192/4446|731/1481||,A|missense_variant|MODERATE|CDK12|ENSG00000167258|transcript|ENST00000581593.1|protein_coding|2/3|c.332G>A|p.Gly111Glu|334/559|332/557|111/184||WARNING_TRANSCRIPT_INCOMPLETE,A|missense_variant|MODERATE|CDK12|ENSG00000167258|transcript|ENST00000584632.5|protein_coding|4/13|c.2189G>A|p.Gly730Glu|2750/4165|2189/3604|730/1200||WARNING_TRANSCRIPT_INCOMPLETE  GT:AD:AF:DP:F1R2:F2R1:FAD:SB	0/0:85,0:0.015:85:24,0:28,0:62,0:52,33,0,0	0/1:48,21:0.289:69:15,7:15,5:39,15:30,18,11,10
# chr17	39494650	.	T	C	.	PASS	AS_FilterStatus=SITE;AS_SB_TABLE=54,150|2,5;DP=217;ECNT=1;GERMQ=93;MBQ=21,28;MFRL=239,301;MMQ=60,60;MPOS=36;NALOD=1.89;NLOD=22.53;POPAF=5.60;TLOD=10.22;ANN=C|missense_variant|MODERATE|CDK12|ENSG00000167258|transcript|ENST00000447079.5|protein_coding|5/14|c.2375T>C|p.Ile792Thr|2961/8287|2375/4473|792/1490||,C|missense_variant|MODERATE|CDK12|ENSG00000167258|transcript|ENST00000430627.6|protein_coding|5/14|c.2375T>C|p.Ile792Thr|2686/5918|2375/4446|792/1481||,C|missense_variant|MODERATE|CDK12|ENSG00000167258|transcript|ENST00000581593.1|protein_coding|3/3|c.515T>C|p.Ile172Thr|517/559|515/557|172/184||WARNING_TRANSCRIPT_INCOMPLETE,C|missense_variant|MODERATE|CDK12|ENSG00000167258|transcript|ENST00000584632.5|protein_coding|5/13|c.2372T>C|p.Ile791Thr|2933/4165|2372/3604|791/1200||WARNING_TRANSCRIPT_INCOMPLETE   GT:AD:AF:DP:F1R2:F2R1:FAD:SB	0/0:97,0:0.013:97:21,0:24,0:75,0:27,70,0,0	0/1:107,7:0.083:114:45,1:24,6:86,7:27,80,2,5