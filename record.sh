#################################################
#  File Name:record.sh
#  Author: Minglu.Xie
#  Mail: minglu.xie.1023@student.uu.se,mingluxie@gmail.com
#  Created Time: Fri 06 Jan 2023 01:07:36 PM CET
#################################################

# create the environment by conda/mamba.
mamba create -n wes
mamba activate wes
mamba install bcftools snpEFF multiqc qualimap

#Download the 'Accession List' from the https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA643767. Here we only want to download the whole exome-sequencing data.
#$cat SRR_Acc_List.txt
#SRR12135911
#SRR12135912

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
