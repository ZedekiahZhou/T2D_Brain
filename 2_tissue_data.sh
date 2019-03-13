#! /bin/bash
#Created in 2018/3/22
#Modified in 2018/9/13
# extract data of each tissue

cd /data/MyProgram/Final_diabrain/0.raw

#gzip -kfd *.gz

cd ../1.clean
#cat ../0.raw/GTEx_Data_20160115_v7_RNAseq_RNASeQCv1.1.8_gene_reads.gct |sed '1,2d'>GTEx_v7_RNAseq_gene_reads.gct

#extract seq_gene list and seq_sample list
#head -1 GTEx_v7_RNAseq_gene_reads.gct|tr '\t' '\n'>RNAseq_sample_12767.list
#awk '{print $1 "\t" $2}' GTEx_v7_RNAseq_gene_reads.gct >id2symbol.tab

#extract sample data for each tissue
tissuename=`ls ./samid/`
tissuename=($tissuename)
n=$[`ls ./samid/|wc -l`-1]
mkdir data temp
echo -e "1\n2" > ./temp/header

#you should mkdir data temp and make the header file with 1&2
for i in `seq 0 $n`
do
	echo "now comes to $i ${tissuename[i]}"
	sed -i 's/\r//g' ./samid/${tissuename[i]}
	cat "./RNAseq_sample_12767.list"|grep -nf ./samid/${tissuename[i]}|cut -d ":" -f 1 >./temp/tail
	cat ./temp/header ./temp/tail >./temp/sam_num_$i;
	tmp=`cat ./temp/sam_num_$i|tr "\n" ">"|sed 's/>/,$/g'|sed 's/^/$/'|sed 's/$$//'|sed 's/,$//'`;
	echo '{print '$tmp'}'>./temp/sam_progfile_$i;
	#gawk -f ./temp/sam_progfile_$i ./GTEx_coding_gene_rpkm.gct > ./data/${tissuename[i]};
	awk -f ./temp/sam_progfile_$i ./GTEx_v7_RNAseq_gene_reads.gct > ./data/${tissuename[i]};
	echo "now finished ${tissuename[i]}";
done
