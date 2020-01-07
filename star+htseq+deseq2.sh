#!/bin/bash

search_dir="/media/tarr/ALMACENAJE/Investigadores/May/RNASeq/"

SMR="/media/tarr/ALMACENAJE/Investigadores/May/RNASeq/ANALISIS/"

merge="/home/tarr/Documentos/Bioinfo/Scripts/merge-paired-reads.sh"

bac16s="/usr/local/bin/sortmerna-2.1/rRNA_databases/silva-bac-16s-id90.fasta"
ibac16s="/usr/local/bin/sortmerna-2.1/index/silva-bac-16s-db"
bac23s="/usr/local/bin/sortmerna-2.1/rRNA_databases/silva-bac-23s-id98.fasta"
ibac23s="/usr/local/bin/sortmerna-2.1/index/silva-bac-23s-db"
arc16s="/usr/local/bin/sortmerna-2.1/rRNA_databases/silva-arc-16s-id95.fasta"
iarc16s="/usr/local/bin/sortmerna-2.1/index/silva-arc-16s-db"
arc23s="/usr/local/bin/sortmerna-2.1/rRNA_databases/silva-arc-23s-id98.fasta"
iarc23s="/usr/local/bin/sortmerna-2.1/index/silva-arc-23s-db"
euk18s="/usr/local/bin/sortmerna-2.1/rRNA_databases/silva-euk-18s-id95.fasta"
ieuk18s="/usr/local/bin/sortmerna-2.1/index/silva-euk-18s-db"
euk28s="/usr/local/bin/sortmerna-2.1/rRNA_databases/silva-euk-28s-id98.fasta"
ieuk28s="/usr/local/bin/sortmerna-2.1/index/silva-euk-28s"
rfam5s="/usr/local/bin/sortmerna-2.1/rRNA_databases/rfam-5s-database-id98.fasta"
irfam5s="/usr/local/bin/sortmerna-2.1/index/rfam-5s-db"
rfam58s="/usr/local/bin/sortmerna-2.1/rRNA_databases/rfam-5.8s-database-id98.fasta"
irfam58s="/usr/local/bin/sortmerna-2.1/index/rfam-5.8s-db"

unmerge="/home/tarr/Documentos/Bioinfo/Scripts/unmerge-paired-reads.sh"

bbdukres="/usr/local/bin/bbmap/resources"

cd ${search_dir}

for file in $(find *R1_001.fastq.gz -type f);
do
	FICHERO="${file}"
	NOMBRE="${FICHERO%_R1*}"
	echo $NOMBRE

	## Analysis of the sequences
	
	fastqc -t 8 --noextract \
	${NOMBRE}_R1_001.fastq.gz ${NOMBRE}_R2_001.fastq.gz
	
	## decompress files
	
	gunzip -c ${NOMBRE}_R1_001.fastq.gz > ${SMR}${NOMBRE}_1.fastq
	gunzip -c ${NOMBRE}_R2_001.fastq.gz > ${SMR}${NOMBRE}_2.fastq
	
	## merge pair reads

	sh ${merge} ${SMR}${NOMBRE}_1.fastq \
	${SMR}${NOMBRE}_2.fastq \
	${SMR}${NOMBRE}.fastq
	
	rm ${SMR}${NOMBRE}_1.fastq
	rm ${SMR}${NOMBRE}_2.fastq
	
	## remove rRNAS in case that there are contaminations
	
	sortmerna --ref \
	${bac16s},${ibac16s}:${bac23s},${ibac23s}:${arc16s},${iarc16s}:${arc23s},${iarc23s}:${euk18s},${ieuk18s}:${euk28s},${ieuk28s}:${rfam5s},${irfam5s}:${rfam58s},${irfam58s} \
	--reads ${SMR}${NOMBRE}.fastq --fastx --num_alignments 1 -m 64000 \
	--aligned ${SMR}${NOMBRE}_rRNA --other ${SMR}${NOMBRE}_non_rRNA \
	-a 32 --paired_in -v --log	

	## unmerge

	sh ${unmerge} ${SMR}${NOMBRE}_non_rRNA.fastq \
	${SMR}${NOMBRE}_1_sort.fastq \
	${SMR}${NOMBRE}_2_sort.fastq 
	
	rm ${SMR}${NOMBRE}_non_rRNA.fastq
	
	## gzip
		
	pigz -f -p 32 ${SMR}${NOMBRE}_1_sort.fastq
	pigz -f -p 32 ${SMR}${NOMBRE}_2_sort.fastq
	
	##trimming adapters
		
	bbduk.sh in1=${SMR}${NOMBRE}_1_sort.fastq.gz \
	in2=${SMR}${NOMBRE}_2_sort.fastq.gz \
	ref=${bbdukres}/adapters.fa \
	out1=${SMR}clean1.SMR.${NOMBRE}.fastq.gz \
	out2=${SMR}clean2.SMR.${NOMBRE}.fastq.gz \
	tbo tpe mink=20 ktrim=r
	
	rm ${SMR}${NOMBRE}_1_sort.fastq.gz
	rm ${SMR}${NOMBRE}_2_sort.fastq.gz

	#Pre-processing files (Q=5, 10 y 20 y m=20)
	cutadapt -q 5,5 --pair-filter=any -m 20 \
	--output=${SMR}out.${NOMBRE}.1_Q5.SMR.fastq.gz \
	--paired-output=${SMR}out.${NOMBRE}.2_Q5.SMR.fastq.gz \
	${SMR}clean1.SMR.${NOMBRE}.fastq.gz \
	${SMR}clean2.SMR.${NOMBRE}.fastq.gz

	cutadapt -q 20,20 --pair-filter=any -m 20 \
	--output=${SMR}out.${NOMBRE}.1_Q20.SMR.fastq.gz \
	--paired-output=${SMR}out.${NOMBRE}.2_Q20.SMR.fastq.gz \
	${SMR}clean1.SMR.${NOMBRE}.fastq.gz \
	${SMR}clean2.SMR.${NOMBRE}.fastq.gz
	
	rm ${SMR}clean1.SMR.${NOMBRE}.fastq.gz
	rm ${SMR}clean2.SMR.${NOMBRE}.fastq.gz
	
done


G="/media/tarr/ALMACENAJE/Investigadores/TARIN/ANALISIS/SMR/RESULTS/GENOME/"

genomedir="/media/tarr/ALMACENAJE/Datos/RNAseq/GENOME/"
GRch38="/media/tarr/ALMACENAJE/Datos/EMBL/GRch38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
gtf="/media/tarr/ALMACENAJE/Datos/RNAseq/Homo_sapiens.GRCh38.96.gtf"

cd ${SMR}

for file in $(find *_Q5.SMR.fastq.gz -type f);
do
	FICHERO="${file}"
	NOMBRE="${FICHERO%_Q5*}"
	echo $NOMBRE
	
	STAR --genomeDir ${genomedir} --readFilesCommand zcat \
	--readFilesIn out.${NOMBRE}.1_Q5.SMR.fastq.gz \
	out.${NOMBRE}.1_Q5.SMR.fastq.gz \
	--readFilesCommand zcat \
	--runThreadN 40 \
	--outSAMstrandField intronMotif --sjdbGTFfile ${gtf} \
	--outFileNamePrefix RESULTS/${NOMBRE}.Q5.STAR \
	--outReadsUnmapped Fastx 
	
	STAR --genomeDir ${genomedir} --readFilesCommand zcat \
	--readFilesIn out.${NOMBRE}.1_Q20.SMR.fastq.gz \
	out.${NOMBRE}.1_Q20.SMR.fastq.gz \
	--readFilesCommand zcat \
	--runThreadN 40 \
	--outSAMstrandField intronMotif --sjdbGTFfile ${gtf} \
	--outFileNamePrefix RESULTS/${NOMBRE}.Q5.STAR \
	--outReadsUnmapped Fastx 

done

cd RESULTS
for file in $(find *.Q5.STARAligned.out.sam -type f);
do
	FICHERO="${file}"
	NOMBRE="${FICHERO%.Q5.STARAligned.out.sam}"
	echo $NOMBRE
	
	samtools view -@ 40 -Sb ${file} > ${NOMBRE}.Q5.bam
	
	samtools sort -@ 40 ${NOMBRE}.Q5.bam \
	-o ${NOMBRE}.Q5.sorted.bam
	
	rm ${G}${NOMBRE}.Q5.bam
	
	samtools index -@ 40 ${NOMBRE}.Q5.sorted.bam
	
	samtools view -@ 40 -Sb ${file} > ${NOMBRE}.Q20.bam
	
	samtools sort -@ 40 ${NOMBRE}.Q20.bam \
	-o ${NOMBRE}.Q20.sorted.bam
	
	rm ${NOMBRE}.Q20.bam
	
	samtools index -@ 40 ${NOMBRE}.Q20.sorted.bam
	
done


cd ${G}
	
done
 
for file in $(find *Q20.sorted.bam -type f);
do
	FICHERO="${file}"
	NOMBRE="${FICHERO%.Q20.sorted.bam}"
	echo $NOMBRE
	
	python -m HTSeq.scripts.count -m intersection-nonempty -f bam -r pos -s no -t exon -i gene_id ${file} ${gtf} > ${G}${NOMBRE}.Q20.txt 

done

for file in $(find *Q5.sorted.bam -type f);
do
	FICHERO="${file}"
	NOMBRE="${FICHERO%.Q5.sorted.bam}"
	echo $NOMBRE
	
	python -m HTSeq.scripts.count -m intersection-nonempty -f bam -r pos -s no -t exon -i gene_id ${file} ${gtf} > ${G}${NOMBRE}.Q5.txt 

done
