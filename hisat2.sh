!#/bin/bash

# Install HISAT2
cd ~/programas
# OS-X
wget -c ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-OSX_x86_64.zip
unzip hisat2-2.1.0-OSX_x86_64.zip
cd /usr/local/bin/
ln -s ~/programas/hisat2-2.1.0/hisat2* .


# Install SAMTOOLS
tar xjf samtools-1.6.tar.bz2 
cd samtools-1.6
./configure
make
cd /usr/local/bin/
ln -s ~/programas/samtools-1.9/samtools samtools

# Download reference index (grch38_snp_tran) 
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38_snp_tran.tar.gz

# Mapping with HISAT2
echo 'Mapping with HISAT2...'
sleep 2
for i in $(ls | grep 'R1')
do
	echo 'Aligning' ${i%_S*}'...'
	R2=$(echo ${i} | sed 's/R1/R2/g')
 	hisat2 -p 8 -x ~/TFM/rnaseq_marcos/input/grch38_snp_tran/genome_snp_tran -1 ./${i} -2 ./$R2 -S ./${i%_S*}.sam 
	echo  ${i%_S*}.sam 'to' ${i%_S*}.bam 'conversion...'
	samtools sort -o ${i%_S*}.bam ./${i%_S*}.sam
	rm ${i%_S*}.sam #ahorro de espacio
	echo 'File' ${i%_S*}.sam 'eliminated'
	echo ${i%_S*} 'alignment completed'
	echo ''
	echo '......................................'
	echo '' 
done;
