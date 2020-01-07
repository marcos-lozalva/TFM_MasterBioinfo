!#/bin/bash

# Install HISAT2
cd ~/programas
# OS-X
wget -c ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/downloads/hisat2-2.1.0-OSX_x86_64.zip
unzip hisat2-2.1.0-OSX_x86_64.zip
cd /usr/local/bin/
ln -s ~/programas/hisat2-2.1.0/hisat2* .

# Linux
wget -c http://ccb.jhu.edu/software/hisat2/dl/hisat2-2.1.0-Linux_x86_64.zip
unzip hisat2-2.1.0-Linux_x86_64.zip

cd ~/bin/
ln -s ~/programas/hisat2-2.1.0/hisat2* .

# Install STRINGTIE
cd ~/programas
# OS-X
wget -c http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.3b.OSX_x86_64.tar.gz
tar xzf stringtie-1.3.3b.OSX_x86_64.tar.gz
cd /usr/local/bin
ln -s ~/programas/stringtie-1.3.3b.OSX_x86_64/stringtie stringtie
# Linux
wget -c http://ccb.jhu.edu/software/stringtie/dl/stringtie-1.3.6.Linux_x86_64.tar.gz
tar xzf stringtie-1.3.6.Linux_x86_64.tar.gz 
cd ~/bin/
ln -s ~/programas/stringtie-1.3.6.Linux_x86_64/stringtie stringtie


# Install SAMTOOLS
tar xjf samtools-1.6.tar.bz2 
cd samtools-1.6
./configure
make
cd /usr/local/bin/
ln -s ~/programas/samtools-1.9/samtools samtools

# Install GFFCOMPARE
cd ~/programas
git clone https://github.com/gpertea/gclib
git clone https://github.com/gpertea/gffcompare
cd gffcompare
make release
cd ~/bin/
ln -s ~/programas/gffcompare/gffcompare

# Download reference index (grch38_snp_tran) 
wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38_snp_tran.tar.gz
# Download annotation reference 
wget ftp://ftp.ensembl.org/pub/release-97/gff3/homo_sapiens/Homo_sapiens.GRCh38.97.chr.gff3.gz
#cd ~/TFM/rnaseq_marcos/raw 
cd /media/marcos/20FA7178FA714B54/raw #por falta de espacio

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



# Assembly with STRINGTIE
mkdir assembly
echo 'STRINGTIE. Assembling mapped reads into transcripts...'
sleep 2
for i in $(ls ./*.bam)
do
	stringtie ${i} -l ${i%\.bam} -p 8 -G ~/TFM/rnaseq_marcos/input/Homo_sapiens.GRCh38.97.chr.gff3 -o assembly/${i%\.bam}_assembled.gtf
done;
echo 'Assembly completed'

# Merge transcripts
echo 'STRINGTIE. Merging transcripts...'
stringtie --merge -p 8 -G ~/TFM/rnaseq_marcos/input/Homo_sapiens.GRCh38.97.chr.gff3 -o  stringtie_merged.gtf ./mergelist.txt
echo 'Merging completed'
# Check out the transcripts
cat stringtie_merged.gtf | head

# Estimation of transcript abundances
cd assembly/
for i in $(ls *.gtf)
do
	echo 'Calculating transcript abundance from' ${i%_assembled.gtf}'...'
	stringtie -e -B -p 8 -G ../stringtie_merged.gtf -o ../ballgown/${i%_assembled.gtf}/${i%_assembled.gtf}.gtf ../${i%_assembled.gtf}.bam
	echo 'Estimation of transcript abundances from' ${i%_assembled.gtf} 'completed'
done;
echo 'Estimation of transcript abundances completed'
cd ../

# GFF compare
gffcompare -r ~/TFM/rnaseq_marcos/input/Homo_sapiens.GRCh38.97.chr.gff3 -G -o merged stringtie_merged.gtf
cat merged.stats

# Generation of count matrices to use them for DE in Deseq2
# Create file sample_lst.txt with 2 columns: sample ID	./ballgown/sample_ID/file.gtf
# Two tables are created, gene_count_matrix.csv and transcript_count_matrix.csv
python2.7 prepDE.py -i sample_lst.txt








