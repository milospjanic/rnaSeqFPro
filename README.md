# rnaSeqFPro

rnaSeqFPro is a script that will do full processing of paired RNA-Seq data starting from fastq.gz files placed in the same folder. Script will sort files and process paired .fastq.gz files. rnaSeqFPro will perform Fastqc quality control, it will map **paired fastq files** to the reference genome hg19, 

#Dependencies

**FastQC**

Instalation (Linux):

<pre>
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
unzip fastqc_v0.11.5.zip
chmod 755 ./FastQC/fastqc
cp ./FastQC/fastqc /usr/local/bin 
</pre>

**STAR**

Instalation (Linux):

<pre>
# Get the latest STAR source
wget https://github.com/alexdobin/STAR/archive/2.5.2b.tar.gz
tar -xzf 2.5.2b.tar.gz
cd STAR-2.5.2b

# Build STAR
make STAR

# If you have a TeX environment, build the documentation
make manual

chmod 755 STAR
cp STAR /usr/local/bin
</pre>

**Reference genome**

Download the reference genome, in this example it is human hg19:

<pre>
mkdir ~/reference_genomes
cd ~/reference_genomes
wget --timestamping 
        'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit ' 
        -O hg19.2bit 
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
chmod 755 twoBitToFa
./twoBitToFa hg19.2bit hg19.fa
</pre>

**Indexing the reference genome**

Use STAR to index the reference genome, use number of core on your machine, e.g. 64.
<pre>
cd ~/reference_genomes
STAR  --runMode genomeGenerate --runThreadN 64 --genomeDir ./ --genomeFastaFiles hg19.fa
</pre>

