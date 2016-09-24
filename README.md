# rnaSeqFPro

rnaSeqFPro is a script that will do full proccessing of RNA-Seq data starting from
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
wget --timestamping 
        'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit ' 
        -O hg19.2bit 
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa
chmod 755 twoBitToFa
./twoBitToFa hg19.2bit hg19.fa
</pre>

        
