# rnaSeqFPro

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
# Get latest STAR source from releases
wget https://github.com/alexdobin/STAR/archive/2.5.2b.tar.gz
tar -xzf 2.5.2b.tar.gz
cd STAR-2.5.2b

# Build STAR
make STAR

# To include STAR-Fusion
git submodule update --init --recursive

# If you have a TeX environment, you may like to build the documentation
make manual
</pre>
