# rnaSeqFPro

rnaSeqFPro is a script that will do full processing of paired RNA-Seq data starting from fastq.gz files placed in the same folder. Script will sort files and process paired .fastq.gz files. rnaSeqFPro will perform Fastqc quality control, it will map **paired fastq files** to the reference genome hg19 using STAR's second pass mapping.

#Dependencies

**Place fastqc.gz in a working folder**

<pre>
mkdir work.folder
cp path-to-files/*fastq.gz work.folder
</pre>

**FastQC**

Instalation (Linux):

<pre>
cd work.folder
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
unzip fastqc_v0.11.5.zip
chmod 755 ./FastQC/fastqc
cp ./FastQC/fastqc /usr/local/bin #may not work run it locally via link:
ln -s ./FastQC/fastqc .
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
mkdir hg19
cd hg19
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

**Download GENCODE transcript annotation**

For example for human hg19 genome:
<pre>
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gtf.gz
</pre>

**Install featureCounts**

Install featureCounts. Downoad Subread binary from Sourceforge.
<pre>
wget https://sourceforge.net/projects/subread/files/subread-1.5.1/subread-1.5.1-Linux-x86_64.tar.gz/download
</pre>
Download fileMulti2TableMod1.
	
wget https://raw.githubusercontent.com/milospjanic/fileMulti2TableMod1/master/fileMulti2TableMod1.awk
	
In R install RGSEPD from Bioconductor.
	
source("https://bioconductor.org/biocLite.R")
biocLite("rgsepd")

In R install DESeq2
		
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

In R install goseq
		
source("https://bioconductor.org/biocLite.R")
biocLite("goseq")

Download Kallisto binary.
	
wget https://github.com/pachterlab/kallisto/releases/download/v0.43.0/kallisto_linux-v0.43.0.tar.gz

Creating meta data table is necessary for RGSEPD to perform analysis using DESeq2/goseq. Table 1. shows an example of a meta data sheet. Note that sample names must be shown without .fastq.gz extension.

file.name.1 A CONDITION1
file.name.2 A CONDITION1
file.name.3 B CONDITION2
file.name.3
B
CONDITION2

rnaSeqFPro is composed of four pipelines that will run a RGSEPD version on either human genome hg19 or mouse genome mm10, using either paired-end (PE) or single-read (SR) sequences. Four additional pipelines are available to run a Kallisto version: PE hg19, SR hg19, PE mm10, and SR mm10. After placing files in the working folder run the script that is suitable for your experiment, e.g: 

chmod 755 rnaSeqFPro.PE.hg19.sh
./rnaSeqFPro.PE.hg19.sh

