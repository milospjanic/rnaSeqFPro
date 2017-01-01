# rnaSeqFPro

rnaSeqFPro is a script that will do full processing of paired RNA-Seq data starting from fastq.gz files placed in the same folder. Script will sort files and process paired .fastq.gz files. rnaSeqFPro will perform Fastqc quality control, it will map **paired fastq files** to the reference genome hg19 using STAR's second pass mapping.

#Dependencies

**Place fastqc.gz in a working folder**

<pre>
mkdir work.folder
cp path-to-files/*fastq.gz work.folder
</pre>

**FastQC**

Instalation (Linux), place FastQC folder in working directory:

<pre>
cd work.folder
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip
unzip fastqc_v0.11.5.zip
chmod 755 ./FastQC/fastqc
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

For the mouse genome:

<pre>
mkdir mm10
cd mm10
wget --timestamping 
        ' http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.2bit' 
        -O mm10.2bit 	
./twoBitToFa mm10.2bit mm10.fa
</pre>

**Indexing the reference genome**

Use STAR to index the reference genome, use number of core on your machine, e.g. 64.
<pre>
cd ~/reference_genomes
STAR  --runMode genomeGenerate --runThreadN 64 --genomeDir ./ --genomeFastaFiles hg19.fa
</pre>

**Download GENCODE transcript annotation to the working folder**

For example for human hg19 genome:
<pre>
cd work.folder
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_25/GRCh37_mapping/gencode.v25lift37.annotation.gtf.gz
</pre>

For mouse mm10 genome:
<pre>
cd work.folder
wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M11/gencode.vM11.annotation.gtf.gz
</pre>

**Install featureCounts**

Install featureCounts. Downoad Subread binary from Sourceforge.
<pre>
wget https://sourceforge.net/projects/subread/files/subread-1.5.1/subread-1.5.1-Linux-x86_64.tar.gz/download
</pre>

**Download fileMulti2TableMod1**

<pre>
wget https://raw.githubusercontent.com/milospjanic/fileMulti2TableMod1/master/fileMulti2TableMod1.awk
</pre>	

**In R install RGSEPD from Bioconductor**

<pre>
source("https://bioconductor.org/biocLite.R")
biocLite("rgsepd")
</pre>

**In R install DESeq2**

<pre>
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
</pre>

**In R install goseq**

<pre>
source("https://bioconductor.org/biocLite.R")
biocLite("goseq")
</pre>

**Download Kallisto binary**

<pre>
wget https://github.com/pachterlab/kallisto/releases/download/v0.43.0/kallisto_linux-v0.43.0.tar.gz
</pre>

**Meta data**

Providing a meta information placed in a file named meta.data is necessary for RGSEPD to perform analysis using DESeq2/goseq. Table 1 shows an example of a meta data sheet for SE reads. **Note that sample names must be shown without .fastq.gz extension.**

<pre>
Sample  Condition       SHORTNAME
file.name.1 A CONDITION1
file.name.2 A CONDITION1
file.name.3 B CONDITION2
file.name.4 B CONDITION2
</pre>

For paired end reads meta data table should contain the names of paired files separated with period (.)

<pre>
Sample  Condition       SHORTNAME
file.name.1.R1.file.name.1.R2 A CONDITION1
file.name.2.R1.file.name.2.R2 A CONDITION1
file.name.3.R1.file.name.3.R2 B CONDITION2
file.name.4.R1.file.name.4.R2 B CONDITION2
</pre>

**Note: the order of the samples in meta.data has to the same as in command: ls -1 for the script to work. RGSEPD will stop if the orders do not match between meta.data and mastertable, therefore create meta.data in the same order as the mastertable is created, using ls -1 hierarchy**

**Meta data needs to be tab separated to avoid errors. Run tabsep on your file if necessary.
https://github.com/milospjanic/tabsep **

**Running**

rnaSeqFPro is composed of four pipelines that will run a RGSEPD version on either human genome hg19 or mouse genome mm10, using either paired-end (PE) or single-read (SR) sequences. 

<pre>
rnaSeqFPro.PE.hg19.sh
rnaSeqFPro.PE.mm10.sh
rnaSeqFPro.SR.hg19.sh
rnaSeqFPro.SR.mm10.sh
</pre>

Four additional pipelines are available to run a Kallisto version: PE hg19, SR hg19, PE mm10, and SR mm10. 

<pre>
rnaSeqFPro.PE.hg19.Kallisto.sh
rnaSeqFPro.PE.mm10.Kallisto.sh
rnaSeqFPro.SR.hg19.Kallisto.sh
rnaSeqFPro.SR.mm10.Kallisto.sh
</pre>

After placing files in the working folder run the script that is suitable for your experiment, e.g: 

<pre>
chmod 755 rnaSeqFPro.PE.hg19.sh
./rnaSeqFPro.PE.hg19.sh
</pre>

**If you are trying to re-run the pipe, e.g. you have already previously ran the pipe but for some reason it stopped, delete the GSEPD_OUTPUT folder as the RGSEPD needs to create this folder itself, it will stop if it encounters this folder already created** 

**Don't forget to place the GENCODE gtf file, FastQC folder and meta.data into the working folder! These are the 3 only requirments neccesary to be in the working directory, in addition to the fastqc.gz files**

**Dont forget that reference genome needs to be in your ~/reference_genomes folder, in case you switch to another user account script may not work because it searches for ~/reference_genomes folder.**
