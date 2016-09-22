#!/bin/bash

#fastqc quality control - requires fastqc installed and placed in PATH

ls -1 *fastq.gz > commands.1
sed -i 's/^/fastqc /g' commands.1

#source commands.1

#mapping with bwa mem - requires bwa mem installed and copied to PATH

files=(*fastq.gz)
for (( i=0; i<${#files[@]} ; i+=2 )) ; do
    echo "${files[i]}" "${files[i+1]}" >> commands.2
done

sed -i 's/^/Reads="/g' commands.2
sed -i 's/$/--readFilesCommand zcat"' commands.2

files=(*fastq.gz)
for (( i=0; i<${#files[@]} ; i+=2 )) ; do
    echo "${files[i]}.${files[i+1]}.sam" >> sam.tmp     
done

files=(*fastq.gz)
for (( i=0; i<${#files[@]} ; i+=2 )) ; do
    mkdir "${files[i]}.${files[i+1]}.STAR"    
done 

files=(*fastq.gz)
for (( i=0; i<${#files[@]} ; i+=2 )) ; do
      cat >> commands.4.${files[i]}.${files[i+1]}.tmp <<EOL
        line 1, GenomeDir=/home/reference_genomes/
        line 2, GenomeFasta=/home/reference_genomes/hg19.fa
        line 3, CommonPars="--runThreadN 64 --outSAMattributes All --genomeLoad NoSharedMemory"
        line 4, echo Proccessing $Reads
        line 5, 
        line 6, # run 1st pass
        line 7, mkdir Pass1
        line 8, cd ${files[i]}.${files[i+1]}.STAR
        line 9, cd Pass1
        line 10, STAR $CommonPars --genomeDir $GenomeDir --readFilesIn $Reads
        line 11, cd ..

        line 12, # make splice junctions database file out of SJ.out.tab, filter out non-canonical junctions
        line 13, mkdir GenomeForPass2
        line 14, cd GenomeForPass2
        line 15, awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} {if($5>0){print $1,$2,$3,strChar[$4]}}' ../Pass1/SJ.out.tab > SJ.out.tab.Pass1.sjdb

        line 16, # generate genome with junctions from the 1st pass
        line 17, STAR --genomeDir ./ --runMode genomeGenerate --genomeFastaFiles $GenomeFasta --sjdbFileChrStartEnd SJ.out.tab.Pass1.sjdb --sjdbOverhang 100 --runThreadN 64
        line 18, cd ..

        line 19 # run 2nd pass with the new genome
        line 20, mkdir Pass2
        line 21, cd Pass2
        line 22, STAR $CommonPars --genomeDir ../GenomeForPass2 --readFilesIn $Reads

        line 23, echo FINISHED $Reads
        line 24, cd ..
        line 25, cd ..
        line 26, done

EOL

awk 'FNR==NR{a[FNR]=$0;next}{ print $0,">",a[FNR]}' sam.tmp commands.2
rm sam.tmp

#source commands.2

#counting with featurecount
