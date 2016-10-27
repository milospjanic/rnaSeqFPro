#!/bin/bash

#fastqc quality control - requires fastqc installed and placed in PATH

ls -1 *fastq.gz > commands.1
sed -i 's/^/.\/FastQC\/fastqc /g' commands.1

source commands.1

#mapping with STAR - requires STAR installed and copied to PATH


files=(*fastq.gz)
for (( i=0; i<${#files[@]} ; i+=2 )) ; do
    mkdir "${files[i]}.${files[i+1]}.STAR"    
done 

files=(*fastq.gz)
for (( i=0; i<${#files[@]} ; i+=2 )) ; do
  cat >> commands.2.${files[i]}.${files[i+1]}.tmp <<EOL
#!/bin/bash
 GenomeDir="~/reference_genomes/"
 GenomeFasta="~/reference_genomes/hg19.fa"
 CommonPars="--runThreadN 64 --outSAMattributes All --genomeLoad NoSharedMemory"
    echo Proccessing `pwd`: ${files[i]} ${files[i+1]}
     
    # run 1st pass
        mkdir Pass1
        cd Pass1
        STAR ${CommonPars} --genomeDir ${GenomeDir} --readFilesIn ${Reads}
        cd ..
    # make splice junctions database file out of SJ.out.tab, filter out non-canonical junctions
        mkdir GenomeForPass2
        cd GenomeForPass2
        awk 'BEGIN {OFS="\t"; strChar[0]="."; strChar[1]="+"; strChar[2]="-";} {if(\$5>0){print \$1,\$2,\$3,strChar[\$4]}}' ../Pass1/SJ.out.tab > SJ.out.tab.Pass1.sjdb
    # generate genome with junctions from the 1st pass
        STAR --genomeDir ./ --runMode genomeGenerate --genomeFastaFiles ${GenomeFasta} --sjdbFileChrStartEnd SJ.out.tab.Pass1.sjdb --sjdbOverhang 100 --runThreadN 64
        cd ..
    # run 2nd pass with the new genome
        mkdir Pass2
        cd Pass2
        STAR ${CommonPars} --genomeDir ../GenomeForPass2 --readFilesIn ${Reads}
        echo FINISHED ${Reads}
        cd ..
        cd ..
        done
EOL
  done

files=(*fastq.gz)
for (( i=0; i<${#files[@]} ; i+=2 )) ; do
    sed -i "8i\\\tcd ${files[i]}.${files[i+1]}.STAR" commands.2.${files[i]}.${files[i+1]}.tmp
done

for (( i=0; i<${#files[@]} ; i+=2 )) ; do
    sed -i "2i\ Reads=\"`pwd`/${files[i]} `pwd`/${files[i+1]} --readFilesCommand zcat\"" commands.2.${files[i]}.${files[i+1]}.tmp
done

for (( i=0; i<${#files[@]} ; i+=2 )) ; do
    source commands.2.${files[i]}.${files[i+1]}.tmp
done

#subscrips

files=(*fastq.gz)
for (( i=0; i<${#files[@]} ; i+=2 )) ; do
    echo "${files[i]}" "${files[i+1]}" >> commands.2
done

touch sam.tmp

files=(*fastq.gz)
for (( i=0; i<${#files[@]} ; i+=2 )) ; do
    echo "${files[i]}.${files[i+1]}.sam" >> sam.tmp     
done

awk 'FNR==NR{a[FNR]=$0;next}{ print $0,">",a[FNR]}' sam.tmp commands.2 > merge.tmp

rm sam.tmp
rm commands.2
rm merge.tmp

#counting with featureCounts

find . -type f -wholename "*Pass2*sam" -exec sh -c '
    for f
        do echo $f
        fileName=$(basename $f);
        filePath=$(dirname $f);
        lastDir=$(basename $filePath);
        prevDir=$(basename $(dirname $filePath));
        echo $prevDir
        echo $lastDir
        echo proccessing  $fileName from $(pwd)/$lastDir into $prevDir.counts.txt;
        featureCounts -a gencode.v25lift37.annotation.gtf.gz -o $prevDir.counts.txt -T 12 -t exon -g gene_id $f
    done' sh {} +
