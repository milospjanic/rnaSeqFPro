#!/bin/bash

#fastqc quality control - requires fastqc installed and placed in PATH

ls -1 *fastqc.gz > commands
sed -i 's/^/fastqc /g' commands
source commands

#mapping with bwa mem - requires bwa mem installed and copied to PATH

files=(*fastq.gz)
for (( i=0; i<${#files[@]} ; i+=2 )) ; do
    echo "${files[i]}" "${files[i+1]} \n" >> commands.2
done

sed  -i 's/^/bwa\ mem\ -t\ 64\ /g' commands.2


