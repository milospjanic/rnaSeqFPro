#!/bin/bash

#fastqc quality control

ls -1 *fastqc.gz > commands
sed -i 's/^/fastqc /g' commands
source commands

#mapping
files=(*fastq.gz)
for (( i=0; i<${#files[@]} ; i+=2 )) ; do
    echo "${files[i]}" "${files[i+1]} \n" >> commands.2
done

sed -i 's/^/bwa\ mem\ -t\ 64\ /g' commands.2
