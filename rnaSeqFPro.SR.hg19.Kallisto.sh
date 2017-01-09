#!/bin/bash

SECONDS=0

###fastqc quality control - requires fastqc installed and placed in PATH

ls -1 *fastq.gz > commands.1
sed -i 's/^/.\/FastQC\/fastqc /g' commands.1

source commands.1

mkdir FastQC_OUTPUT
mv *zip FastQC_OUTPUT
mv *html FastQC_OUTPUT

###pseudomapping with Kallisto - requires Kallisto installed and copied to PATH

files=(*fastq.gz)
for (( i=0; i<${#files[@]} ; i+=1 )) ; do
    mkdir "${files[i]}.kallisto"    
done 

#check if Kallisto index is present, if not create one
FILE=GENCODE_transcripts_human 

if [ ! -f $FILE ]
then
kallisto index -i GENCODE_transcripts_human gencode.v25lift37.transcripts.fa.gz
fi

#pseudo-mapping with kallisto

GenomeDir='~/reference_genomes/hg19/'
GenomeFasta='~/reference_genomes/hg19/hg19.fa'

files=(*fastq.gz)
for (( i=0; i<${#files[@]} ; i+=1 )) ; do

echo $(pwd)/${files[i]}

Reads="$(pwd)/"${files[i]}" "
Index="-i $(pwd)/GENCODE_transcripts_human"
Parameters='--single -l 200 -s 20'
WorkDir="$(pwd)"

echo $Reads

  cat >> commands.2.${files[i]}.tmp <<EOL
#!/bin/bash
    
    echo Proccessing `pwd`: ${files[i]}
    
    # enter the correct folder
	cd ${files[i]}.kallisto
    # run Kallisto
        kallisto quant $Index -o ${files[i]}.output $Parameters $Reads
        cd ..
        echo FINISHED $WorkDir/${files[i]} 
EOL
  done


for (( i=0; i<${#files[@]} ; i+=1 )) ; do
    sed -i "2i\ Reads=\"`pwd`/${files[i]} \"" commands.2.${files[i]}.tmp
done

for (( i=0; i<${#files[@]} ; i+=1 )) ; do
    source commands.2.${files[i]}.tmp
done

# Find all files with abundance.tsv extension and cut the first and $2 column, save it as .cut file

find -name '*abundance.tsv' | xargs -I % sh -c 'cut -f1,4 % | sed "s/^[^|]*|//g" | sed "s/|.*|.*|.*|.*|//g" > %.cut1;'

# remove header
find -name '*abundance.tsv.cut1' | xargs -I % sh -c 'tail -n+2 % > %.cut2;'

#select the top GENCODE isoform according to estimated read counts

find -name '*abundance.tsv.cut1.cut2' | xargs -I % sh -c 'awk '\''{for(i=1;i<=NF;i++) t+=$i; print t"\t"$0; t=0}'\'' % | sort -k2,2 -k1,1nr | awk '\''!a[$2]++'\'' | cut -f2- > %.cleaned;'
find -name '*abundance.tsv.cut1.cut2.cleaned' | xargs -I % sh -c 'tabsep %;'

# download fileMulti2TableMod1.awk

wget https://raw.githubusercontent.com/milospjanic/fileMulti2TableMod1/master/fileMulti2TableMod1.awk

# Find .file.cut files and call fileMulti2TableMod1.awk script to create master table

filescut=$(find -name *.cut1.cut2.cleaned | sort | tr '\n' ' ')  
awk -f fileMulti2TableMod1.awk $(echo $filescut)> mastertable

#clean up mastertable
sed -e 's/\.[0-9]*//g' mastertable | sed 's/_[0-9]*//g' > mastertable.2
mv mastertable.2 mastertable

# add header to mastertable

find -name *.cut1.cut2 | sort | sed 's/.fastq.gz.output.abundance.tsv.cut1.cut2//g' | sed -e 's/.*.fastq.gz.kallisto.//g' | tr '\n' ' ' > header
awk '{$1=" "$1}1' header > header2
cat header2 mastertable > mastertable.2

mv mastertable.2 mastertable

# remove .cut files
find -name *.cut1 | xargs rm
find -name *.cut1.cut2 | xargs rm

# remove fileMulti2TableMod1.awk
rm fileMulti2TableMod1.awk

#remove header
rm header
rm header2


###write R script for ID conversion, needs biomaRt

echo "#!/usr/bin/Rscript
library(biomaRt)
listMarts(host=\"grch37.ensembl.org\")

ensembl = useMart(\"ENSEMBL_MART_ENSEMBL\",dataset=\"hsapiens_gene_ensembl\", host=\"grch37.ensembl.org\")

id_merge_mrna = getBM(attributes=c(\"ensembl_gene_id\",\"refseq_mrna\"),mart=ensembl)
write.table(id_merge_mrna, file=\"id_merge.mrna.txt\", sep = \"\t\", quote =F, col.names=F, row.names=F)

id_merge_ncrna = getBM(attributes=c(\"ensembl_gene_id\",\"refseq_ncrna\"),mart=ensembl)
write.table(id_merge_ncrna, file=\"id_merge.ncrna.txt\", sep = \"\t\", quote =F, col.names=F, row.names=F)

" > script.r

#run R script

chmod 775 script.r
./script.r

awk < id_merge.mrna.txt > id_merge.mrna.short.txt 'NF>1'
awk < id_merge.ncrna.txt > id_merge.ncrna.short.txt 'NF>1'
mv id_merge.mrna.short.txt id_merge.mrna.txt
mv id_merge.ncrna.short.txt id_merge.ncrna.txt
cat id_merge.mrna.txt id_merge.ncrna.txt > id_merge.txt

tabsep id_merge.txt
tabsep mastertable

#Use awk to append gene names

awk 'NR==FNR {h[$1] = $1; h2[$1] = $2; next} {if (h2[$1]) print h2[$1], $0}' id_merge.txt mastertable > mastertable.genename
tabsep mastertable.genename
cut -f1,3- mastertable.genename >mastertable.genename.2
mv mastertable.genename.2 mastertable.genename
head -n1 mastertable > header
sed -i 's/.fastq.gz//g'  header
sed -i 's/.STAR//g' header
cat header mastertable.genename > mastertable.genename.2
mv mastertable.genename.2 mastertable.genename
rm header
tabsep mastertable.genename

#reassign GENCODE counted reads on RefSeq identifiers - solution to ambiguous isoform assignment error

awk '{for(i=1;i<=NF;i++) t+=$i; print t"\t"$0; t=0}' mastertable.genename | sort -k2,2 -k1,1nr | awk '!a[$2]++' | cut -f2- > mastertable.genename.cleaned
mv mastertable.genename.cleaned mastertable.genename
tabsep mastertable.genename

#remove temporary files

rm id_merge.mrna.txt
rm id_merge.ncrna.txt
#rm id_merge.txt
rm script.r


###create  R script

touch script.R

echo "#!/usr/bin/Rscript" > script.R
echo "library(rgsepd)" >>script.R
echo "data<-read.delim(\"mastertable.genename\", header=T, row.names = 1, check.names=F)" >>script.R
echo "meta<-read.delim(\"meta.data\", header=T)" >>script.R

echo "G <- GSEPD_INIT(Output_Folder=\"GSEPD OUTPUT\", finalCounts=data, sampleMeta=meta, COLORS=c(blue=\"#4DA3FF\",black=\"#000000\",gold=\"#FFFF4D\"))" >> script.R

echo "G <- GSEPD_ChangeConditions( G, c(\"A\",\"B\"))" >>script.R
echo "G <- GSEPD_Process( G )" >> script.R

chmod 775 script.R
./script.R
rm script.R

#create R script for DESeq

touch script.deseq.R

echo "#!/usr/bin/Rscript" > script.deseq.R
echo "library(DESeq)" >> script.deseq.R
echo "data<-read.delim(\"mastertable.genename\", header=T, row.names = 1, check.names=F)" >>script.deseq.R
echo "meta<-read.delim(\"meta.data\", header=T)" >>script.deseq.R
###
echo "conds<-factor(meta\$Condition)
sampleTable<-data.frame(sampleName=colnames(data), condition=conds)
countsTable = data
#rownames(countsTable) <- countsTable$Geneid
#countsTable <- countsTable[,-1]
cds <- newCountDataSet( countsTable, conds)
cds <- estimateSizeFactors( cds )
sizeFactors( cds )
head(counts(cds))
head(counts(cds,normalized=TRUE))
cds = estimateDispersions( cds, method=\"blind\", sharingMode=\"fit-only\" )
str( fitInfo(cds) )
plotDispEsts( cds )

res = nbinomTest( cds, \"A\", \"B\" )
head(res)
plotMA(res)
hist(res\$pval, breaks=100, col=\"skyblue\", border=\"slateblue\", main=\"\")
resSig = res[ res\$padj < 0.1, ]
write.csv( res, file=\"Result_Table.csv\" )
write.csv( resSig[ order( resSig\$foldChange, -resSig\$baseMean ), ] , file=\"DownReg_Result_Table.csv\" )
write.csv( resSig[ order( -resSig\$foldChange, -resSig\$baseMean ), ], file=\"UpReg_Result_Table.csv\" )

cdsBlind = estimateDispersions( cds, method=\"blind\" )
vsd = varianceStabilizingTransformation( cdsBlind )
library(\"RColorBrewer\")
library(\"gplots\")

select = order(rowMeans(counts(cds)), decreasing=TRUE)[1:250]
hmcol = colorRampPalette(brewer.pal(9, \"GnBu\"))(100)
heatmap.2(exprs(vsd)[select,], col = hmcol, trace=\"none\", margin=c(10, 6))

select = order(rowMeans(counts(cds)), decreasing=TRUE)[1:500]
heatmap.2(exprs(vsd)[select,], col = hmcol, trace=\"none\", margin=c(10, 6))

select = order(rowMeans(counts(cds)), decreasing=TRUE)[1:1000]
heatmap.2(exprs(vsd)[select,], col = hmcol, trace=\"none\", margin=c(10, 6))

print(plotPCA(vsd, intgroup=c(\"condition\")))

" >> script.deseq.R

chmod 775 script.deseq.R
./script.deseq.R
rm script.deseq.R

duration=$SECONDS
echo "$(($duration / 60)) minutes and $(($duration % 60)) seconds elapsed."
