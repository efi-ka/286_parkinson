# Libraries
library(GenomicRanges)
library(ExomeDepth)
# library(tidyverse)
# library(foreach)
# library(doParallel)
# library(ggplot2)



### CNV coordinates and IDs in hg38
hg38CSVs="/scratch/Reference_sequences/Genomes/hg38/DGV/dgv_cvs_prefixed.bed"
hg38CSVs=read.csv(hg38CSVs,header=TRUE,sep="\t")
Ncsvs=nrow(hg38CSVs)

### Exon coordinates in hg38
# hg38exons="/scratch/Reference_sequences/Genomes/hg38/MANE/1.2/mane_exons.bed"
hg38exons="/scratch/Reference_sequences/Genomes/hg38/MANE/1.2/mane_exons_autosomesonly.bed"
# hg38exons="/scratch/Reference_sequences/Genomes/hg38/MANE/1.2/mane_exons_autosomesonly_testset.bed"
# hg38exons="/scratch/Reference_sequences/Genomes/hg38/MANE/1.2/mane_exons_autosomesonly_testset_noprefix.bed"
hg38exons=read.csv(hg38exons,header=TRUE,sep="\t")
Nexons=nrow(hg38exons)



### CNV identity annotation object
hg38CSVsGR=GRanges	(
						seqnames=hg38CSVs$chromosome,
						ranges=IRanges(hg38CSVs$start,end=hg38CSVs$stop,names=as.character(1:Ncsvs)),
						names=hg38CSVs$name
					)

### Gene identity annotation object
hg38exonsGR=GRanges	(
						seqnames=hg38exons$chromosome,
						ranges=IRanges(hg38exons$start,end=hg38exons$end,names=as.character(1:Nexons)),
						names=hg38exons$name
					)



### Loading BAMs
bamspath="/home/joel/Projects/Parkinson/CNV/input/bam_realpaths.txt"
# bamspath="/home/joel/Projects/Parkinson/CNV/input/testbampaths.txt"
bams=readLines(bamspath)
samplenames=sapply(bams,function(x) paste('X',tail(unlist(strsplit(x,"/")),n=1),sep=''),USE.NAMES=F) # Top kek

### Count reads per every exon interval
message('');message('Counting reads in all bams...')
readcounts=getBamCounts	(
							bed.frame=hg38exons,
							bam.files=bams,
							include.chr=F
						)
readcounts=as(readcounts,'data.frame') # Just to cast it as a df. Perhaps pointless?

# The reference set calculations require a matrix with only counts...
readcounts_mat=as.matrix(readcounts[,5:ncol(readcounts)])

### Loop starts here, probably
for (i in 1:length(readcounts_mat)) {

		# Constructing the reference set
		message('');message('Constructing reference...')

		sample=readcounts_mat[,i]
		ref_counts=readcounts_mat[,-i]

		ref_set=select.reference.set	(
											test.counts=sample,
											reference.counts=ref_counts,
											bin.length=(readcounts$end-readcounts$start)/1000,
											n.bins.reduced=10000
										)

		ref_selected=as.matrix(readcounts[,ref_set$reference.choice,drop=FALSE])
		ref_selected=apply(X=ref_selected,MAR=1,FUN=sum)

		# Testing the thing
		message('');message('Testing...')
		test=new	(
						'ExomeDepth',
						test=sample,
						reference=ref_selected,
						formula='cbind(test,reference)~1'
					)

		message('');message('Calling CNVs...')
		cnvcalls=CallCNVs	(
								x=test,
								transition.probability=10^-4,
								chromosome=readcounts$chromosome,
								start=readcounts$start,
								end=readcounts$end,
								name=readcounts$exon
							)

		# Annotation
		message('');message('Annotation gene names...')
		cnvcalls=AnnotateExtra(x=cnvcalls,
		reference.annotation=hg38exonsGR,
		min.overlap=0.01,
		column.name='hg38.exons')

		message('');message('Annotation CNV IDs...')
		cnvcalls=AnnotateExtra(x=cnvcalls,
		reference.annotation=hg38CSVsGR,
		min.overlap=0.5,
		column.name='DGV_CSVs')

		# Output file
		message('');message('Writing output...')
		sample=substr(samplenames[i],2,9999);sample=substr(sample,1,nchar(sample)-4) # Get rid of prefix and suffix
		sample=paste('output',sample,sep='/');sample=paste(sample,'CNV_calls.tsv',sep='_')
		write.table(cnvcalls@CNV.calls,file=sample,quote=F,sep='\t',row.names=F)
}
# Loop end