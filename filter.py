import sys,argparse

############ Input and such ############
desc    ='''This program filters the ExomeDepth output TSV files by gene symbol. Genes associated with
PD according to HPO are included, others excluded. Sorts on BF (call quality).

python3 filter.py infile outfile
'''
parser  =argparse.ArgumentParser(description=desc,formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('f',metavar='ED',help='Path to ED tsv file.')
parser.add_argument('o',metavar='O',help='Path to output tsv file.')
args    =parser.parse_args()
########################################


with open('../input/hpo_genes_parkinsonism.tsv') as f:
	hpo={line.strip() for line in f}

# Columns:
# start.p    end.p    type    nexons  start   end     chromosome      id      BF      reads.expected  reads.observed  reads.ratio     hg38.exons      DGV_CSVs

keptlines=list()
with open(args.f) as f:
	header=next(f) # save header
	for line in f:
		fields=line.split()
		genes={x.split('_')[0] for x in fields[12].split(',')}
		if hpo&genes:
			keptlines.append(line)
keptlines.sort(key=lambda x: float(x.split()[8]),reverse=True)

if keptlines:
	with open(args.o,'w') as f:
		f.write(header)
		for line in keptlines:
			f.write(line)