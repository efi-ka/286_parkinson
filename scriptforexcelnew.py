#!/usr/bin/python

import sys
import pandas as pd



def process_file(input_file, output_file):

    # Lists for columns
    data = {
        'Chr': [], 'Position': [], 'Ref allele': [], 'Alt allele': [], 'Gene': [], 'Gene_name': [],
        'Existing_variation': [], 'Feature': [], 'Feature_type': [], 'Consequence': [], 'cDNA position': [],
        'CDS_position': [], 'Protein_position': [], 'Amino_acids': [], 'Codons': [], 'Impact': [],
        'HGVSc': [], 'HGVSp': [], '1000 Genomes Phase 3 (EUR)': [], 'gnomAD_genomes_NFE': [],
        'gnomAD_exomes_NFE': [], 'Exon': [], 'Intron': [], 'SIFT': [], 'Polyphen': [], 'LRT_score_prediction': [],
        'MutationTaster_score_prediction': [], 'MutationAssessor_score_prediction': [], 'FATHMM_score_prediction': [],
        'PROVEAN_score_prediction': [], 'MetaLR': [], 'REVEL_score': [], 'CADD_phred': [], 'DANN': [],
        'GERP': [], 'SpliceAI': [], 'ClinVar clinical significance': [], 'ClinVar trait/disease': [],
        'GT': [], 'GQ': [], 'DP': []
    }



    with open(input_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue  # Skip metadata lines



            l = line.rstrip()
            s = l.split('\t')   # Split the record line
            inf = l.split(';')  # Info column is semicolon-separated
            ss = l.split('|')   # Split the record line
            dp = s[9].split(':')[3]
            nr = s[9].split(':')[1]


            # Populate lists

            data['Chr'].append(s[0])
            data['Position'].append(s[1])
            data['Ref allele'].append(s[3])
            data['Alt allele'].append(s[4])
            data['Gene'].append(ss[4])
            data['Gene_name'].append(ss[5])
            data['Existing_variation'].append(ss[3])
            data['Feature'].append(ss[9])
            data['Feature_type'].append(ss[8])
            data['Consequence'].append(ss[1])
            data['Impact'].append(ss[2])
            data['cDNA position'].append(ss[11])
            data['CDS_position'].append(ss[12])
            data['Protein_position'].append(ss[16])
            data['Amino_acids'].append(ss[7])
            data['Codons'].append(ss[6])
            data['HGVSc'].append(ss[10])
            data['HGVSp'].append(ss[15])
            data['1000 Genomes Phase 3 (EUR)'].append(ss[37])
            data['gnomAD_genomes_NFE'].append(ss[44])
            data['gnomAD_exomes_NFE'].append(ss[41])
            data['Exon'].append(ss[13])
            data['Intron'].append(ss[14])
            data['SIFT'].append(ss[34])
            data['Polyphen'].append(ss[35])
            data['LRT_score_prediction'].append(ss[52] + '|' + ss[53] if ss[52] else '')
            data['MutationTaster_score_prediction'].append(ss[54] + '|' + ss[55] if ss[54] else '')
            data['MutationAssessor_score_prediction'].append(ss[56] + '|' + ss[57] if ss[56] else '')
            data['FATHMM_score_prediction'].append(ss[58] + '|' + ss[59] if ss[58] else '')
            data['PROVEAN_score_prediction'].append(ss[60] + '|' + ss[61] if ss[60] else '')
            data['MetaLR'].append(ss[65])
            data['REVEL_score'].append(ss[70])
            data['CADD_phred'].append(ss[71])
            data['DANN'].append(ss[72])
            data['GERP'].append(ss[73])
            data['SpliceAI'].append(inf[-1])
            data['ClinVar clinical significance'].append(ss[74])
            data['ClinVar trait/disease'].append(ss[75])
            data['GT'].append(s[9].split(':')[0])  # Genotype
            data['GQ'].append(s[9].split(':')[6])  # Genotype quality
            data['DP'].append(dp + '[' + nr + ']')  # Read depth


    # Convert the data dictionary to a pandas DataFrame
    df = pd.DataFrame(data)


    # Write to Excel
    df.to_excel(output_file, index=False)



if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python3 scriptforexcelnew.py <input_file> <output_file>")
        sys.exit(1)



    input_file = sys.argv[1]
    output_file = sys.argv[2]


    process_file(input_file, output_file)
