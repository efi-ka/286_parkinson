#!/bin/python3
import sys
def process_files(genes_file, input_file):
    genes = set()

    try:
        with open(genes_file, 'r') as genes_f:
            for lines in genes_f:
                lines = lines.rstrip()
                genes.add(lines)
    except FileNotFoundError:
        print(f"Error: The file {genes_file} was not found")
        sys.exit(1)        #Exits the program if the file is not found


    with open(input_file, 'r') as input_f:
        for line in input_f:
            line = line.rstrip()
            if line.startswith('#'):
                print(line)
            else:
                fields = line.split('|')
                gen = fields[5]
                if gen in genes:
                    print(line)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 scriptforgenesfun.py <input_file>")
        sys.exit(1)                        #Exits the script if the arguments are incorrect

    genes_file = '/home/efthymia/snakemake_park/parkgenelist.txt'    #Path to the txt file with the PD genes
    input_file = sys.argv[1]

    process_files(genes_file, input_file)
