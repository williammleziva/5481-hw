import argparse
from operator import indexOf
from pathlib import Path

parser = argparse.ArgumentParser(prog='AAconvert', #create arg variables
                                    usage ='%(prog)s [-s]')
parser.add_argument('-s', '--sequence', help = "path to sequence (.fasta)", required=True, type=str)
parser.add_argument('-o', '--outputfile', help = "path to the output file (.txt)", required=True)
args = parser.parse_args()

lines = Path(args.sequence).read_text().splitlines() #read arg lines into local strings
title = lines[0]
sequence = lines[1]

def translate_dna(seq):  
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',    
    }
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i : i+3]
            protein += table[codon]
            
    return protein

seqprotein = translate_dna(sequence)
#print(seqprotein)

with open(args.outputfile,"w") as fout: 
    fout.write(title+' Amino Acids\n')
    fout.write(seqprotein+'\n')