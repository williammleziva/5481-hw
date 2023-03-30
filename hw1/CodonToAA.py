import sys
import csv



def main():
    #Get names and check for file types
    print("starting")
    args = sys.argv[1:] 
    inputfname = args[0];
    outputfname = args[1];

    codonlist = []
    codoncount= []

    if(inputfname[-4:]!=".csv"):                                #read in csv to list
        print("Input file must be '.csv'\n")
    if(outputfname[-4:]!=".csv"):
        print("Output file must be '.csv'\n")

    with open(inputfname, 'r') as fin:
        f = csv.DictReader(fin, fieldnames=['codon','count'])
        for col in f:
            codonlist.append(col['codon'])
            codoncount.append(int(col['count']))

    print(codonlist)
    print(codoncount)

        #convert to AA
    AAtable = {
        'ATA':'Ile', 'ATC':'Ile', 'ATT':'Ile', 'ATG':'Met',
        'ACA':'Thr', 'ACC':'Thr', 'ACG':'Thr', 'ACT':'Thr',
        'AAC':'Asn', 'AAT':'Asn', 'AAA':'Lys', 'AAG':'Lys',
        'AGC':'Ser', 'AGT':'Ser', 'AGA':'Arg', 'AGG':'Arg',
        'CTA':'Leu', 'CTC':'Leu', 'CTG':'Leu', 'CTT':'Leu',
        'CCA':'Pro', 'CCC':'Pro', 'CCG':'Pro', 'CCT':'Pro',
        'CAC':'His', 'CAT':'His', 'CAA':'Gln', 'CAG':'Gln',
        'CGA':'Arg', 'CGC':'Arg', 'CGG':'Arg', 'CGT':'Arg',
        'GTA':'Val', 'GTC':'Val', 'GTG':'Val', 'GTT':'Val',
        'GCA':'Ala', 'GCC':'Ala', 'GCG':'Ala', 'GCT':'Ala',
        'GAC':'Asp', 'GAT':'Asp', 'GAA':'Glu', 'GAG':'Glu',
        'GGA':'Gly', 'GGC':'Gly', 'GGG':'Gly', 'GGT':'Gly',
        'TCA':'Ser', 'TCC':'Ser', 'TCG':'Ser', 'TCT':'Ser',
        'TTC':'Phe', 'TTT':'Phe', 'TTA':'Leu', 'TTG':'Leu',
        'TAC':'Tyr', 'TAT':'Tyr', 'TAA':'Stp', 'TAG':'Stp',
        'TGC':'Cys', 'TGT':'Cys', 'TGA':'Stp', 'TGG':'Trp',    
    }
    
    AAlist = []
    AAcount= []

    for i in range(len(codonlist)):
        AA = AAtable[codonlist[i]]
        if AA in AAlist:                                        #if in list add codoncount
            AAcount[AAlist.index(AA)]+=codoncount[i]
        else:
            AAlist.append(AA)                                   #if not in list add to list and codoncount
            AAcount.append(codoncount[i])

    print(AAlist)
    print(AAcount)


    with open(outputfname,'w') as fout:
        writer =csv.writer(fout)
        for i in range(len(AAlist)):
            row = [AAlist[i]] + [AAcount[i]]
            writer.writerow(row)

main();