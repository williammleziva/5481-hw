#William Mleziva - mlezi006@umn.edu
#CSCI 5481 HW 1
import sys
import csv

def main():
    #Get filenames
    args = sys.argv[1:] 
    inputfname = args[0];
    outputfname = args[1];
    #check file types
    if(inputfname[-4:]!=".fna"):
        print("Input file must be '.fna'\n")
        exit()
    if(outputfname[-4:]!=".csv"):
        print("Output file must be '.csv'\n")
        exit()
    #create empty lists
    codonlist = [] #holds codon names
    codoncount= [] #holds occurance count for codon in codonlist[] with same index 

    ###open file, iterate through lines
    with open(inputfname) as fin:
        for line in fin:
            
            if (line[0] != '>'): #if line isnt title/header
                sequence = line[:-1]  #remove '\n'

                #remove extra character from sequence end
                extra_char_cnt = (len(sequence))%3  
                if extra_char_cnt!=0:
                    sequence = sequence[:-1*extra_char_cnt] 
                
                #iterate through sequence in sets of 3 chars
                for i in range(0,len(sequence),3):
                    codon = sequence[i : i+3]
                    #if in list add 1 to count
                    if codon in codonlist:                      
                        codoncount[codonlist.index(codon)]+=1
                    #if not in list add to list 
                    else:
                        codonlist.append(codon)                 
                        codoncount.append(1)

    ####write lists to csv
    #initialize writer/open file
    with open(outputfname,'w') as fout:
        writer =csv.writer(fout)
        #iterate through lists, writing rows of [codon,count]
        for i in range(len(codonlist)):
            row = [codonlist[i]] + [codoncount[i]]
            writer.writerow(row)

    ###print lists if uncommented
    #print(codonlist)
    #print(codoncount)

#calls main on file execution
main();