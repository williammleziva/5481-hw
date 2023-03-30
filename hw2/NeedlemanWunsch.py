import argparse
from operator import indexOf
from pathlib import Path
import numpy as np
import sys

#python3 NeedlemanWunsch.py -q fake_seq1.fasta -r fake_seq2.fasta -o fake_output.txt
#python3 NeedlemanWunsch.py -q hw2_files/pfizer_mrna.fna -r hw2_files/sars_spike_protein.fna -o pfizer_sars_nogap.txt
#python3 NeedlemanWunsch.py -q sars_spike_protein.aa -r pfizer_mrna.aa -o SarsPfizerAAalignment.txt

GAP_PENALTY = -2 #penalties har coded
MISMATCH_PENALTY = -1
MATCH_SCORE = 1

#SECTION 1 - argument parser
parser = argparse.ArgumentParser(prog='NeedlemanWunsch', #create arg variables
                                    usage ='%(prog)s [-q] [-r] [-o]')
parser.add_argument('-q', '--sequence1', help = "path to first sequence (.fasta)", required=True, type=str)
parser.add_argument('-r', '--sequence2', help = "path to second sequence (.fasta)", required=True)
parser.add_argument('-o', '--outputfile', help = "path to the output file (.txt)", required=True)
parser.add_argument('-g', '--gappenalty', action='store_true', help= "penalize start/end gap, default off")
            # if -g included -> start&end gaps penalized
args = parser.parse_args()
###print(args)

#SECTION 2 - read in sequences
lines = Path(args.sequence1).read_text().splitlines() #read arg lines into local strings
title1 = lines[0]
seq1 = lines[1]
len1 = len(seq1)
###print(seq1,len1)

lines = Path(args.sequence2).read_text().splitlines()
title2 = lines[0]
seq2 = lines[1]
len2 = len(seq2)
###print(seq2,len2)

#SECTION 3 - Initializing Matrix
Matrix = np.zeros((len1+1,len2+1,2)) #x,y,0 holds scores, x,y,1 holds directions
if (args.gappenalty): # if start and end gap penalties 
    for i in range(len1+1):
        Matrix[i,0,0] = i * GAP_PENALTY #init gap penalties
    for i in range(len2+1):
        Matrix[0,i,0] = i * GAP_PENALTY 
Matrix[:,0,1] = 10 # init directions : 10=left, 100=above
Matrix[0,:,1] = 100

###print(Matrix[:,:,0],'\n', Matrix[:,:,1])

#SECTION 4 - Fill Matrix
for i in range(len1):
    for j in range (len2): #loop through matrix
        if seq1[i] == seq2[j]: #if match 
            diag = Matrix[i,j,0] + MATCH_SCORE 
        else: #not match
            diag = Matrix[i,j,0] + MISMATCH_PENALTY #get all scores to feed target box(i+1,j+1)
        above = Matrix[i+1,j,0] + GAP_PENALTY 
        left = Matrix[i,j+1,0] + GAP_PENALTY
        best = Matrix[i+1,j+1,0] = max(diag,above,left) #get best, save to matrix
        if (best == diag): Matrix[i+1,j+1,1] += 1  #save direction: 1=diagonal, 10=left, 100=abvove
        if (best == left): Matrix[i+1,j+1,1] += 10 
        if (best == above): Matrix[i+1,j+1,1] += 100 

###print(Matrix[:,:,0],'\n', Matrix[:,:,1])

#SECTION 5 - Find Path
i = 0; j =0 #init
gap_seq1 = []; gap_seq2 = []
alignment_score = 0

if (args.gappenalty): # if start and end gap penalties 
    i = len1; j = len2 #start at bottom right corner
    alignment_score = Matrix[i,j,0] # save score
    while i>0 or j>0: #end at 0,0 corner
        dir = Matrix[i,j,1] # traceback path through matrix, filling lists
        if dir%10 == 1:
            gap_seq1+=seq1[i-1]; gap_seq2+=seq2[j-1]
            i-=1;j-=1
        elif dir%100 > 1: 
            gap_seq1+=seq1[i-1]; gap_seq2+='_'
            i-=1
        elif dir > 11: 
            gap_seq1+='_'; gap_seq2+=seq2[j-1]
            j-=1
else: #no start and end gap penalties 
    bestj = max(Matrix[-1,:,0]) #start at largest endge
    besti = max(Matrix[:,-1,0])
    if bestj >= besti:
        templist = Matrix[-1,:,0].tolist() # get last occurance index for largest at edge
        templist.reverse()
        j = len2 - indexOf(templist,bestj) 
        i = len1
        alignment_score = Matrix[i,j+1,0] #save score
        for k in range (len2,j,-1): #fill lists for 'skipped' edge section
            gap_seq1+='_'
            gap_seq2+=seq2[k-1]
    else:
        templist = Matrix[:,-1,0].tolist() #Same^^ for other edge
        templist.reverse()
        i = len1 - indexOf(templist,besti)
        j = len2
        alignment_score = Matrix[i,j,0]
        for k in range (len1,i,-1):
            gap_seq2+='_'
            gap_seq1+=seq1[k-1]
    ###print(i,besti,j,bestj)

    while i>0 and j>0: #stop at first edge
        ###print('i:',i,'j:',j)
        dir = Matrix[i,j,1] # traceback path through matrix, filling lists
        if dir%10 == 1:
            gap_seq1+=seq1[i-1]; gap_seq2+=seq2[j-1]
            i-=1;j-=1
        elif dir%100 > 1: 
            gap_seq1+=seq1[i-1]; gap_seq2+='_'
            i-=1
        elif dir > 11: 
            gap_seq1+='_'; gap_seq2+=seq2[j-1]
            j-=1
    while j>0: #finish top/side edge for 'skipped' filling lists
        gap_seq1+='_'
        gap_seq2+=seq2[j-1]
        j-=1
    while i>0:
        gap_seq2+='_'
        gap_seq1+=seq1[i-1]
        i-=1

#SECTION 6 - Prepare data
gap_seq1.reverse() #flip 
gap_seq2.reverse()
###print (gap_seq1, '\n', gap_seq2)

matchlist=[]
for k in range(len(gap_seq1)): #fill list for matching illustration
    if(gap_seq1[k]==gap_seq2[k]): matchlist+='|'
    elif(gap_seq1[k]== '_' or gap_seq2[k] == '_'): matchlist+=' '
    else: matchlist+='x'
###print(matchlist)

#convert to string
align1 = "".join(gap_seq1)
matches = "".join(matchlist)
align2 = "".join(gap_seq2)

#SECTION 7 - write data as specified
with open(args.outputfile,"w") as fout: 
    fout.write(str(alignment_score)+'\n')
    fout.write(title1+'\n')
    fout.write(align1+'\n')
    fout.write(matches+'\n')
    fout.write(align2+'\n')
    fout.write(title2+'\n')

