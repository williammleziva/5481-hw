import argparse
from operator import indexOf
import numpy as np
import pandas as pd
from skbio import TreeNode as tn
from copy import deepcopy
import random

#python3 NeiSaitou.py -i data/example/example.fna
#Rscript data/hw3-plot-edges.r edges.txt data/example/example-tip-labels.txt bootstrap.txt tree-bootstrap.pdf
#Rscript data/hw3-plot-edges.r edges.txt data/example/example-tip-labels.txt tree.pdf
#Rscript data/hw3-plot-newick.r tree.tre data/example/example-tip-labels.txt tree-newick.pdf


##Define arguement parser
parser = argparse.ArgumentParser(prog='Nei-Saitou neighbor joining', usage ='%(prog)s [-i]')
parser.add_argument('-i', '--input', help = "input fasta file", required=True, type=str)
args = parser.parse_args()

#function to count sequences+ length in file.
def sequenceCount(file):
    with open(file,"r") as fin:
        line_count = sum(1 for line in fin)
    with open(file,"r") as fin:
        fin.readline()
        line_length = len(fin.readline())
    return int(line_count/2),line_length

#function to calculate dissimilarity between sequences, returns dissimilarity
def similarity(str1, str2):
    return (len(str1)-sum(char1 == char2 for char1,char2 in zip(str1,str2)))/(len(str1)-1)

#call seq_count, initialize matrix
seq_count,seq_length = sequenceCount(args.input)
d_matrix = np.zeros((seq_count,seq_count))


names = [] #init lists to hold seq names and indices
indexList = []
indexNum = 0
with open(args.input,"r") as fin1: #open file once, iterate through
    row = 0
    for str1 in fin1:
        if str1[0] == '>': #if name line
            names.append(str1[1:-1]) #add name+index to lists
            indexList.append(indexNum) 
            indexNum+=1
        else:   #if seq line
            with open(args.input,"r") as fin2: #open file second time
                col = 0
                for str2 in fin2: 
                    if str2[0] != '>': #if not name 
                        d_matrix[row,col] = similarity(str1,str2) #calculate similarity, add to matrix
                        #print("c: ",comp,"s1: ",str1[0:20]," s2: ",str2[0:20])
                        col+=1
            row+=1

disSimDF = pd.DataFrame(d_matrix,columns=names,index=names) #load into pandas df with labels          
#print(disSimDF)
#write genetic distances
with open('genetic-distances.txt','w') as genDistFile:
    #remove initial whitespace and convert spaces to tab
    genDistFile.write((disSimDF.to_string()).lstrip().replace('    ','\t').replace('  ','\t'))
print( '( - genetic-distances.txt WRITTEN - )')

def Q_Matrix(d_matrix): #function to caluclate qmatrix form dmatrix, returns matrix locations
    minRow = minCol = 0 #init values
    minVal = float('inf') 
    q_matrix = np.zeros(d_matrix.shape)#init matrix
    for i in range(d_matrix.shape[0]): #iterate through matrix
        for j in range(d_matrix.shape[1]):
            if i != j:
                q_matrix[i][j] = (d_matrix.shape[0]-2)* d_matrix[i][j] - sum(d_matrix[i]) - sum(d_matrix[j]) #calculate q val, add to qmatrix
                if q_matrix[i][j] <= minVal: #save smallest value+ positions
                    minRow=i; minCol=j; minVal=q_matrix[i][j]
    #print("minVal: ",minVal,"  minRow: ",minRow,"  minCol: ",minCol) 
    #print(q_matrix)
    return minRow,minCol

def PairMemberDistanceToNewNode(d_matrix,minRow,minCol): #function to get min pair distance to new node, return edges
    RowEdge = (d_matrix[minRow][minCol]/2.0) + ((sum(d_matrix[minRow])-sum(d_matrix[minCol]))/(2.0*(d_matrix.shape[0]-2)))
    #ColEdge = d_matrix[minRow][minCol] - RowEdge   
    ColEdge = (d_matrix[minCol][minRow]/2.0) + ((sum(d_matrix[minCol])-sum(d_matrix[minRow]))/(2.0*(d_matrix.shape[0]-2)))
    #print("rowEdge: ",RowEdge,"  ColEdge: ",ColEdge)
    return RowEdge,ColEdge

def UpdateD_Matrix(d_matrix,minRow,minCol): #function to return new dmatrix
    newNodeDistRow = [] #init new row
    for t in range(d_matrix.shape[0]): #calculate new row/col vals
        newNodeDistRow.append((d_matrix[minRow][t]+d_matrix[minCol][t]-d_matrix[minRow][minCol])/2.0)
    #prep new row/col to be added
    newNodeDistCol = newNodeDistRow.copy()
    newNodeDistCol.append(0.0) 
    newNodeDistRow = np.asarray(newNodeDistRow)
    newNodeDistCol = np.reshape(newNodeDistCol,(-1,1))
    #add new col and row
    d_matrix = np.vstack([d_matrix,newNodeDistRow])
    d_matrix = np.hstack([d_matrix,newNodeDistCol])
    #remove previous minimum pair row+col
    d_matrix = np.delete(d_matrix,[minRow,minCol],axis=0)
    d_matrix = np.delete(d_matrix,[minRow,minCol],axis=1)

    return d_matrix

def UpdateIndexList(minRow,minCol,indexList,indexNum): #function to keep track of noded index
    #remove old min indexes
    if minRow>minCol: 
       rowindex= indexList.pop(minRow)
       colindex= indexList.pop(minCol)
    else:
       colindex= indexList.pop(minCol)
       rowindex= indexList.pop(minRow)
    #add new index to list
    indexList.append(indexNum)
    #print("list: ",indexList,"  row: ",rowindex,"  col: ",colindex)
    return indexList, rowindex, colindex



def NeighborJoin(d_matrix,indexList,indexNum):
    relationList=[] #init relation loop

    while d_matrix.shape[0] > 2: #main nj algorithm loop

        minRow, minCol = Q_Matrix(d_matrix) #get min pair

        rowEdge,colEdge = PairMemberDistanceToNewNode(d_matrix,minRow,minCol) #get min pair distances to new node

        d_matrix = UpdateD_Matrix(d_matrix,minRow,minCol) #get new dmatrix w/o min pair and with new node

        indexList, rowIndex, colIndex = UpdateIndexList(minRow,minCol,indexList,indexNum) #update list and get index

        #add node to relationship list 
        node1 = [rowIndex,rowEdge,indexNum] #child,dist,parent
        node2 = [colIndex,colEdge,indexNum]
        if rowEdge<=colEdge:
            relationList.append(node1.copy())
            relationList.append(node2.copy())
        else:
            relationList.append(node2.copy())
            relationList.append(node1.copy())
        indexNum+=1#increment node index

    #add last node to relationlist
    indexList, colIndex, rowIndex = UpdateIndexList(0,1,indexList,indexNum)#last node 
    root = [colIndex,d_matrix[0][1],rowIndex] 
    relationList.append(root.copy())

    #flip list
    relationList.reverse()
    #print(relationList)
    return relationList,indexList

relationList,indexList = NeighborJoin(d_matrix,indexList,indexNum)#call neighbor join

#keep a copy without names changed
relationIndexList = deepcopy(relationList)

#replace leaf index with seq names in relation list
for node in relationList:
    if node[0] < len(names):
        node[0] = names[node[0]]
    if node[2] < len(names):
        node[2] = names[node[2]]

#root node for tree structure
tr_names = tn(str(relationList[0][2]))
#attach nodes from relation list to tree
for node in relationList:
    n = tr_names.find(str(node[2]))
    n.append(tn(name = str(node[0]),length = node[1]))

tr_names.write('tree.tre') #write newick file
print( '( - tree.tre WRITTEN - )')

def convertIndex(seq_count,index): #convert a single node index to what is intended
    index = int(index)
    if index < seq_count:
        index+=1
    else: 
        index = seq_count+(2*seq_count-index)-2
        #index+=1
    return str(index) 

def constructIndexTree(seq_count,relationIndexList): #makes tree from relation list with correct indices
    tr_index = tn(convertIndex(seq_count,str(relationIndexList[0][2])))
    #attach nodes fromn relation list to tree
    for node in relationIndexList:
        n = tr_index.find(convertIndex(seq_count,str(node[2])))
        n.append(tn(name = convertIndex(seq_count,str(node[0])),length = node[1],support = convertIndex(seq_count,str(node[2]))))
    return tr_index

tr_index = constructIndexTree(seq_count,relationIndexList) # get tree with index

print(tr_names.ascii_art())
#print(tr_index.ascii_art())

with open("edges.txt","w") as fout: #write to edges.txt file
    for node in tr_index.preorder(): #navigate as preorder
        if node.name != convertIndex(seq_count,str(relationList[0][2])): #skip first because no descendants
            fout.write(node.support+"\t"+node.name+"\t"+str(node.length)+"\n")

def get_tips(tree): #returns array of list of tips for each node 
    leaflist = []
    for node in tree.preorder(): #for each node
        tiplist = []
        for n in node.tips(): #save name
            tiplist.append(n.name)
        tiplist.sort() #sort to better check equality later
        if len(tiplist) > 0:
            leaflist.append([node.name,tiplist.copy()])#dont save leaf nodes
        #print ("node: ",node.name,"  tips: ",tiplist)
    leaflist=np.array(leaflist.copy()) #convert to array
    return leaflist 

mainTips = get_tips(tr_index) # get tip array for real sequences
#print(mainTips)

## bootstrap
def getcolumn(selectedCol,file): #returns list of characters, 1 from each sequence at location
    #column = np.chararray((1,seq_count))
    col = []
    with open(file,"r") as fin: #open sequence file
        for line in fin: #for each sequence, get character at location
            if line[0]!='>':
                col += line[selectedCol] # add character to list
    return col

def constructBootstrapDF(seq_length,seq_count,file): #make dataframe of bootstrap data with seq_count=row, and seq_length=col
    bootstrapSequences = np.chararray((seq_length,seq_count))
    for i in range(seq_length): #for each character in sequence
        selectCol = random.randint(0, seq_length-1) #choose random location
        #print(selectCol)
        column = getcolumn(selectCol,file) #get coluimn at location
        # print column
        for j in range(seq_count): #assign characters from list to sequence
            bootstrapSequences[i,j] = column[j]
    bootstrapDF = pd.DataFrame(bootstrapSequences) #convert to pandas DataFrame
    return bootstrapDF

def colToString(DataFrame, colIndex):   #convert from DataFrame column to string 
    seqstr = DataFrame[colIndex].to_list()
    seqstr=[x.decode('utf-8') for x in seqstr]
    seqstr = ''.join(seqstr)
    return seqstr

def bootstrapDMatrix(seq_count,bootstrapDF): #create dmatrix form bootstrap DataFrame
    d_matrix = np.zeros((seq_count,seq_count)) #init dmatrix
    indexList = []
    for str1Index in range (seq_count): #loop through sequences once
        str1 = colToString(bootstrapDF, str1Index) #get string from DFcolumn
        indexList.append(str1Index) 
        for str2Index in range (seq_count): #loop through sequences again
            str2 = colToString(bootstrapDF, str2Index)
            d_matrix[str1Index,str2Index] = similarity(str1,str2) #calculate similarity, add to matrix
            #print("c: ",comp,"s1: ",str1[0:20]," s2: ",str2[0:20])
    indexNum = seq_count
    return d_matrix,indexList,indexNum

def compareTips(mainTips,bootstrapTips,TipSums): #compare a bootstrap tip array to main keep totals if identical
    equalityList = []
    newTipSums = []
    for node1 in mainTips: #iterate through real seq tips 
        appendval = 0.0
        for node2 in bootstrapTips: #iterate through bootstrap tips
            if node1[1] == node2[1]: #if identical list
                appendval+=1.0
        equalityList.append(appendval) #list of 1 for same or 0 for not 
    for x in range (0, len (equalityList)):  
        newTipSums.append( TipSums[x] + equalityList[x])  #add to sum total
    #print(newTipSums)
    return newTipSums

print( '( - BEGINNING BOOTSTRAP - )')

TipSums=[]
for i in range(len(mainTips)):#fill tip array with zeros to keep Idnetical tip sums 
    TipSums.append(0)

for i in range(100): #100 bootstrap sequences 
    print("Bootstrap Count: ",i)
    bootstrapDF = constructBootstrapDF(seq_length,seq_count,args.input )#get dataframe 

    d_matrix,indexList,indexNum = bootstrapDMatrix(seq_count,bootstrapDF) #get dmatrix

    relationList,indexList = NeighborJoin(d_matrix,indexList,indexNum) #run NJ

    tr_index = constructIndexTree(seq_count,relationList) #build tree
    #print(tr_index.ascii_art())
    bootstrapTips = get_tips(tr_index)#get tip array
    #print(bootstrapTips)
    TipSums = compareTips(mainTips,bootstrapTips,TipSums) #check tip similarity and keep totals

with open("bootstrap.txt","w") as fout: #wirte to bootstrap.txt
    i =0
    for tip in mainTips: # ! Note --> subtracting 0.01 here as 1 in the bootstrap file was not being visualaized in the Rscript
        fout.write((str(tip[0]))+"\t"+str((int(TipSums[i])/100.0)-.01)+"\n") #dividing by bootstrap count (100)
        i+=1
print( '( - bootstrap.txt WRITTEN - )\n\n')

print( '( - COMPLETE - )')