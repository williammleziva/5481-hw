# Needleman-Wunsch Alignment Agorithm 
## Written in python 3.7
---
### To view help info and usage, Run 

    python3 NeedlemanWunsch.py -h 


### Example usage with gap penalties for the start and end: 

    python3 NeedlemanWunsch.py -q seq1.fasta -r seq2.fasta -o output.txt -g 

### Example usage with no gap penalties for the start and end: 

    python3 NeedlemanWunsch.py -q seq1.fasta -r seq2.fasta -o output.txt 