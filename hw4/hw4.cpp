#define _GLIBCXX_USE_CXX11_ABI 0/1

#include <iostream>
#include <fstream>
#include <string>


using namespace std;

int main () {
    string line; 
    ifstream seqfile ("Homework4-seqs-with-primers.fna");

    int seqCount, seqLength = 0; 
    if (seqfile.is_open()){
        while ( getline (seqfile,line) )
        {
        if (line[0]  != '>'){
            seqLength = line.size();
            seqCount++;
        }
        }
    }
    seqfile.close();

    int countArray[seqLength][5] = {0};

    for( int i = 0; i<seqLength; i++){
        if (i%100 == 0) cout << "\n" << i;
        cout << "." << std::flush;
        ifstream seqfile ("Homework4-seqs-with-primers.fna");
        while ( getline (seqfile,line) )
        {
            if (line[i] == 'A'){
                countArray[i][0]+=1;
            }
            else if (line[i] == 'C'){
                countArray[i][1]+=1;
            }
            else if (line[i] == 'G'){
                countArray[i][2]+=1;
            }
            else if (line[i] == 'T'){
                countArray[i][3]+=1;
            }
            else if (line[i] == '-'){
                countArray[i][4]+=1;
            }
        }
        seqfile.close();
    }

    ofstream sol1 ("solution-problem-1.txt");
    for( int i = 0; i<seqLength; i++){
        cout << "A" << countArray[i][0] << "  C" << countArray[i][1] << "  G" << countArray[i][2] << "  T" << countArray[i][3] << "  -" << countArray[i][4] << "\n";
        float k = max(max(countArray[i][0],countArray[i][1]),max(countArray[i][2],countArray[i][3]))/seqLength;
        sol1.write(k);
    }

    return 0;
}