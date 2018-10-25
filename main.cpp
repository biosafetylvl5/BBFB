#include <iostream>
#include <fstream>
#include <set> 
#include <bitset>
#include <string>
#include <time.h>
#include <math.h> 
#include <thread>
using namespace std;

#define OVERLAP_THRESHOLD 3 // Threshold for what counts as "Binding"
#define OVERLAP_SKIP 4 // Number of basepairs to stop comparing at for frameshifts - eg. 4 won't compare any overlaps for oligos of length 5

const unsigned int A = 0b001;
const unsigned int T = 0b011;
const unsigned int G = 0b101;
const unsigned int C = 0b111;
const unsigned int X = 0b110;
const unsigned int U = 0b000;

// Prints binary representation of arg1
void brint(unsigned int var){
    cout << bitset<sizeof(var)*8>(var) << endl;
}

// Prints out strand variables for humans
//     arg1: strand, arg2: length of strand
//     arg3: print three bits along with charecter representation
//     arg4: offset for printing, eg
//           unsigned int strand1 = getStrand({A,A}, 2);
//           printStrand(strand1, 2, false, 0);
//           printStrand(strand1, 2, false, 1);
//           would print:
//           AA
//            AA                 (with this being the offset one)
void printStrand(unsigned int strand, int length, bool printBin=false, int offset=0){
    if(offset>=0){
        int k = printBin?5:1;
        for(int j=0;j<offset*k;j++) cout << " ";
    }
    for(int i=0;i<length;i++){
        unsigned int base = (strand & (7<<i*3))>>i*3;
        if(base==A) cout << "A";
        if(base==T) cout << "T";
        if(base==G) cout << "G";
        if(base==C) cout << "C";
        if(base==U) cout << "U";
        if(base==X) cout << "X";
        if(printBin) cout << bitset<3>(base) << " ";
    }
    cout << endl;
}

/*
bind takes two 4-byte numbers and returns their binding "score" 

    This function takes a very conservative and naieve approach to binding, 
    and does not take order of the DNA strand into effect.
    For each base match, it will increment and eventually return a counter 
    based on the two bases, the results are defined in the results array:

    results[]:
    000 - Same Base
    001 - Match!
    010 - Match!
    011 - Match!
    100 - Bad Mismatch
    101 - Match!
    110 - Normal Mismatch
    111 - Match!
    1000 - U match
*/

#define MATCH 1
#define NMATCH -1
#define BMATCH -2
#define UMATCH 0


float results[9] = {NMATCH,MATCH,MATCH,MATCH,BMATCH,MATCH,NMATCH,MATCH,UMATCH};
float bindBaseSlow(const unsigned int b1, const unsigned int b2){
    if((b1&b2)&1) return results[b1^b2];
    if(b2==U||b1==U) return UMATCH;
    if(b2==X||b1==X) return MATCH;
    return -1;
}

float baseArray[8][8];
float twoBaseArray[0b111111][0b111111];

void generateTwoBaseArray(){
    for(unsigned int i=0;i<0b1000000;i++){
        for(unsigned int j=0;i<0b1000000;i++){
	    twoBaseArray[i][j] = bindBaseSlow(i,j);
        }
    }

}
void generateBaseArray(){
    for(unsigned int i=0;i<0b1000;i++){
        for(unsigned int j=0;i<0b1000;i++){
	    baseArray[i][j] = bindBaseSlow(i,j);
        }
    }
    generateTwoBaseArray();

}

inline int bind(const unsigned int s1, const unsigned int s2, const unsigned int numBases, const unsigned int offset=0){
    int result = 0;
    unsigned int mask = 7;
    for(int i=offset;i<numBases;i++){
        result += baseArray[(s1&mask)>>i*3][(s2&mask)>>i*3];
        mask = mask << 3;
    }
    return result;
}
inline bool bindStrand(const unsigned int s1, const unsigned int s2, const unsigned int length){
    int counter = 0;
    // binding is compared in three seperate sections to lessen the variable size needed for s1, s2 (data lost by bit shifts isn't a problem this way)

    // zed test: binding head on
    int result = bind(s1, s2, length);
    if(result>0){ // non binding modes don't contribute negatively 
        counter += result;
    }
    //first test: 3'->5' shifts, aka right shifts
    unsigned int s2s = s2 >> 3*OVERLAP_SKIP;
    for(int i=1+OVERLAP_SKIP;i<length;i++){
        s2s = s2s >> 3;
        result = bind(s1, s2s, length-i);
        if(result>0){ // non binding modes don't contribute negatively 
            counter += result;
        }
    }
    //second test: 5'->3' shifts, aka left shifts
    s2s = s2 >> 3*OVERLAP_SKIP;
    for(int i=1+OVERLAP_SKIP;i<length;i++){
        s2s = s2s << 3;
        result = bind(s1, s2s, length-i);
        if(result>0){ // non binding modes don't contribute negatively 
            counter += result;
        }
    }
    return counter>=OVERLAP_THRESHOLD;
}

unsigned int getStrand(int bases[], int length){
    unsigned int strand = 0;
    for(int i=0;i<length;i++){
        unsigned int base = bases[i];
        base = base << i*3;
        strand = strand | base;
    }
    return strand;
}

// getComp returns the DNA complement of a strand, eg. A~T, G~C.
// Note: this does *invert* the direction of the strand, 
//       returning a strand that would bind the given strand
unsigned int getComp(unsigned int strand, int length){
    unsigned int newStrand = 0;
    for(int i=0;i<length;i++){
        unsigned int base = (strand & (7<<i*3))>>i*3;
        if(base==A) newStrand = newStrand | T<<(length-i-1)*3;
        if(base==T) newStrand = newStrand | A<<(length-i-1)*3;
        if(base==G) newStrand = newStrand | C<<(length-i-1)*3;
        if(base==C) newStrand = newStrand | G<<(length-i-1)*3;
        if(base==U) newStrand = newStrand | U<<(length-i-1)*3;
        if(base==X) newStrand = newStrand | X<<(length-i-1)*3;
    }
    return newStrand;
}

// generateStrands take in a length Arg1 and returns a set of all strands of that length
// with base pairs A,G,T,C but not X or U
set<int> generateStrands(int length){
    if(length*3 > sizeof(unsigned int)*8){
        cout << "length to large, got " << length << " for length but unsigned int only ";
        cout << sizeof(unsigned int) << "bytes large!" << endl;
        throw -1;
    }
    if(length==1) return {A,T,G,C};
    set<int> shorterStrands = generateStrands(length-1);
    set<int> strands;
    for(const unsigned int &strand : shorterStrands){
        unsigned int newStrand = strand<<3;
        strands.insert(newStrand|A);
        strands.insert(newStrand|G);
        strands.insert(newStrand|T);
        strands.insert(newStrand|C);
    }
    return strands;
}

void generateRelations(const unsigned int genLen){
  set<int> strands = generateStrands(genLen);
  ofstream file; 
  file.open("binding.bin", ios::out|ios::binary|ios::trunc);
  clock_t start, end;
  long double timeUsed;
  start = clock();
  if (file.is_open()){
  unsigned long i = 0;
  unsigned long projectedIterations = pow(pow(4,genLen), 2);
  unsigned long twentyith = projectedIterations/20;
  int size = sizeof(unsigned int);
  for(const unsigned int &strand1 : strands){
    file.write((char*) &strand1, size);
      for(const unsigned int &strand2 : strands){
        bool bound = bindStrand(strand1, strand2, genLen);
        if(bound) file.write((char*) &strand2, size);
        if((i%twentyith)==0){
            end = clock();
              timeUsed = ((long double) (end - start)) / CLOCKS_PER_SEC;
            cout << "Completed " << i << " comparisons, " << projectedIterations-i;
            cout << " to go. "<< (i/projectedIterations)*100 << "\% done, taking a total of ";
            cout << timeUsed/60 << " minutes, with " << ((timeUsed)*(projectedIterations-i))/(i*60);
            cout << " predicted to go." << endl;
        }
         i++;
    }    
  }
  cout << "Done! Closing file ... " << endl << endl;
  file.close();
  end = clock();
  timeUsed = ((long double) (end - start)) / CLOCKS_PER_SEC;
  cout << "Completed calculating and exporting binding possibilities for all strands of length ";
  cout << genLen << " for a total of  " << pow(pow(4,genLen),2) << " matching calculations in ";
  cout << timeUsed/60 << " minutes." << endl;
  }else{
    cout << "Could not open file" << endl;
  }
}

// testSpeed runs bindStrand for all oligos of length Arg1 Arg2 number of times
// and reports the individual and average time taken for each run
// eg. testSpeed(5,10) would compare all oligos of length 5 10 different times
void testSpeed(const unsigned int genLen, const unsigned int numRun){
  set<int> strands = generateStrands(genLen);
  long double totalTimeUsed = 0;
  clock_t start, end;
  long double timeUsed;
  for(int q=0;q<numRun;q++){
      start = clock();
      for(const unsigned int &strand1 : strands){
        for(const unsigned int &strand2 : strands){
          bindStrand(strand1, strand2, genLen);
        }    
      }
      end = clock();
      timeUsed = ((long double) (end - start)) / CLOCKS_PER_SEC;
      totalTimeUsed += timeUsed;
      cout << genLen << "bp Test CPU Time Used Run "<< q <<": " << timeUsed << endl;
  }
  cout << genLen << "bp Average CPU Time Used: " << totalTimeUsed/numRun << endl;
}

int main()
{
  int length = 6;
  generateBaseArray();
  //for human testing, we can create strand vars like this:
  int strandy1[length] = {A,G,C,T,T,C};
  int strandy2[length] = {T,C,G,T,A,G};
  unsigned int strand1 = getStrand(strandy1, length); // strand1 3'->5'
  unsigned int strand2 = getStrand(strandy2, length); // strand2 5'->3'
  bindStrand(strand1, strand2, length); // to test binding

  testSpeed(4, 100);
}

