#include <iostream>
#include <fstream>
#include <set> 
#include <bitset>
#include <string>
#include <time.h>
#include <math.h> 
using namespace std;

#define OVERLAP_THRESHOLD 3

const unsigned int A = 0b001;
const unsigned int T = 0b011;
const unsigned int G = 0b101;
const unsigned int C = 0b111;
const unsigned int X = 0b110;
const unsigned int U = 0b000;

void brint(unsigned int var){
	cout << bitset<sizeof(var)*8>(var) << endl;
}

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
bind takes two 4-byte numbers that represent DNA strands and returns True if the "bind" or False if they will not

This function takes a very conservative and naieve approach to binding, and does not take order of the DNA strand into effect.
For each base match, it will increment and eventually return a counter based on the two bases.
A~T,G~C      => +1
A||T~G||C    => -1
A||T||G||C~X => 1
A||T||G||C~U => 0

*/
/*
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

int results[9] = {-1,1,1,1,-2,1,-1,1,0};
int bindBase(unsigned int b1, unsigned int b2){
	if(!((b1|b2)&1==0)) return results[b1^b2];
	if(b1==U||b2==U) return results[8];
	if(b1==X||b2==X) return results[1];
	cout << "Error -1" << endl;
	throw -1;
}

int bind(unsigned int s1, unsigned int s2, unsigned int numBases, unsigned int offset=0){
        // numBases can't be larger than sizeof(unsigned int)
        if(numBases*3 > sizeof(unsigned int)*8) throw 3;
	int result = 0;
	unsigned int mask = 7;
	for(int i=0;i<numBases;i+=1){
		if(i>=offset){
		unsigned int s1m = s1 & mask;
		unsigned int s2m = s2 & mask;
		int comp = bindBase(s1m>>i*3, s2m>>i*3);
                result += comp;
                mask = mask << 3;
		}
	}
	return result;
}

bool bindStrand(unsigned int s1, unsigned int s2, unsigned int length){
        int counter = 0;
	int result;
        //binding head on
	result = bind(s1, s2, length);
        if(result>0){ // non binding modes don't contribute negatively 
		counter += result;
	}
	//first test 3'->5' shifts, aka right shifts
	unsigned int s2s = s2;
	for(int i=1;i<length;i++){
		s2s = s2s >> 3;
                result = bind(s1, s2s, length-i);
                if(result>0){ // non binding modes don't contribute negatively 
			counter += result;
		}
	}
	//second test 5'->3' shifts, aka left shifts
	s2s = s2;
	for(int i=1;i<length;i++){
		s2s = s2s << 3;
                result = bind(s1, s2s, length-i);
                if(result>0){ // non binding modes don't contribute negatively 
			counter += result;
		}
	}
	//second test 5'->3' shifts, aka left shifts
	/*
	s2s = s2;
	for(int i=1;i<numBases;i++){
		s2s = s2s << i*3;
	}
	*/
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


set<int> generateStrands(int length){
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

int main()
{
  // A=000, T=001, G=010, C=011, X=100
  // in this case leading 2 bits are empty
  int length = 6;
  int strandy1[length] = {A,G,C,T,T,C};
  int strandy2[length] = {T,C,G,T,A,G};
  unsigned int strand1 = getStrand(strandy1, length);//strand1 3'->5'
  unsigned int strand2 = getStrand(strandy2, length);//strand2 5'->3'
  
  bindStrand(strand1, strand2, length);

  set<int> strands; 
  int genLen=8;
  strands = generateStrands(genLen);
  
  
  /*
  long double totalTimeUsed = 0;
  int numRun = 100;
  for(int q=0;q<numRun;q++){
  clock_t start, end;
  long double timeUsed;
  start = clock();
  for(const unsigned int &strand1 : strands){
  	for(const unsigned int &strand2 : strands){
		bindStrand(strand1, strand2, genLen);
	}	
  }
  end = clock();
  timeUsed = ((long double) (end - start)) / CLOCKS_PER_SEC;
  totalTimeUsed += timeUsed;
  cout << "CPU Time Used: " << timeUsed << endl;
  }
  cout << "Average CPU Time Used: " << totalTimeUsed/numRun << endl;
  */

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
		file.write((char*) &strand2, size);
		bool bound = bindStrand(strand1, strand2, genLen);
		file.write((char*) &bound, 1);
		if((i%twentyith)==0){
			end = clock();
  			timeUsed = ((long double) (end - start)) / CLOCKS_PER_SEC;
			cout << "Completed " << i << " comparisons, " << projectedIterations-i << " to go. "<< (i/projectedIterations)*100 << "\% done, taking a total of ";
			cout << timeUsed/60 << " minutes, with " << ((timeUsed)*(projectedIterations-i))/(i*60)<< " predicted to go." << endl;
		}
 		i++;
	}	
  }
  cout << "Done! Closing file ... " << endl << endl;
  file.close();
  end = clock();
  timeUsed = ((long double) (end - start)) / CLOCKS_PER_SEC;
  cout << "Completed calculating and exporting binding possibilities for all strands of length " << genLen << " for a total of  " << pow(pow(4,genLen),2) << " matching calculations in " << timeUsed/60 << " minutes." << endl;
  }else{
	cout << "Could not open file" << endl;
  }
}
