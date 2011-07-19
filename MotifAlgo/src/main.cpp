//============================================================================
// Name        : median.cpp
// Author      : josue
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <stdlib.h>
#include <vector>
#include "medianString.h"
#include "OligoAnalisys.h"
using namespace std;

XNA processSet(XNA seq);

int main() {

	//char* finFSA = "cr25s.fsa";
	//char* finFSA = "example.fasta";
	char* finFSA = "oligo.fasta";
	XNA s1;
	OligoAnalisys olig;
	s1.readFasta(finFSA);

	s1 = processSet(s1);
	olig.performAnalysis(s1, 3000, 6);
	return 0;
}

XNA processSet(XNA seq) {
	int numberOfSequences = seq.chain.size();
	int length = 0;
	int lengthAux = 0;
	length = seq.chain[0].length();

	for (int i = 1; i < numberOfSequences; i++) {
		lengthAux = seq.chain[i].length();
		if (lengthAux < length) {
			length = lengthAux;
		}
	}

	for (int j = 0; j < numberOfSequences; j++) {
		seq.chain[j].seq.erase(seq.chain[j].seq.begin() + length,
				seq.chain[j].seq.end());
	}
	return seq;
}
