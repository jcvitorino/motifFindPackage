/*
 * medianString.h
 *
 *  Created on: 19/07/2011
 *      Author: josue
 */

#ifndef MEDIANSTRING_H_
#define MEDIANSTRING_H_
#include <string>
#include <cmath>
#include "XNA.h"

using namespace std;

struct ret {
	XNAseq s;
	int i;
};

struct Freq{
	XNAseq olig;
	double freq;

};


int hammingDistance(string s, string t);
int totalDistance(XNAseq prefix, XNA seq, int sizeP);
struct ret bypass(XNAseq a, int i);
struct ret nextVertex(XNAseq a, int i, int L);
XNAseq medianString(XNA seq, int L, int mis);
XNAseq nextLeaf(XNAseq a, int L);
void allLeaves(int L);
struct Freq computeFrequency(XNA setSeq, XNAseq seq);

#endif /* MEDIANSTRING_H_ */
