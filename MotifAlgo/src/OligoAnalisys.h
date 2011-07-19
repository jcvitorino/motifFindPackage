/*
 * OligoAnalisys.h
 *
 *  Created on: 19/07/2011
 *      Author: josue
 */

#ifndef OLIGOANALISYS_H_
#define OLIGOANALISYS_H_

#include <vector>
#include "medianString.h"

using namespace std;

struct Fe_b {
	struct Freq F;
	double T;
	double Eocccb;
	double Poccb_eq;
	double Poccb_great;
	double sig;

};

class OligoAnalisys {
public:
	OligoAnalisys();
	virtual ~OligoAnalisys();
	vector<struct Fe_b> performAnalysis(XNA seq, double n, int len);
	double probOccBEq(double Eocc, double T, double n);
	double probOccBGreater(double Eocc, double T, double n);
	double sigValue(double Poccb, int w);
	double fatorial(double x) {
		if (x == 0 || x == 1)
			return 1;
		else
			return (x * fatorial(x - 1));
	}


};

#endif /* OLIGOANALISYS_H_ */
