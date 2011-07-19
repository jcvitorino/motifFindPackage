/*
 * OligoAnalisys.cpp
 *
 *  Created on: 19/07/2011
 *      Author: josue
 */

#include "OligoAnalisys.h"

OligoAnalisys::OligoAnalisys() {
	// TODO Auto-generated constructor stub

}

OligoAnalisys::~OligoAnalisys() {
	// TODO Auto-generated destructor stub
}


vector<struct Fe_b> OligoAnalisys::performAnalysis(XNA seq, double n, int len) {
	vector<struct Fe_b> Feb;
	struct  Fe_b Fq;
	int S = seq.chain.size();
	int L = seq.chain[0].length();
	cout << L <<endl;
	//	for (int j = 1; j < 10; ++j) {
		XNAseq s;
		for (int k = 0; k < len; ++k) {
			s.seq = s.seq + 'A';
		}
		int i = 0;
		int cont = pow(4.0, len);
		while (i < cont) {
			Fq.F = computeFrequency(seq, s);
			if (Fq.F.freq != 0) {
				Fq.T = 2 * S * (L - Fq.F.olig.seq.length() + 1);
				Fq.Eocccb = Fq.F.freq * Fq.T;
				Fq.Poccb_eq = probOccBEq(Fq.Eocccb, Fq.T, n);
				Fq.Poccb_great = probOccBGreater(Fq.Eocccb, Fq.T, n);
				Fq.sig = sigValue(Fq.Poccb_great, Fq.F.olig.seq.length());
				cout << Fq.F.olig.seq << ' ' << Fq.T << ' ' << Fq.F.freq << ' '
								<< Fq.Eocccb << ' ' << Fq.Poccb_eq << ' ' << Fq.Poccb_great
								<< ' ' << Fq.sig << endl;
				Feb.push_back(Fq);
			}
			s = nextLeaf(s, len);
			i++;
		}
//	}
	return Feb;
}

double OligoAnalisys::probOccBEq(double Eocc, double T, double n) {
	double TFactorial = fatorial(T);
	double T_minus_n = fatorial((T - n));
	double fatorial_n = fatorial(n);

	double p_Eocc = pow(Eocc, n);
	double p2_Eocc = pow((1 - Eocc), (T - n));

	double partialResul = TFactorial / (T_minus_n * fatorial_n);
	double resul = partialResul * p_Eocc * p2_Eocc;
	return resul;
}

double OligoAnalisys::probOccBGreater(double Eocc, double T, double n) {
	double sum = 0;
	for (int j = n; j <= T; ++j) {
		sum = sum + probOccBEq(Eocc, T, j);
	}
	return (sum);
}

double OligoAnalisys::sigValue(double Poccb, int w) {
	double D = pow(4.0, w) - ((pow(4.0, w) - pow(4.0, (w / 2)))/2);
	double sig = -log10(Poccb * D);
	return sig;
}
