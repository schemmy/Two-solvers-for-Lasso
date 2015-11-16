#ifndef PROBLEMHELPER_H_
#define PROBLEMHELPER_H_

#include "DWTtrans.h"


void computeObjectiveValue(std::vector<double> &b, std::vector<double> &x, std::vector<double> &fx,
                           unsigned int &n, double &lambda, double &obj, vector<unsigned int> &indexKeep);
void computeGradientFx(std::vector<double> &b, std::vector<double> &x, std::vector<double> &fx,
                       unsigned int &n, std::vector<double> &gradient, vector<unsigned int> &indexKeep);
void computeGradientFx_fast(std::vector<double> &ATb, std::vector<double> &x, unsigned int &n, std::vector<double> &gradient);
void AxMinusb(std::vector<double> &b, std::vector<double> &x,
              std::vector<double> &sol, unsigned int &n, vector<unsigned int> &indexKeep);


void computeObjectiveValue(std::vector<double> &b, std::vector<double> &x, std::vector<double> &fx,
                           unsigned int &n, double &lambda, double &obj, vector<unsigned int> &indexKeep) {

	double normL2f = cblas_dnrm2(n, &fx[0], 1);
	double normL1g = cblas_dasum(n, &x[0], 1);

	obj = normL2f * normL2f + lambda * normL1g;

}

void computeGradientFx(std::vector<double> &b, std::vector<double> &x, std::vector<double> &fx,
                       unsigned int &n, std::vector<double> &gradient, vector<unsigned int> &indexKeep) {

	forwardDWT(fx, gradient, n, indexKeep);
	cblas_dscal(n, 2.0, &gradient[0], 1);

}

void computeGradientFx_fast(std::vector<double> &ATb, std::vector<double> &x, std::vector<double> &fx,
                       unsigned int &n, std::vector<double> &gradient, vector<unsigned int> &indexKeep) {
	
	for (unsigned int i = 0; i < n; i++)
		gradient[i] = 2.0 * (x[i] - ATb[i]);

}


void AxMinusb(std::vector<double> &b, std::vector<double> &x,
              std::vector<double> &sol, unsigned int &n, vector<unsigned int> &indexKeep) {

	inverseDWT(x, sol, n, indexKeep);
	cblas_daxpy(n, -1.0, &b[0], 1, &sol[0], 1);

}


#endif