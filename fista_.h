#ifndef FISTA_H_
#define FISTA_H_

#include <cmath>
using namespace std;
#include "DWTtrans.h"
#include "problemHelper.h"


void solverFISTA(std::vector<double> &b, std::vector<double> &x, double &Lipschisz,
                 double &lambda, double &eta, vector<unsigned int> &indexKeep, int limit);
void updateProximalStep(std::vector<double> &b, double &Lipschisz, double &lambda, 
						std::vector<double> &x, std::vector<double> &fx, std::vector<double> &fGradient, std::vector<double> &sol,
						unsigned int &n, vector<unsigned int> &indexKeep);
void computeApproxQValue(std::vector<double> &b, std::vector<double> &x, std::vector<double> &y,
                         std::vector<double> &fy, std::vector<double> &fGradient, unsigned int &n, double &lambda, double &obj, double &Lipschisz,
                         vector<unsigned int> &indexKeep);


void solverFISTA(std::vector<double> &b, std::vector<double> &x, double &Lipschisz,
                 double &lambda, double &eta, vector<unsigned int> &indexKeep, int limit) {

	double tk = 1.0;
	double tkp = 1.0;
	double newLipschisz = Lipschisz;
	unsigned int n = x.size();
	std::vector<double> y(n);
	std::vector<double> oldX(n);
	std::vector<double> deltaX(n);
	std::vector<double> potentialX(n);
	std::vector<double> fy(n);
	std::vector<double> fx(n);
	std::vector<double> yGradient(n);
	cblas_dcopy(n, &x[0], 1, &y[0], 1);

	std::vector<double> ATb(n);
	forwardDWT(b, ATb, n, indexKeep);
	double objF;
	double objQ;
	double obj_re;
	for (unsigned int iter = 0; iter < limit; iter++) {

		AxMinusb(b, y, fy, n, indexKeep);
		computeGradientFx(b, y, fy, n, yGradient, indexKeep);

		for (unsigned int i = 0; i < 500; i++) {

			double fractional = pow(eta, i);
			newLipschisz = fractional * Lipschisz;
			//compute objective
			updateProximalStep(b, newLipschisz, lambda, y, fy, yGradient, potentialX, n, indexKeep);

			AxMinusb(b, potentialX, fx, n, indexKeep);
			computeObjectiveValue(b, potentialX, fx, n, lambda, objF, indexKeep);			//for (unsigned int j = 0; j < n; j++) cout<<potentialX[j]<<"   ";
			//compute model objective
			computeApproxQValue(b, potentialX, y, fy, yGradient, n, lambda, objQ, newLipschisz, indexKeep);
			//cout <<iter<<"   " << i << "   " << objF << "  " << objQ << "   " << newLipschisz << endl;
	
			if (objF <= objQ) {
				// x_k = P_L(y_k)
				cblas_dcopy(n, &potentialX[0], 1, &x[0], 1);

				cblas_dcopy(n, &x[0], 1, &deltaX[0], 1);
				cblas_daxpy(n, -1.0, &oldX[0], 1, &deltaX[0], 1);

				tkp = 0.5 * ( 1.0 + sqrt(1.0 + 4.0 * tk * tk) );

				// update y
				double alpha = (tk - 1.0) / tkp; //cout<<alpha<<endl;
				cblas_dcopy(n, &x[0], 1, &y[0], 1);
				cblas_daxpy(n, alpha, &deltaX[0], 1, &y[0], 1);		
				cout<<iter<<"   "<<i<<"   "<<objF<<"  "<<objQ<<"   "<<"   "<<(obj_re - objF)/obj_re<<"   "<<newLipschisz<<endl;
				break;
			}

		}
		
		Lipschisz = newLipschisz-0.05;
		tk = tkp;
		cblas_dcopy(n, &x[0], 1, &oldX[0], 1); 		//for (unsigned int j = 0; j < n; j++) cout<<j<<"  "<<potentialX[j]<<"   "<<y[j]<<endl;;
		if (abs((obj_re - objF)/obj_re) < 1e-3 )  
			break;
		obj_re = objF;
	}

}

void updateProximalStep(std::vector<double> &b, double &Lipschisz, double &lambda, 
	std::vector<double> &x, std::vector<double> &fx, std::vector<double> &fGradient, std::vector<double> &sol, 
	unsigned int &n, vector<unsigned int> &indexKeep) {

	std::vector<double> innerOperator(n);
	double beta = lambda / Lipschisz;
	cblas_dcopy(n, &x[0], 1, &innerOperator[0], 1);
	cblas_daxpy(n, -1.0 / Lipschisz, &fGradient[0], 1, &innerOperator[0], 1);

	unsigned int j = 0;	
	for (unsigned int idx = 0; idx < n; idx++) {
		if (innerOperator[idx] >= 0)
			sol[idx] = max(innerOperator[idx] - beta, 0.0);
		else
			sol[idx] = -max(-innerOperator[idx] - beta, 0.0);	
//cout<<idx<<"   "<<sol[idx]<<"  "<<b[idx]<<endl;
	}

}

void computeApproxQValue(std::vector<double> &b, std::vector<double> &x, std::vector<double> &y,
                         std::vector<double> &fy, std::vector<double> &yGradient, unsigned int &n,
                         double &lambda, double &obj, double &Lipschisz, vector<unsigned int> &indexKeep) {

	std::vector<double> xMinusY(n);
	double dotProduct;

	double normL2f = cblas_dnrm2(n, &fy[0], 1);

	cblas_dcopy(n, &x[0], 1, &xMinusY[0], 1);
	cblas_daxpy(n, -1.0, &y[0], 1, &xMinusY[0], 1);
	dotProduct = cblas_ddot(n, &xMinusY[0], 1, &yGradient[0], 1);

	double normL2XMinusY = cblas_dnrm2(n, &xMinusY[0], 1);

	double normL1g = cblas_dasum(n, &x[0], 1);

	obj = normL2f * normL2f + dotProduct + 0.5 * Lipschisz * normL2XMinusY * normL2XMinusY + lambda * normL1g;

}


#endif