#ifndef GPSRBASIC_H_
#define GPSRBASIC_H_


#include <cmath>
using namespace std;
#include "problemHelper.h"

void solverGPSRBasic(std::vector<double> &b, std::vector<double> &x, double &lambda,
			vector<unsigned int> &indexKeep, double &mu, double &backtrack) {

	double tol = 1e-3;
	int iter = 0;
	int iterOuter = 0;
	double objF;
	double objFNew;
	double objFOld;
	double stopObjTest;
	double stepSize;
	double stepSizeInitial;

	unsigned int n = x.size();
	std::vector<double> funGradient(n);
	std::vector<double> fx(n);
	std::vector<double> ones(n, 1.0);
	std::vector<double> uGrad(n);
	std::vector<double> vGrad(n);
	std::vector<double> u(n);
	std::vector<double> v(n);
	std::vector<double> uOld(n);
	std::vector<double> vOld(n);
	std::vector<double> uCondGrad(n);
	std::vector<double> vCondGrad(n);
	std::vector<double> Auv(n);
	std::vector<double> du(n);
	std::vector<double> dv(n);
	std::vector<double> uNew(n);
	std::vector<double> vNew(n);
	std::vector<double> xNew(n);
	std::vector<double> dx(n);
	std::vector<double> uvMin(n);
	//std::vector<double> ATb(n);
	//forwardDWT(b, ATb, n, indexKeep);
	// compute initial objective value
	AxMinusb(b, x, fx, n, indexKeep);
	computeObjectiveValue(b, x, fx, n, lambda, objF, indexKeep);

	// start iteration
	while (stopObjTest >= tol) {
		// compute gradient
		//computeGradientFx_fast(ATb, x, fx, n, funGradient, indexKeep);
		computeGradientFx(b, x, fx, n, funGradient, indexKeep);
		cblas_dcopy(n, &funGradient[0], 1, &uGrad[0], 1);
		cblas_daxpy(n, lambda, &ones[0], 1, &uGrad[0], 1);
		cblas_dcopy(n, &funGradient[0], 1, &vGrad[0], 1);
		cblas_dscal(n, -1.0, &vGrad[0], 1);
		cblas_daxpy(n, lambda, &ones[0], 1, &vGrad[0], 1);
		//compute initial stepSize, reference Matlab code line 432
		cblas_dcopy(n, &u[0], 1, &uOld[0], 1);
		cblas_dcopy(n, &v[0], 1, &vOld[0], 1);
		cblas_dcopy(n, &uGrad[0], 1, &uCondGrad[0], 1);
		cblas_dcopy(n, &vGrad[0], 1, &vCondGrad[0], 1);
		for (unsigned int idx = 0; idx < n; idx++) {
			if ( uOld[idx] <= 0 && uGrad[idx] >= 0 )
				uCondGrad[idx] = 0;			
			if ( vOld[idx] <= 0 && vGrad[idx] >= 0 )
				vCondGrad[idx] = 0;
		}

		cblas_daxpy(n, -1.0, &vCondGrad[0], 1, &uCondGrad[0], 1);
		inverseDWT(uCondGrad, Auv, n, indexKeep);
		double denom = cblas_dnrm2(n, &Auv[0], 1);
		denom = denom * denom;
		double nom1 = cblas_ddot(n, &uCondGrad[0], 1, &uGrad[0], 1);
		double nom2 = cblas_ddot(n, &vCondGrad[0], 1, &vGrad[0], 1);

		stepSizeInitial = (nom1 + nom2) / denom;
		stepSize = stepSizeInitial;

		iter = 0;
		while (1){
			// calculate step for this stepSize
			for (unsigned int idx = 0; idx < n; idx++) {
				du[idx] = max(u[idx] - stepSize * uGrad[idx], 0.0) - u[idx];
				dv[idx] = max(v[idx] - stepSize * vGrad[idx], 0.0) - v[idx];
			}
			cblas_dcopy(n, &u[0], 1, &uNew[0], 1);
			cblas_daxpy(n, 1.0, &du[0], 1, &uNew[0], 1);			
			cblas_dcopy(n, &v[0], 1, &vNew[0], 1);
			cblas_daxpy(n, 1.0, &dv[0], 1, &vNew[0], 1);
			cblas_dcopy(n, &du[0], 1, &dx[0], 1);
			cblas_daxpy(n, -1.0, &dv[0], 1, &dx[0], 1);
			cblas_dcopy(n, &x[0], 1, &xNew[0], 1);
			cblas_daxpy(n, 1.0, &dx[0], 1, &xNew[0], 1);
			// compute new objective
			AxMinusb(b, xNew, fx, n, indexKeep);
			computeObjectiveValue(b, xNew, fx, n, lambda, objFNew, indexKeep);
			// test sufficient decreasing condition
			double temp1 = cblas_ddot(n, &du[0], 1, &uGrad[0], 1);
			double temp2 = cblas_ddot(n, &dv[0], 1, &vGrad[0], 1);
			if (objFNew <= objF + mu * (temp1 + temp2) ){
				//cout << iter<<"  "<< objFNew <<"  " << objF + mu * (temp1 + temp2) <<endl;
				//cout<< stepSizeInitial<< "   "<< stepSize<<endl;
				break;
			}
			stepSize *= backtrack;
			iter++;
		}
		cblas_dcopy(n, &uNew[0], 1, &u[0], 1);
		cblas_dcopy(n, &vNew[0], 1, &v[0], 1);
		objFOld = objF;
		objF = objFNew;
		for (unsigned int idx = 0; idx < n; idx++) {
			uvMin[idx] = min(u[idx], v[idx]); 
		}
		cblas_daxpy(n, -1.0, &uvMin[0], 1, &u[0], 1);
		cblas_daxpy(n, -1.0, &uvMin[0], 1, &v[0], 1);
		cblas_dcopy(n, &u[0], 1, &x[0], 1);
		cblas_daxpy(n, -1.0, &v[0], 1, &x[0], 1);

		stopObjTest = abs((objF - objFOld) / objFOld);
		cout << iterOuter<<"  "<<iter<<"  "<< objF <<"  " << stopObjTest <<"  "<< stepSizeInitial<< "  "<< stepSize<<endl;
		iterOuter++;

	}


}



#endif