#include <math.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <vector>
#include <iomanip>      // std::setprecision
#include <set>
#include <sstream>
using namespace std;

#include <boost/timer/timer.hpp>
#include <gsl/gsl_blas.h> // for BLAS

#include "fista_.h"
#include "gpsrbasic.h"
#include "readImage.h"
#include "DWTtrans.h"
#include "parameters.h"

int main(int argc, char *argv[]) {

	char* inputFile; //= "test5.bmp";
	int method = 1;
	unsigned int n = 10000;
	unsigned int width = 0;
	unsigned int height = 0;
	double tau = 0.7; // percentage of pixels to blur
	double lambda = 1; // in the objective function
	double L = 1.9;	// Lipschiz for fista
	double eta = 1.2; // multiplier for fista

	double mu = 0.5; // sufficient decreasing parameter for GPSR
	double backtrack = 0.7; // multiplier for backtrack line search
	vector<double> dataR(n);
	vector<double> dataG(n);
	vector<double> dataB(n);
	parseDistributedOptions(inputFile, method, tau, lambda, argc, argv);

	readImage(inputFile, dataR, dataG, dataB, width, height); //cout<<width<<"  "<<height;
	unsigned int N = width * height;
	unsigned int numBlackPixel = floor(N * tau);
	unsigned int numPixel = N - numBlackPixel;
	vector<unsigned int> indexBlur(numBlackPixel); //cout<<numBlackPixel<<endl;
	vector<unsigned int> indexKeep(numPixel); //cout<<numPixel<<endl;
	disturbImage(dataR, dataG, dataB, N, numBlackPixel, indexBlur, indexKeep);
	saveBlurImage(inputFile, dataR, dataG, dataB, width, height);

	vector<double> xR(N);
	vector<double> xG(N);
	vector<double> xB(N);
	// vector<double> AxR(N);
	// vector<double> AxG(N);
	// vector<double> AxB(N);
	//inverseDWT3(xR, xG, xB, AxR, AxG, AxB, N, indexBlur);
	boost::timer::cpu_timer timer;
	if (method ==1) {
		solverFISTA(dataR, xR, L, lambda, eta, indexKeep, 1000);
		solverFISTA(dataG, xG, L, lambda, eta, indexKeep, 1000);
		solverFISTA(dataB, xB, L, lambda, eta, indexKeep, 1000);
		PureInverseDWT3(xR, xG, xB, N);
		saveImage(inputFile, tau, method, xR, xG, xB, width, height);
	}
	else if (method == 2) {
		solverGPSRBasic(dataR, xR, lambda, indexKeep, mu, backtrack);
		solverGPSRBasic(dataG, xG, lambda, indexKeep, mu, backtrack);
		solverGPSRBasic(dataB, xB, lambda, indexKeep, mu, backtrack);
		PureInverseDWT3(xR, xG, xB, N);
		saveImage(inputFile, 0, method, xR, xG, xB, width, height);
	}
	else if (method == 3) {
		int iter = 25;
		for (int i = 0; i < 6; i++) {
			PureForwardDWT3(xR, xG, xB, N);
			solverFISTA(dataR, xR, L, lambda, eta, indexKeep, iter);
			solverFISTA(dataG, xG, L, lambda, eta, indexKeep, iter);
			solverFISTA(dataB, xB, L, lambda, eta, indexKeep, iter);
			PureInverseDWT3(xR, xG, xB, N);
			int rec = (i+1) * iter;
			saveImage(inputFile, rec, method, xR, xG, xB, width, height);
		}
	}
	boost::timer::cpu_times elapsed = timer.elapsed();
	std::cout << " WALLCLOCK TIME: " << elapsed.wall / 1e9 << " seconds" << std::endl;

	//PureInverseDWT3(xR, xG, xB, N);
	//saveImage(inputFile, xR, xG, xB, width, height);


	return 0;
}
