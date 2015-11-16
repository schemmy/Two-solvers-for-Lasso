#ifndef READIMAGE_H
#define READIMAGE_H

#include "CImg.h"
using namespace cimg_library;
#include <gsl/gsl_sort.h>
#include <string>

void readImage(const char* fileName, std::vector<double> &dataR, std::vector<double> &dataG,
	std::vector<double> &dataB, unsigned int &width, unsigned int &height){

	string line1 = fileName;    
    string line2 = ".bmp";
	string totalLine = line1 + line2;
	const char *cstr = totalLine.c_str();
	CImg<double> image(cstr);
	width = image.width();
	height = image.height();
	unsigned int N = width * height;
	dataR.resize(N);
	dataG.resize(N);
	dataB.resize(N);
	for (int i = 0; i < width; i++){
		for (int j = 0; j < height; j++){
			unsigned int idx = i * height + j;
			dataR[idx] = image(i, j, 0, 0);// / 256;
			dataG[idx] = image(i, j, 0, 1);// / 256;
			dataB[idx] = image(i, j, 0, 2);// / 256;
		}
	}

}


void saveImage(const char* fileName, double iter, int method, std::vector<double> &dataR, std::vector<double> &dataG, 
	std::vector<double> &dataB, unsigned int &width, unsigned int &height){

	CImg<double> imageToSave(width, height, 1, 3, 0);
	for (int i = 0; i < width; i++){
		for (int j = 0; j < height; j++){
			unsigned int idx = i * height + j;
			imageToSave(i, j, 0, 0) = dataR[idx];// * 256;
			imageToSave(i, j, 0, 1) = dataG[idx];// * 256;
			imageToSave(i, j, 0, 2) = dataB[idx];// * 256;
		}
	}
	imageToSave.normalize(0, 255);
	CImgDisplay dark_disp (imageToSave, "New Image", 0);
	string line1 = fileName; 
	string line11 = "_"; 
	string line12 = to_string(method);
    string line2 = "_result_";
    string line3 = to_string(iter);
    string line4 = ".bmp";
	string totalLine = line1 + line11 + line12 + line2 + line3 + line4;
	const char *cstr = totalLine.c_str();
    imageToSave.save(cstr);
}

void saveBlurImage(const char* fileName, std::vector<double> &dataR, std::vector<double> &dataG, 
	std::vector<double> &dataB, unsigned int &width, unsigned int &height){

	CImg<double> imageToSave(width, height, 1, 3, 0);
	for (int i = 0; i < width; i++){
		for (int j = 0; j < height; j++){
			unsigned int idx = i * height + j;
			imageToSave(i, j, 0, 0) = dataR[idx];// * 256;
			imageToSave(i, j, 0, 1) = dataG[idx];// * 256;
			imageToSave(i, j, 0, 2) = dataB[idx];// * 256;
		}
	}
	imageToSave.normalize(0, 255);
	CImgDisplay dark_disp (imageToSave, "New Image", 0);
	string line1 = fileName;    
    string line2 = "_blur";
    string line3 = ".bmp";
	string totalLine = line1 + line2 + line3;
	const char *cstr = totalLine.c_str();
    imageToSave.save(cstr);

}

void disturbImage(std::vector<double> &dataR, std::vector<double> &dataG,
	std::vector<double> &dataB, unsigned int &N, unsigned int &numBlackPixel,
	std::vector<unsigned int> &indexBlur, std::vector<unsigned int> &indexKeep){

	double *randd = new double[N];
	double *index = new double[N];
	double *luckyIndex = new double[numBlackPixel];
	double *unluckyIndex = new double[N - numBlackPixel];

	for (unsigned int i = 0; i < N; i++){
		randd[i] = rand() / (RAND_MAX + 0.0);
		index[i] = i;
	}
	
	gsl_sort2(randd, 1, index, 1, size_t(N));

	for (unsigned int i = 0; i < numBlackPixel; i++){
		luckyIndex[i] = index[i];
	}
	for (unsigned int i = numBlackPixel; i < N; i++){
		unluckyIndex[i - numBlackPixel] = index[i];
	}	
	gsl_sort(luckyIndex, 1, size_t(numBlackPixel));
	for (unsigned int i = 0; i < numBlackPixel; i++){
		indexBlur[i] = (unsigned int)luckyIndex[i];
		dataR[indexBlur[i]] = 0.0;
		dataG[indexBlur[i]] = 0.0;
		dataB[indexBlur[i]] = 0.0;
	}
	for (unsigned int i = 0; i < N - numBlackPixel; i++)
		indexKeep[i] = (unsigned int)unluckyIndex[i];
}



#endif

// void rescalePixel(std::vector<double> &xR, std::vector<double> &xG, 
// 	std::vector<double> &xB, unsigned int &N){

// 	vector<double> small(1);
// 	vector<double> big(1);

//  	gsl_sort_smallest(&small[0], 1, &xR[0], 1, N);
//  	gsl_sort_largest(&big[0], 1, &xR[0], 1, N);

// 	for (unsigned int i = 1; i < N; i++)
// 		xR[i] = (xR[i] - small[0]) / (big[0] - small[0]) * 256;
 	
//  	gsl_sort_smallest(&small[0], 1, &xG[0], 1, N);
//  	gsl_sort_largest(&big[0], 1, &xG[0], 1, N);

// 	for (unsigned int i = 1; i < N; i++)
// 		xG[i] = (xG[i] - small[0]) / (big[0] - small[0]) * 256;
 	
// 	gsl_sort_smallest(&small[0], 1, &xB[0], 1, N);
//  	gsl_sort_largest(&big[0], 1, &xB[0], 1, N);

// 	for (unsigned int i = 1; i < N; i++)
// 		xB[i] = (xB[i] - small[0]) / (big[0] - small[0]) * 256;

// }