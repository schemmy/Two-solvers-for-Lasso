#ifndef DWTTRANS_H
#define DWTTRANS_H

#include <gsl/gsl_wavelet.h>
#include <gsl/gsl_wavelet2d.h>


void inverseDWT3(std::vector<double> &xR, std::vector<double> &xG, std::vector<double> &xB,
	std::vector<double> &AxR, std::vector<double> &AxG, std::vector<double> &AxB,
	unsigned int &N, vector<unsigned int> &indexBlur){
	gsl_wavelet *w;
	gsl_wavelet_workspace *work;

	w = gsl_wavelet_alloc (gsl_wavelet_daubechies, 4);
	work = gsl_wavelet_workspace_alloc (N);

	gsl_wavelet_transform_inverse (w, &xR[0], 1, N, work);
	gsl_wavelet_transform_inverse (w, &xG[0], 1, N, work);
	gsl_wavelet_transform_inverse (w, &xB[0], 1, N, work);
	//gsl_wavelet2d_transform_inverse (w, x, 1, width, height, work); width == height

	unsigned int j = 0;	
	// can write the following in sparse matrix vector mulplication
	for (unsigned int i = 0; i < N; i++){
		if (i != indexBlur[j]){
			AxR[i] = xR[i]; 
			AxG[i] = xG[i]; 
			AxB[i] = xB[i]; 
			//cout<<i<<"   ";
		}
		else{
			j++;
			AxR[i] = 0; 
			AxG[i] = 0; 
			AxB[i] = 0; 
		}
	}


	gsl_wavelet_free (w);
	gsl_wavelet_workspace_free (work);

}

void inverseDWT(std::vector<double> &x, std::vector<double> &Ax, 
	unsigned int &N, vector<unsigned int> &indexKeep){
	gsl_wavelet *w;
	gsl_wavelet_workspace *work;

	w = gsl_wavelet_alloc (gsl_wavelet_daubechies, 4);
	work = gsl_wavelet_workspace_alloc (N);
	unsigned int n = sqrt(N);
	std::vector<double> temp(N);
	
	cblas_dcopy(N, &x[0], 1, &temp[0], 1);

	//gsl_wavelet_transform_inverse (w, &Ax[0], 1, N, work);
	gsl_wavelet2d_transform_inverse (w, &temp[0], n, n, n, work); //width == height

	//unsigned int j = 0;	
	// can write the following in sparse matrix vector mulplication
	// for (unsigned int i = 0; i < N; i++){
	// 	if (i == indexBlur[j]){	
	// 		j++;
	// 		Ax[i] = 0; 
	// 	}
	// }
	cblas_dscal(N, 0, &Ax[0], 1);
	for (unsigned int i = 0; i < indexKeep.size(); i++)
		Ax[indexKeep[i]] = temp[indexKeep[i]];

	gsl_wavelet_free (w);
	gsl_wavelet_workspace_free (work);

}

void forwardDWT(std::vector<double> &x, std::vector<double> &Ax, 
	unsigned int &N, vector<unsigned int> &indexKeep){
	gsl_wavelet *w;
	gsl_wavelet_workspace *work;

	w = gsl_wavelet_alloc (gsl_wavelet_daubechies, 4);
	work = gsl_wavelet_workspace_alloc (N);
	unsigned int n = sqrt(N);
	std::vector<double> temp(N);

	//unsigned int j = 0;	
	// can write the following in sparse matrix vector mulplication
	// for (unsigned int i = 0; i < N; i++){
	// 	if (i != indexBlur[j]){
	// 		Ax[i] = x[i]; 
	// 	}
	// 	else{
	// 		j++;
	// 		Ax[i] = 0; 
	// 	}
	// }
	cblas_dcopy(N, &x[0], 1, &temp[0], 1);
	
	cblas_dscal(N, 0, &Ax[0], 1);
	for (unsigned int i = 0; i < indexKeep.size(); i++)
		Ax[indexKeep[i]] = temp[indexKeep[i]];
	//gsl_wavelet_transform_forward(w, &Ax[0], 1, N, work);
	gsl_wavelet2d_transform_forward (w, &Ax[0], n, n, n, work); //width == height


	gsl_wavelet_free (w);
	gsl_wavelet_workspace_free (work);

}


void PureInverseDWT(std::vector<double> &x, unsigned int &N){
	gsl_wavelet *w;
	gsl_wavelet_workspace *work;

	w = gsl_wavelet_alloc (gsl_wavelet_daubechies, 4);
	work = gsl_wavelet_workspace_alloc (N);
	unsigned int n = sqrt(N);

	//gsl_wavelet_transform_inverse (w, &x[0], 1, N, work);
	gsl_wavelet2d_transform_inverse (w, &x[0], n, n, n, work); //width == height

	gsl_wavelet_free (w);
	gsl_wavelet_workspace_free (work);

}

void PureInverseDWT3(std::vector<double> &xR, std::vector<double> &xG, std::vector<double> &xB,
		 unsigned int &N){
	gsl_wavelet *w;
	gsl_wavelet_workspace *work;

	w = gsl_wavelet_alloc (gsl_wavelet_daubechies, 4);
	work = gsl_wavelet_workspace_alloc (N);
	unsigned int n = sqrt(N);

	//gsl_wavelet_transform_inverse (w, &x[0], 1, N, work);
	gsl_wavelet2d_transform_inverse (w, &xR[0], n, n, n, work); //width == height
	gsl_wavelet2d_transform_inverse (w, &xG[0], n, n, n, work); //width == height
	gsl_wavelet2d_transform_inverse (w, &xB[0], n, n, n, work); //width == height

	gsl_wavelet_free (w);
	gsl_wavelet_workspace_free (work);

}

void PureForwardDWT3(std::vector<double> &xR, std::vector<double> &xG, std::vector<double> &xB,
		 unsigned int &N){
	gsl_wavelet *w;
	gsl_wavelet_workspace *work;

	w = gsl_wavelet_alloc (gsl_wavelet_daubechies, 4);
	work = gsl_wavelet_workspace_alloc (N);
	unsigned int n = sqrt(N);

	//gsl_wavelet_transform_inverse (w, &x[0], 1, N, work);
	gsl_wavelet2d_transform_forward (w, &xR[0], n, n, n, work); //width == height
	gsl_wavelet2d_transform_forward (w, &xG[0], n, n, n, work); //width == height
	gsl_wavelet2d_transform_forward (w, &xB[0], n, n, n, work); //width == height

	gsl_wavelet_free (w);
	gsl_wavelet_workspace_free (work);

}
#endif

// void inverseDWT3(std::vector<double> &xR, std::vector<double> &xG, std::vector<double> &xB,
// 	std::vector<double> &AxR, std::vector<double> &AxG, std::vector<double> &AxB,
// 	unsigned int &N, unsigned int &numPixel, vector<unsigned int> &indexBlur){
// 	gsl_wavelet *w;
// 	gsl_wavelet_workspace *work;

// 	w = gsl_wavelet_alloc (gsl_wavelet_daubechies, 4);
// 	work = gsl_wavelet_workspace_alloc (N);

// 	gsl_wavelet_transform_inverse (w, &xR[0], 1, N, work);
// 	gsl_wavelet_transform_inverse (w, &xG[0], 1, N, work);
// 	gsl_wavelet_transform_inverse (w, &xB[0], 1, N, work);
// 	//gsl_wavelet2d_transform_inverse (w, x, 1, width, height, work); width == height

// 	unsigned int j = 0;	
// 	unsigned int idx = 0;	
// 	// can write the following in sparse matrix vector mulplication
// 	for (unsigned int i = 0; i < N; i++){
// 		if (i != indexBlur[j]){
// 			AxR[idx] = xR[i]; 
// 			AxG[idx] = xG[i]; 
// 			AxB[idx] = xB[i]; 
// 			idx++; //cout<<i<<"   ";
// 		}
// 		else{
// 			j++;
// 		}
// 	}


// 	gsl_wavelet_free (w);
// 	gsl_wavelet_workspace_free (work);

// }