#include <R.h>

void subsampling(double *x, double *y, double *zxy, int *n, double *pi, double *mu, double *sigma, double *nu, double *tau, int *g) {
	int i, j, k;
	double temp;
	for(k=0;k<(*n);k++){
		temp = 0;
		for(i=0;i<(*g);i++){
			for(j=0;j<(*g);j++){
				temp = temp + pi[i*(*g)+j] * ( exp(0-0.5*(x[k]-mu[i])*(x[k]-mu[i])/sigma[i])/sqrt(2*PI*sigma[i]) )*( exp(0-0.5*(y[k]-nu[j])*(y[k]-nu[j])/tau[j])/sqrt(2*PI*tau[j]) );
			}
		}
		for(i=0;i<(*g);i++) {
			for(j=0;j<(*g);j++){
				zxy[(j*(*g)+i)*(*n)+k] = zxy[(j*(*g)+i)*(*n)+k] + pi[i*(*g)+j] * ( exp(0-0.5*(x[k]-mu[i])*(x[k]-mu[i])/sigma[i])/sqrt(2*PI*sigma[i]) )*( exp(0-0.5*(y[k]-nu[j])*(y[k]-nu[j])/tau[j])/sqrt(2*PI*tau[j]) )/temp;
			}
		}
	}

}
