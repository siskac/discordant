#include <R.h>
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

void em_normal_partial_concordant(double *x, double *y, double *zxy, int *n, double *pi, double *mu, double *sigma, double *nu, double *tau, int *g, double *loglik, double *tol, int *restriction, int *constrain, int *iteration, int *convergence) {
	int i, j, k, flag, iter;
	double temp, loglik2;


	/*theoretical restriction on the null*/
	if((*restriction)>0){
	/*initialization*/
		flag=1;
		(*loglik) = 0;
		for(k=0;k<*n;k++){	
			(*loglik) = (*loglik) - 0.5*x[k]*x[k] - 0.5*log(2*PI) - 0.5*y[k]*y[k] - 0.5*log(2*PI);
		}
		/*iteration*/
		
		iter = 0;
		while(flag>0 & iter<(*iteration)){
		/*M step for pi, mu and sigma*/
			for(i=0;i<(*g);i++){
				pi[i] = 0;
				for(k=0;k<(*n);k++){
					for(j=0;j<(*g);j++){
						pi[i] = pi[i] + zxy[(j*(*g)+i)*(*n)+k];
					}
				}
				mu[i] = 0;
				for(k=0;k<(*n);k++){
					for(j=0;j<(*g);j++){
						mu[i] = mu[i] + zxy[(j*(*g)+i)*(*n)+k]*x[k];
					}
				}
				mu[i] = mu[i] / pi[i];
				sigma[i] = 0;
				for(k=0;k<(*n);k++){
					for(j=0;j<(*g);j++){
						sigma[i] = sigma[i] + zxy[(j*(*g)+i)*(*n)+k]*(x[k]-mu[i])*(x[k]-mu[i]);
					}
				}
				sigma[i] = sigma[i] / pi[i];
			}
			for(j=0;j<(*g);j++){
				pi[j] = 0;
				for(k=0;k<(*n);k++){
					for(i=0;i<(*g);i++){
						pi[j] = pi[j] + zxy[(j*(*g)+i)*(*n)+k];
					}
				}
				nu[j] = 0;
				for(k=0;k<(*n);k++){
					for(i=0;i<(*g);i++){
						nu[j] = nu[j] + zxy[(j*(*g)+i)*(*n)+k]*y[k];
					}
				}
				nu[j] = nu[j] / pi[j];
				tau[j] = 0;
				for(k=0;k<(*n);k++){
					for(i=0;i<(*g);i++){
						tau[j] = tau[j] + zxy[(j*(*g)+i)*(*n)+k]*(y[k]-nu[j])*(y[k]-nu[j]);
					}
				}
				tau[j] = tau[j] / pi[j];
			}
			for(i=0;i<(*g);i++){
				for(j=0;j<(*g);j++){
					pi[i*(*g)+j] = 0;
					for(k=0;k<(*n);k++){
						pi[i*(*g)+j] = pi[i*(*g)+j] + zxy[(j*(*g)+i)*(*n)+k];
					}
					pi[i*(*g)+j] = pi[i*(*g)+j] / (double)(*n);
				}
			}
			/*constrain*/
			for(i=0;i<(*g);i++){
				if(constrain[i]==0-1){
					if(mu[i]>0.0){
						mu[i]=0.0;
					}
					if(nu[i]>0.0){
						nu[i]=0.0;
					}
				}
				if(constrain[i]==0){
					mu[i] = 0.0;
					sigma[i] = 1.0;
					nu[i] = 0.0;
					tau[i] = 1.0;
				}
				if(constrain[i]==0+1){
					if(mu[i]<0.0){
						mu[i]=0.0;
					}
					if(nu[i]<0.0){
						nu[i]=0.0;
					}
				}
			}
			/*check convergence*/
			loglik2 = (*loglik);
			(*loglik) = 0;
			for(k=0;k<(*n);k++){
				temp = 0;
				for(i=0;i<(*g);i++){
					for(j=0;j<(*g);j++){
						temp = temp + pi[i*(*g)+j] * ( exp(0-0.5*(x[k]-mu[i])*(x[k]-mu[i])/sigma[i])/sqrt(2*PI*sigma[i]) )*( exp(0-0.5*(y[k]-nu[j])*(y[k]-nu[j])/tau[j])/sqrt(2*PI*tau[j]) );
					}
				}
				(*loglik) = (*loglik) + log(temp);
			}
			if(fabs((*loglik)-loglik2)<(*tol)){
				flag = 0;
			}
			/*E step for z*/
			for(k=0;k<(*n);k++){
				temp = 0;
				for(i=0;i<(*g);i++){
					for(j=0;j<(*g);j++){
						temp = temp + pi[i*(*g)+j] * ( exp(0-0.5*(x[k]-mu[i])*(x[k]-mu[i])/sigma[i])/sqrt(2*PI*sigma[i]) )*( exp(0-0.5*(y[k]-nu[j])*(y[k]-nu[j])/tau[j])/sqrt(2*PI*tau[j]) );
					}
				}
				for(i=0;i<(*g);i++){
					for(j=0;j<(*g);j++){
						zxy[(j*(*g)+i)*(*n)+k] = pi[i*(*g)+j] * ( exp(0-0.5*(x[k]-mu[i])*(x[k]-mu[i])/sigma[i])/sqrt(2*PI*sigma[i]) )*( exp(0-0.5*(y[k]-nu[j])*(y[k]-nu[j])/tau[j])/sqrt(2*PI*tau[j]) )/temp;
					}
				}
			}
			iter = iter + 1;
		}
		(*convergence) = 0;
		if(iter<(*iteration)){
			(*convergence) = 1;
		}
	}
	/*no restriction on the null*/
	else{
	/*initialization*/
		flag=1;
		(*loglik) = 0;
		for(k=0;k<(*n);k++){
			(*loglik) = (*loglik) - 0.5*x[k]*x[k] - 0.5*log(2*PI) - 0.5*y[k]*y[k] - 0.5*log(2*PI);
		}
		/*iteration*/
		iter = 0;
		while(flag>0 & iter<(*iteration)){
		/*M step for pi, mu and sigma*/
			for(i=0;i<(*g);i++){
				pi[i] = 0;
				for(k=0;k<(*n);k++){
					for(j=0;j<(*g);j++){
						pi[i] = pi[i] + zxy[(j*(*g)+i)*(*n)+k];
					}
				}
				mu[i] = 0;
				for(k=0;k<(*n);k++){
					for(j=0;j<(*g);j++){
						mu[i] = mu[i] + zxy[(j*(*g)+i)*(*n)+k]*x[k];
					}
				}
				mu[i] = mu[i] / pi[i];
				sigma[i] = 0;
				for(k=0;k<(*n);k++){
					for(j=0;j<(*g);j++){
						sigma[i] = sigma[i] + zxy[(j*(*g)+i)*(*n)+k]*(x[k]-mu[i])*(x[k]-mu[i]);
					}
				}
				sigma[i] = sigma[i] / pi[i];
			}
			for(j=0;j<(*g);j++){
				pi[j] = 0;
				for(k=0;k<(*n);k++){
					for(i=0;i<(*g);i++){
						pi[j] = pi[j] + zxy[(j*(*g)+i)*(*n)+k];
					}
				}
				nu[j] = 0;
				for(k=0;k<(*n);k++){
					for(i=0;i<(*g);i++){
						nu[j] = nu[j] + zxy[(j*(*g)+i)*(*n)+k]*y[k];
					}
				}
				nu[j] = nu[j] / pi[j];
				tau[j] = 0;
				for(k=0;k<(*n);k++){
					for(i=0;i<(*g);i++){
						tau[j] = tau[j] + zxy[(j*(*g)+i)*(*n)+k]*(y[k]-nu[j])*(y[k]-nu[j]);
					}
				}
				tau[j] = tau[j] / pi[j];
			}
			for(i=0;i<(*g);i++){
				for(j=0;j<(*g);j++){
					pi[i*(*g)+j] = 0;
					for(k=0;k<(*n);k++){
						pi[i*(*g)+j] = pi[i*(*g)+j] + zxy[(j*(*g)+i)*(*n)+k];
					}
					pi[i*(*g)+j] = pi[i*(*g)+j] / (double)(*n);
				}
			}
			/*check convergence*/
			loglik2 = (*loglik);
			(*loglik) = 0;
			for(k=0;k<(*n);k++){
				temp = 0;
				for(i=0;i<(*g);i++){
					for(j=0;j<(*g);j++){
						temp = temp + pi[i*(*g)+j] * ( exp(0-0.5*(x[k]-mu[i])*(x[k]-mu[i])/sigma[i])/sqrt(2*PI*sigma[i]) )*( exp(0-0.5*(y[k]-nu[j])*(y[k]-nu[j])/tau[j])/sqrt(2*PI*tau[j]) );
					}
				}
				(*loglik) = (*loglik) + log(temp);
			}
			if(fabs((*loglik)-loglik2)<(*tol)){
				flag = 0;
			}
			/*E step for z*/
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
			iter = iter + 1;
		}
		(*convergence) = 0;
		if(iter<(*iteration)){
			(*convergence) = 1;
		}
	}
}

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

R_CallMethodDef callMethods[]  = {
  {"em_normal_partial_concordant", (DL_FUNC) &em_normal_partial_concordant, 16},
  {"subsampling", (DL_FUNC) &subsampling, 11},
  {NULL, NULL, 0}
};

void R_init_discordant(DllInfo *info) {
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

