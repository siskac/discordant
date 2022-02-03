#include <Rcpp.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <iostream>
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::List em_normal_partial_concordant_cpp(Rcpp::NumericVector x,
                                        Rcpp::NumericVector y,
                                        Rcpp::NumericVector zxy,
                                        int n,
                                        Rcpp::NumericVector pi,
                                        Rcpp::NumericVector mu,
                                        Rcpp::NumericVector sigma,
                                        Rcpp::NumericVector nu,
                                        Rcpp::NumericVector tau,
                                        int g,
                                        double loglik,
                                        double tol,
                                        int restriction,
                                        Rcpp::NumericVector constrain,
                                        int iteration,
                                        int convergence) {
    int i, j, k, flag, iter;
    double temp, loglik2;
    
    /*initialization*/
    flag = 1;
    loglik = 0;
    for (k = 0; k < n; k++) {
        loglik = loglik - 0.5*x[k]*x[k] - 0.5*log(2*M_PI) - 
            0.5*y[k]*y[k] - 0.5*log(2*M_PI);
    }
    
    /*iteration*/
    iter = 0;
    while ((flag > 0) & (iter < iteration)) {

        /*M step for pi, mu and sigma*/
        for (i = 0; i < g; i++) {
            
            pi[i] = 0;
            for (k = 0; k < n; k++) {
                for (j = 0; j < g; j++){
                    pi[i] = pi[i] + zxy[(j*g+i)*n+k];
                }
            }
            
            mu[i] = 0;
            for (k = 0; k < n; k++) {
                for (j = 0; j < g; j++) {
                    mu[i] = mu[i] + zxy[(j*g+i)*n+k] * x[k];
                }
            }
            
            if(pi[i] == 0) { stop("Divide-by-zero error. Increase number of features " 
                                  "or reduce number of components. \n  If subsampling=TRUE, "
                                  "you may need to set subsampling=FALSE.\n"
                                  "  %s (Line %d)", __FILE__, __LINE__); }
            
            mu[i] = mu[i] / pi[i];
            
            sigma[i] = 0;
            
            for (k = 0; k < n; k++) {
                for (j = 0; j < g; j++) {
                    sigma[i] = sigma[i] + zxy[(j*g+i)*n+k] * 
                        (x[k]-mu[i]) * (x[k]-mu[i]);
                }
            }

            if(pi[i] == 0) { stop("Divide-by-zero error. Increase number of features " 
                                  "or reduce number of components. \n  If subsampling=TRUE, "
                                  "you may need to set subsampling=FALSE.\n"
                                  "  %s (Line %d)", __FILE__, __LINE__); }
            sigma[i] = sigma[i] / pi[i];

        }

        for (j = 0; j < g; j++) {

            pi[j] = 0;
            for (k = 0; k < n; k++) {
                for (i = 0; i < g; i++) {
                    pi[j] = pi[j] + zxy[(j*g+i)*n+k];
                }
            }
            
            nu[j] = 0;
            for (k = 0; k < n; k++) {
                for (i = 0; i < g; i++) {
                    nu[j] = nu[j] + zxy[(j*g+i)*n+k] * y[k];
                }
            }

            if(pi[j] == 0) { stop("Divide-by-zero error. Increase number of features " 
                                  "or reduce number of components. \n  If subsampling=TRUE, "
                                  "you may need to set subsampling=FALSE.\n"
                                  "  %s (Line %d)", __FILE__, __LINE__); }
            nu[j] = nu[j] / pi[j];

            tau[j] = 0;
            for (k = 0; k < n; k++) {
                for (i = 0; i < g; i++){
                    tau[j] = tau[j] + zxy[(j*g+i)*n+k] * 
                        (y[k]-nu[j]) * (y[k]-nu[j]);
                }
            }

            if(pi[j] == 0) { stop("Divide-by-zero error. Increase number of features " 
                                  "or reduce number of components. \n  If subsampling=TRUE, "
                                  "you may need to set subsampling=FALSE.\n"
                                  "  %s (Line %d)", __FILE__, __LINE__); }
            tau[j] = tau[j] / pi[j];
        }

        for (i = 0; i < g; i++) {
            for (j = 0; j < g; j++) {
                pi[i*g+j] = 0;
                for (k = 0; k < n; k++) {
                    pi[i*g+j] = pi[i*g+j] + zxy[(j*g+i)*n+k];
                }
                pi[i*g+j] = pi[i*g+j] / static_cast<double>(n);
            }
        }

        /* check convergence */
        loglik2 = loglik;
        loglik = 0;
        for (k = 0; k < n; k++) {
            temp = 0;
            for (i = 0; i < g; i++) {
                for (j = 0; j < g; j++) {
                    if(sigma[i] == 0 || tau[j] == 0) { 
                        stop("Divide-by-zero error. Increase number of features " 
                             "or reduce number of components. \n  If subsampling=TRUE, "
                             "you may need to set subsampling=FALSE.\n"
                             "  %s (Line %d)", __FILE__, __LINE__); }
                    if(sigma[i] < 0 || tau[j] < 0) { 
                        stop("Non-real value error. Increase number of features " 
                             "or reduce number of components. \n  If subsampling=TRUE, "
                             "you may need to set subsampling=FALSE.\n"
                             "  %s (Line %d)", __FILE__, __LINE__); }
                    temp = temp + pi[i*g+j] * 
                        (exp(0-0.5*(x[k]-mu[i])*(x[k]-mu[i])/sigma[i])/sqrt(2*M_PI*sigma[i])) * 
                        (exp(0-0.5*(y[k]-nu[j])*(y[k]-nu[j])/tau[j])/sqrt(2*M_PI*tau[j]));
                }
            }
            loglik = loglik + log(temp);
        }
        
        if (fabs(loglik-loglik2) < tol){
            flag = 0;
        }
        
        /*E step for z*/
        for (k = 0; k < n; k++){
            temp = 0;
            for (i = 0; i < g; i++) {
                for (j = 0; j < g; j++) {
                    if(sigma[i] == 0 || tau[j] == 0) { 
                        stop("Divide-by-zero error. Increase number of features " 
                             "or reduce number of components. \n  If subsampling=TRUE, "
                             "you may need to set subsampling=FALSE.\n"
                             "  %s (Line %d)", __FILE__, __LINE__); }
                    if(sigma[i] < 0 || tau[j] < 0) { 
                        stop("Non-real value error.. Increase number of features " 
                             "or reduce number of components. \n  If subsampling=TRUE, "
                             "you may need to set subsampling=FALSE.\n"
                             "  %s (Line %d)", __FILE__, __LINE__); }
                    temp = temp + pi[i*g+j] * 
                        (exp(0-0.5*(x[k]-mu[i])*(x[k]-mu[i])/sigma[i])/sqrt(2*M_PI*sigma[i])) *
                        (exp(0-0.5*(y[k]-nu[j])*(y[k]-nu[j])/tau[j])/sqrt(2*M_PI*tau[j]));
                }
            }
            for (i = 0; i < g;i++) {
                for (j = 0; j < g; j++) {
                    if(sigma[i] == 0 || tau[j] == 0) { 
                        stop("Divide-by-zero error. Increase number of features " 
                             "or reduce number of components. \n  If subsampling=TRUE, "
                             "you may need to set subsampling=FALSE.\n"
                             "  %s (Line %d)", __FILE__, __LINE__); }
                    if(sigma[i] < 0 || tau[j] < 0) { 
                        stop("Non-real value error.. Increase number of features " 
                             "or reduce number of components. \n  If subsampling=TRUE, "
                             "you may need to set subsampling=FALSE.\n"
                             "  %s (Line %d)", __FILE__, __LINE__); }
                    zxy[(j*g+i)*n+k] = zxy[(j*g+i)*n+k] + pi[i*g+j] * 
                        (exp(0-0.5*(x[k]-mu[i])*(x[k]-mu[i])/sigma[i])/sqrt(2*M_PI*sigma[i])) *
                        (exp(0-0.5*(y[k]-nu[j])*(y[k]-nu[j])/tau[j])/sqrt(2*M_PI*tau[j])) / temp;
                }
            }
        }
        iter = iter + 1;
    
    }

    convergence = 0;
    if(iter < iteration) {
        convergence = 1;
    }
    
    Rcpp::List rtn = Rcpp::List::create(x, y, zxy, n, pi, mu, sigma, nu,
                                        tau, g, loglik, tol, restriction,
                                        constrain, iteration, convergence);
    
    return rtn;
}

// [[Rcpp::export]]
Rcpp::List subsampling_cpp(Rcpp::NumericVector x, 
                       Rcpp::NumericVector y,
                       Rcpp::NumericVector zxy,
                       int n, 
                       Rcpp::NumericVector pi, 
                       Rcpp::NumericVector mu, 
                       Rcpp::NumericVector sigma, 
                       Rcpp::NumericVector nu, 
                       Rcpp::NumericVector tau, 
                       int g) {
    int i, j, k;
    double temp;
    
    for (k = 0; k < n; k++) {
        temp = 0;
        for (i = 0; i < g; i++) {
            for (j = 0; j < g; j++) {
                if(sigma[i] == 0 || tau[j] == 0) { stop("Divide-by-zero error (Line 186)."); }
                if(sigma[i] < 0 || tau[j] < 0) { stop("Non-real error (Line 187)."); }
                temp = temp + pi[i*g+j] * 
                    (exp(0-0.5*(x[k]-mu[i])*(x[k]-mu[i])/sigma[i])/sqrt(2*M_PI*sigma[i])) *
                    (exp(0-0.5*(y[k]-nu[j])*(y[k]-nu[j])/tau[j])/sqrt(2*M_PI*tau[j]));
            }
        }
        for (i = 0; i < g; i++) {
            for (j = 0; j < g; j++) {
                if(sigma[i] == 0 || tau[j] == 0) { stop("Divide-by-zero error (Line 195)."); }
                if(sigma[i] < 0 || tau[j] < 0) { stop("Non-real error (Line 196)."); }
                zxy[(j*g+i)*n+k] = zxy[(j*g+i)*n+k] + pi[i*g+j] * 
                    (exp(0-0.5*(x[k]-mu[i])*(x[k]-mu[i])/sigma[i])/sqrt(2*M_PI*sigma[i])) *
                    (exp(0-0.5*(y[k]-nu[j])*(y[k]-nu[j])/tau[j])/sqrt(2*M_PI*tau[j]) ) / temp;
            }
        }
    }
    
    Rcpp::List rtn = Rcpp::List::create(x, y, zxy, n, pi, mu, 
                                        sigma, nu, tau, g);
    
    return rtn;
}
