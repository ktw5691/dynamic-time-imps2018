#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>


/**
 * This function takes a double and gives back a double
 * It computes the logistic function (i.e. the inverse of the logit link function)
 * @param x, the double value e.g. a normally distributed number
 * @return logistic(x), the double value e.g. a number between 0 and 1
 */
double mathfunction_logistic(const double x){
	double value = 1.0/(1.0 + exp(-x));
	return value;
}

/**
 * This function takes a gsl_vector and modifies its second argument (another gsl_vector)
 * It computes the softmax function (e.g. for multinomial logistic regression)
 * @param x, vector of double values e.g. a vector of normally distributed numbers
 * @param result, softmax(x), e.g. a vector of numbers between 0 and 1 that sum to 1
 */
void mathfunction_softmax(const gsl_vector *x, gsl_vector *result){
	/* Elementwise exponentiation */
	size_t index=0;
	for(index=0; index < x->size; index++){
		gsl_vector_set(result, index, exp(gsl_vector_get(x, index)));
	}
	
	/* Sum for the scaling coeficient */
	double scale = 0.0;
	for(index=0; index < x->size; index++){
		scale += gsl_vector_get(result, index);
	}
	
	/* Multiply all elements of result by 1/scale */
	gsl_blas_dscal(1/scale, result);
}


void function_dynam(const double tstart, const double tend, size_t regime, const gsl_vector *xstart, double *param, size_t n_gparam, const gsl_vector *co_variate, void (*g)(double, size_t, const gsl_vector *, double *, size_t, const gsl_vector *, gsl_vector *), gsl_vector *x_tend){
	gsl_matrix *Amatrix = gsl_matrix_calloc(2, 2);


	gsl_blas_dgemv(CblasNoTrans, 1.0, Amatrix, xstart, 0.0, x_tend);

	gsl_matrix_free(Amatrix);

	gsl_vector *intVector = gsl_vector_calloc(2);


	gsl_vector_add(x_tend, intVector);
	gsl_vector_free(intVector);

}


void function_jacob_dynam(const double tstart, const double tend, size_t regime, const gsl_vector *xstart, double *param, size_t num_func_param, const gsl_vector *co_variate, void (*g)(double, size_t, double *, const gsl_vector *, gsl_matrix *), gsl_matrix *Jx){

}



void function_measurement(size_t t, size_t regime, double *param, const gsl_vector *eta, const gsl_vector *co_variate, gsl_matrix *Ht, gsl_vector *y){


		gsl_vector *intVector = gsl_vector_calloc(6);

	gsl_matrix_set(Ht, 0, 0, 1);
	gsl_matrix_set(Ht, 1, 0, param[0]);
	gsl_matrix_set(Ht, 2, 0, param[1]);
	gsl_matrix_set(Ht, 3, 1, 1);
	gsl_matrix_set(Ht, 4, 1, param[2]);
	gsl_matrix_set(Ht, 5, 1, param[3]);

	gsl_vector_set(intVector, 0, param[4]);
	gsl_vector_set(intVector, 1, param[5]);
	gsl_vector_set(intVector, 2, param[6]);
	gsl_vector_set(intVector, 3, param[7]);
	gsl_vector_set(intVector, 4, param[8]);
	gsl_vector_set(intVector, 5, param[9]);
 
	gsl_blas_dgemv(CblasNoTrans, 1.0, Ht, eta, 0.0, y);
 
	gsl_vector_add(y, intVector);
 	gsl_vector_free(intVector);
 
}



void function_noise_cov(size_t t, size_t regime, double *param, gsl_matrix *y_noise_cov, gsl_matrix *eta_noise_cov){


	gsl_matrix_set(eta_noise_cov, 0, 0, param[10]);
	gsl_matrix_set(eta_noise_cov, 1, 0, param[11]);
	gsl_matrix_set(eta_noise_cov, 0, 1, param[11]);
	gsl_matrix_set(eta_noise_cov, 1, 1, param[12]);

	gsl_matrix_set(y_noise_cov, 0, 0, param[13]);
	gsl_matrix_set(y_noise_cov, 1, 1, param[14]);
	gsl_matrix_set(y_noise_cov, 2, 2, param[15]);
	gsl_matrix_set(y_noise_cov, 3, 3, param[16]);
	gsl_matrix_set(y_noise_cov, 4, 4, param[17]);
	gsl_matrix_set(y_noise_cov, 5, 5, param[18]);
 
}



void function_initial_condition(double *param, gsl_vector **co_variate, gsl_vector **pr_0, gsl_vector **eta_0, gsl_matrix **error_cov_0, size_t *index_sbj){
	
	gsl_vector *Pvector = gsl_vector_calloc(1);
	gsl_vector *Pintercept = gsl_vector_calloc(1);
	gsl_vector *Padd = gsl_vector_calloc(1);
	gsl_vector *Preset = gsl_vector_calloc(1);
	gsl_vector_set(Pvector, 0, 1);
	gsl_vector_add(Padd, Pvector);
	gsl_vector_add(Padd, Pintercept);
	gsl_vector_add(Preset, Pvector);
	gsl_vector_add(Preset, Pintercept);
	gsl_vector *eta_local = gsl_vector_calloc(2);
	size_t num_regime=pr_0[0]->size;
	size_t dim_latent_var=error_cov_0[0]->size1;
	size_t num_sbj=(eta_0[0]->size)/(dim_latent_var);
	size_t i;
	size_t regime;
	for(regime=0; regime < num_regime; regime++){
		for(i=0; i < num_sbj; i++){
			gsl_vector_set(eta_0[regime], i*dim_latent_var+0, gsl_vector_get(eta_local, 0));
			gsl_vector_set(eta_0[regime], i*dim_latent_var+1, gsl_vector_get(eta_local, 1));
			gsl_vector_set_zero(eta_local);
		}
	}
	for(i=0; i < num_sbj; i++){
		mathfunction_softmax(Padd, pr_0[i]);
	}
	gsl_vector_free(Pvector);
	gsl_vector_free(Pintercept);
	gsl_vector_free(Padd);
	gsl_vector_free(Preset);
	gsl_vector_free(eta_local);
}


/**
 * This function modifies some of the parameters so that it satisfies the model constraint.
 * Do not include parameters in noise_cov matrices 
 */
void function_transform(double *param){
}



void function_regime_switch(size_t t, size_t type, double *param, const gsl_vector *co_variate, gsl_matrix *regime_switch_mat){
	gsl_matrix_set_identity(regime_switch_mat);
}

