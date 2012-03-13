# include "struct.h"

/************************************/
/************************************/
/************************************/
/************************************/
/************************************/
int  unnorm_cross (double *x, double *y, double * v) {

    /* v is the output */
    int sum_to_zero (double *a, double *b);

    if (sum_to_zero (x, y) ) return 1;
    
    v[0] = x[1]*y[2] -  x[2]*y[1];
    v[1] = x[2]*y[0] -  x[0]*y[2]; 
    v[2] = x[0]*y[1] -  x[1]*y[0];

    return 0;
    
}
/************************************/
int  normalized_cross (double *x, double *y, double * v, double *norm_ptr) {

    /* v is the output */
    double norm = 0;
    int i;
    int sum_to_zero (double *a, double *b);

    if (sum_to_zero (x, y) ) return 1;
    
    v[0] = x[1]*y[2] -  x[2]*y[1];
    norm += v[0]*v[0];
    v[1] = x[2]*y[0] -  x[0]*y[2]; 
    norm += v[1]*v[1];
    v[2] = x[0]*y[1] -  x[1]*y[0];
    norm += v[2]*v[2];
    norm = sqrt (norm);

    for (i=0; i<3; i++ ) {
	v[i] /= norm;
    }

    if ( norm_ptr) *norm_ptr = norm;

    return 0;
    
}

/************************************/
int unnorm_dot (double *x, double *y, double * dot) {


    double cosine = 0;
    int i;
    
    for (i=0; i<3; i++ ) {
	cosine += x[i]*y[i];
    }

    if (cosine > 1.0 ) 
	cosine = 1.0; /* this should be numerical */
   
    *dot = cosine;
   
    return 0;
    
}

/***************************************/
int rotate(double **Ynew, int NY, double **R, double ** Y) {

    int n, i,j;

    for (n=0; n<NY; n++ ) {

	for (i=0; i<3; i++ ) {
	    Ynew[n][i] = 0;
	    for (j=0; j<3; j++ ) {
		Ynew[n][i] += R[i][j]*Y[n][j];
	    }
	}
    }

    return 0;
}

/***************************************/
/***************************************/
/***************************************/
/*****************************************/
int vec_out (double *vec, int dim,  char * name ) {
    int i;
    double norm = 0;
    
    printf (" %5s  ", name);
    for (i=0; i<dim; i++) {
	printf ("%6.3lf  ", vec[i]);
	norm +=  vec[i]* vec[i];
    }
    printf ("   norm: %7.4lf  ", norm);
    printf ("\n");
    return 0;
}


/************************************/
int  sum_to_zero (double *a, double *b ) {
    int i;
    double sum,aux;

    sum = 0;
    for (i=0; i<3; i++ ) {
	aux = a[i] + b[i];
	sum += aux*aux;
    }
    return (sum < 0.001);
}
