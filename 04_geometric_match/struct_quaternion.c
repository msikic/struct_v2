# include "struct.h"

# define SVD
/**********************************************************/
/**********************************************************/
/**********************************************************/

int quat_to_R ( double quat[4], double **R ){

    
    double q0 = quat[0];
    double q1 = quat[1];
    double q2 = quat[2];
    double q3 = quat[3];


    R[0][0] = 1 -2*q2*q2 -2*q3*q3;
    R[1][1] = 1 -2*q3*q3 -2*q1*q1 ;
    R[2][2] = 1 -2*q1*q1 -2*q2*q2 ;

    R[0][1] = 2*q1*q2 - 2*q0*q3;
    R[1][0] = 2*q1*q2 + 2*q0*q3;
    
    R[0][2] = 2*q1*q3 + 2*q0*q2;
    R[2][0] = 2*q1*q3 - 2*q0*q2;
    
    R[1][2] = 2*q2*q3 - 2*q0*q1;
    R[2][1] = 2*q2*q3 + 2*q0*q1;

    return 0;

}
/**********************************************************/
/**********************************************************/
/**********************************************************/
double q_dotprod (double *quat1, double *quat2) {

    int i;
    double product = 0;
    
     for (i=0; i<=3; i++) {
	product +=  quat1[i]*quat2[i];
    }
     
    return fabs(product);
}

/**********************************************************/
/**********************************************************/
/**********************************************************/
int multiply (double *quat1, double *quat2_in,
	      int conjugate,  double *product) {

    double quat2[4];
    int i;
    
    memset ( product, 0,  4*sizeof (double) );
    memcpy ( quat2, quat2_in,  4*sizeof (double) );
    
    if ( conjugate ) {
	for (i=1; i<=3; i++) quat2[i] = -quat2[i];
    }
    product[0] = quat1[0]*quat2[0];
    for (i=1; i<=3; i++) product[0] -= quat1[i]*quat2[i];

    for (i=1; i<=3; i++) product[i] +=  quat1[0]*quat2[i];
    
    for (i=1; i<=3; i++) product[i] +=  quat1[i]*quat2[0];

    /* the vector mult part */
    product[1] +=  quat1[2]*quat2[3] -  quat1[3]*quat2[2];
    product[2] +=  quat1[3]*quat2[1] -  quat1[1]*quat2[3];
    product[3] +=  quat1[1]*quat2[2] -  quat1[2]*quat2[1];
    
    return 0;
}
/**********************************************************/
/**********************************************************/
/**********************************************************/
/************************************/
int  qmap_bisec (double *x0, double *x1, double *y0, double *y1, double * quat){

    
    double q1[4], q2[4];
    double v[3], cosine = 0, theta =0, sine;
    double bisec_x[3], bisec_y[3];
    double norm_x[3], norm_y[3];
    double bisec_x_prime[3];
    double norm;
    int i,j;
    
    int  normalized_cross (double *x, double *y, double * v, double *norm_ptr);
    int  unnorm_dot (double *x, double *y, double * dot);

    /* 0) find normals and bisectors */
    if ( normalized_cross (x0, x1, norm_x, NULL) )return 1;
    if ( normalized_cross (y0, y1, norm_y, NULL) )return 1;

    norm = 0;
    for (i=0; i<3; i++ ) {
	bisec_x[i] = (x0[i] + x1[i])/2;
	norm += bisec_x[i]*bisec_x[i];
    }
    norm += sqrt(norm);
    if (norm) for (i=0; i<3; i++ ) bisec_x[i] /= norm;
    
    norm = 0;
    for (i=0; i<3; i++ ) {
	bisec_y[i] = (y0[i] + y1[i])/2;
	norm += bisec_y[i]*bisec_y[i];
    }
    norm += sqrt(norm);
    if (norm) for (i=0; i<3; i++ ) bisec_y[i] /= norm;
    
    /* 1) find any quat that will match the normals */
    /* check that it is not the same vector already: */
    unnorm_dot (norm_x, norm_y, &cosine);
    if ( 1-cosine > 0.001 ) {
	
	double ** R;
	
	if ( normalized_cross (norm_x, norm_y, v,  NULL) )return 1;

	theta = acos (cosine);

	q1[0] = cos (theta/2);
	sine = sin(theta/2);
	for (i=0; i<3; i++ ) {
	    q1[i+1] = sine*v[i];
	}

	if ( ! (R=dmatrix(3,3) ) ) return 1; /* compiler is bugging me otherwise */
	quat_to_R (q1, R);

	for (i=0; i<3; i++ ) {
	    bisec_x_prime[i] = 0;
	    for (j=0; j<3; j++ ) {
		bisec_x_prime[i] += R[i][j]*bisec_x[j];
	    }
	}
	free_dmatrix (R);
	
    } else {
	memcpy (bisec_x_prime, bisec_x, 3*sizeof(double) );
	memset (q1, 0,  4*sizeof(double) );
	q1[0] = 1.0;
    }

    
    /* 2) find quat thru norm_y which maps bisectors one onto another */
   
    unnorm_dot (bisec_x_prime, bisec_y, &cosine);
    /* determining theta: angle btw planes = angle btw the normals */
    if ( 1-cosine > 0.001 ) {
	
	
	theta = acos (cosine);
	
	/*still need to determine the sign of rotation */
	if (normalized_cross (bisec_x_prime, bisec_y, v, NULL)) return 1;
	
	unnorm_dot (v, norm_y, &cosine);
	if ( cosine > 0.0 ) { 
	    sine =  sin(theta/2);
	} else {
	    sine =  -sin(theta/2);
	}
	for (i=0; i<3; i++ ) {
	    v[i]    = norm_y[i];
	}
	q2[0] = cos (theta/2);
	for (i=0; i<3; i++ ) {
	    q2[i+1] = sine*v[i];
	}

    } else {
	memset (q2, 0,  4*sizeof(double) );
	q2[0] = 1.0;
    }

    /* multiply the two quats */
    multiply (q2, q1, 0, quat);
    
    return 0;
       
}
/************************************/
/************************************/
int  qmap (double *x0, double *x1, double *y0, double *y1, double * quat){

    double q1[4], q2[4];
    double v[3], v2[3], cosine = 0, theta =0, sine;
    double x_prime[3];
    int i,j;
    
    int  normalized_cross (double *x, double *y, double * v, double *norm_ptr);
    int  unnorm_dot (double *x, double *y, double * dot);
    
   
    /* 1) find any quat for  x0 -->  y0 */

    /* check that it is not the same vector already: */
    unnorm_dot (x0, y0, &cosine);
    if ( 1-cosine > 0.001 ) {
	
	double ** R;
	
	if ( normalized_cross (x0, y0, v, NULL) )return 1;

	theta = acos (cosine);

	q1[0] = cos (theta/2);
	sine = sin(theta/2);
	for (i=0; i<3; i++ ) {
	    q1[i+1] = sine*v[i];
	}

	if ( ! (R=dmatrix(3,3) ) ) return 1; /* compiler is bugging me otherwise */
	quat_to_R (q1, R);

	for (i=0; i<3; i++ ) {
	    x_prime[i] = 0;
	    for (j=0; j<3; j++ ) {
		x_prime[i] += R[i][j]*x1[j];
	    }
	}
	free_dmatrix (R);
	
    } else {
	memcpy (x_prime, x1, 3*sizeof(double) );
	memset (q1, 0,  4*sizeof(double) );
	q1[0] = 1.0;
    }

    
    /* 2) find quat thru y0 which maps
       x' to the plane <y0,y1> */
   
    unnorm_dot (x_prime, y1, &cosine);
    /* determining theta: angle btw planes = angle btw the normals */
    if ( 1-cosine > 0.001 ) {
	if ( normalized_cross (x_prime, y0,v, NULL)) return 1;
	if ( normalized_cross (y1, y0, v2, NULL)) return 1;
	unnorm_dot (v, v2, &cosine);
	
	theta = acos (cosine);
	
	/*still need to determine the sign
	  of rotation */
	if (normalized_cross (x_prime, y1,v, NULL)) return 1;
	unnorm_dot (v, y0, &cosine);
	if ( cosine > 0.0 ) { 
	    sine =  sin(theta/2);
	} else {
	    sine =  -sin(theta/2);
	}
	for (i=0; i<3; i++ ) {
	    v[i]    = y0[i];
	}
	q2[0] = cos (theta/2);
	for (i=0; i<3; i++ ) {
	    q2[i+1] = sine*v[i];
	}

    } else {
	memset (q2, 0,  4*sizeof(double) );
	q2[0] = 1.0;
    }

    /* multiply the two quats */
    multiply (q2, q1, 0, quat);
    

    return 0;
       
}
/************************************/

int random_q ( double exp_s[4],  double theta_step ) {
    
    double abs_s;
    double s[3];
    double theta, sin_theta;
    int i;
    
    /* pick s from a cube of side length 2N */
    abs_s = 0.0;
    for (i=0; i<3; i++) {
	s[i] = (0.5 - drand48())*2;
	abs_s += s[i]*s[i];
    }
    abs_s = sqrt (abs_s);
    for (i=0; i<3; i++) {
	s[i] /= abs_s;
    }
    //vec_out (s, 3, "s");
	
    /* calculate exponent (? this is not the exponent ...)*/
    theta = (0.5 - drand48())*2*theta_step;
    exp_s[0]  = cos(theta);
    sin_theta = sin(theta);
    for (i=1; i<4; i++) {
	exp_s[i] = sin_theta*s[i-1];
    }
    //vec_out (exp_s, 4, "exp_s");

    return 0;
}

/************************************/
int print_quat ( FILE *fptr, double * quat, char *name ) {

    fprintf (fptr, " %10s  %6.3lf  %6.3lf  %6.3lf  %6.3lf \n", name,
	     quat[0], quat[1],quat[2], quat[3]) ;
    return 0;
}


double qnorm (double *quat) {
    double norm = 0;
    int i;
    
    for (i=0; i<4; i++) {
	norm += quat[i]*quat[i];
    }
    return sqrt(norm);
}
