# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include <string.h>
//# include "struct.h"

# define DELTA 0.5
# define NR_POINTS 100



int main ( int argc, char * argv[]) {

    int nr_points = NR_POINTS;
    int a, b, a_steps, b_steps;
    double alpha_step;
    double beta_step;
    double cos_alpha, cos_beta;
    double avg_sq = 0.0;
    double delta  = DELTA;
    double x[2][3] = {{0.0}};
    double y[2][3] = {{0.0}};
    double moment (int moment_n, double x[][3], double y[][3], double delta);

    if ( argc < 3 ) {
	fprintf (stderr, "Usage: %s <nr points> <delta>\n", argv[0]);
	exit (1);
	
    } else {
	nr_points = atoi (argv[1] );
	delta     = atof (argv[2] );
    }

    a_steps = b_steps = nr_points;
    alpha_step = 2.0/a_steps;
    beta_step  = 2.0/b_steps;
   

    printf ("%% points %6d \n", nr_points);
    printf ("%% delta  %6.3lf \n", delta);
    
    for (a =0; a<= a_steps; a++) {
 	
	cos_alpha = 1 -  a*alpha_step;
	memset (y[0], 0.0, 2*3*sizeof(double));
	y[0][2] = 1;
	y[1][2] = cos_alpha;
	y[1][0] = sqrt(1-cos_alpha*cos_alpha);
	    
	for (b=a; b<= b_steps; b++) {
	    
	    cos_beta =  1 -  b*beta_step;
	    memset (x[0], 0.0, 2*3*sizeof(double));
	    
	    x[0][2] = 1;
	    x[1][2] = cos_beta;
	    x[1][0] = sqrt(1-cos_beta*cos_beta);
    
	    avg_sq =  moment (2, x, y, delta);
	    if ( avg_sq > 1.e-7) {
		printf (" %3d  %3d   %8.3lf  %8.3lf    %14.10lf  \n",
		           a,   b, cos_alpha, cos_beta, avg_sq);
		fflush (stdout);
	    }
	}
    }


    return 0;
}


/******************************************/
double moment ( int moment_n, double x[][3], double y[][3], double delta) {

    int n_steps     = 100;
    int phi_n_steps = 100;
    int i, j, n;
    int k, t, p;
    double value = 0.0;
    double delta_sq = delta*delta;
    double kappa, theta, phi;
    double kappa_step, theta_step, phi_step;
    double cosk, sink, costh, sinth;
    double surface_element, phi_contributn;
    
    double ATA     [4][4] = {{0.0}};
    double prev_ATA[4][4] = {{0.0}};
    double ATA_sum [4][4] = {{0.0}};
    double a[3] = {0.0}, b[3] = {0.0};
    double q[4] = {0.0};
    /**************/
    int construct_ATA (double ATA[4][4], double a[3],double b[3]);
    int add_matrices  (double matrix1[4][4],double matrix2[4][4],double result[4][4]);
    double braket     (double ATA[4][4], double q[4]);
    
    /* sum of A^TA matrices */
    for (n=0; n<moment_n; n++ ) {
	
	for (i=0; i<3; i++ ) {
	    a[i] = y[n][i] + x[n][i];
	    b[i] = y[n][i] - x[n][i];
	}
	construct_ATA (ATA, a, b);
	add_matrices (prev_ATA, ATA, ATA_sum);
	memcpy (prev_ATA[0], ATA_sum[0], 4*4*sizeof(double));
    }

    if (0) { /*check */

	q[0] = 1;
	for (i=0; i<4; i++ ) {
	    for (j=0; j<4; j++ ) {
		printf ("%8.3lf ", ATA_sum[i][j]);
	    }
	    printf ("\n");
	}
	printf ("\n");
	double aux, sum_sq = 0.0;;
	for (i=0; i<4; i++ ) {
	    aux = x[1][i]-y[1][i];
	    sum_sq +=  aux*aux;
	}
	printf ( " 00: %8.4le   %8.4le  %8.4le \n",
		 ATA_sum[0][0], braket(ATA_sum, q), sum_sq);
	exit(1);
    }

    
    value = 0.0;
    kappa_step = M_PI/n_steps;
    theta_step = M_PI/n_steps;
    phi_step   = 2*M_PI/phi_n_steps;
    
    /* for angles kappa, theta phi */
    for (k=1; k<=n_steps; k++) {
	kappa = ( (double)k-0.5)*kappa_step;
	cosk = cos(kappa);
	sink = sin(kappa);
	
	q[0] = cosk;
	
	for (t=1; t<=n_steps; t++) {
	    theta = ( (double)t-0.5)*theta_step;
	    sinth = sin (theta);
	    costh = cos (theta);

	    q[3] = sink*costh;
	    surface_element = sink*sink*sinth;
	    phi_contributn  = 0;
	    
	    for (p=1; p<=phi_n_steps; p++) {
		phi  = ((double)p-0.5)*phi_step;

		q[1] = sink*sinth*cos(phi); 
		q[2] = sink*sinth*sin(phi); 

		phi_contributn += exp(-braket(ATA_sum, q)/delta_sq);
		//value += sink*sink*sinth;
		//printf ( "%3d  %3d  %3d  %8.4le  %8.2le   --  %8.4le \n", k, t, p,
		//braket(ATA, q), exp(-braket(ATA, q)/delta_sq), phi_contributn);
	    }

	    value += phi_contributn*surface_element;
	}
    }
    
    value *= kappa_step*theta_step*phi_step; /*integration steps*/

    value /= 2*M_PI*M_PI; /* the function should return the average*/
    
    return value;
    
}

/************************************************/
double braket  (double ATA[4][4],double  q[4]){
    int i,j;
    double value = 0.0;

    for (i=0; i<4; i++ ) {
	for (j=0; j<4; j++ ) {
	    value += q[i]*ATA[i][j]*q[j];
	}
    }
   
    return value;
}
 
/************************************************/
int construct_ATA (double ATA[4][4], double a[3], double  b[3]){

    int i,j,k;
    double A[4][4] = {{ 0.0, -b[0], -b[1], -b[2]},
		      {b[0],   0.0, -a[2],  a[1]},
		      {b[1],  a[2],   0.0, -a[0]},
		      {b[2], -a[1],  a[0],   0.0}};
    
    for (i=0; i<4; i++ ) {
	for (j=0; j<4; j++ ) {
	    ATA[i][j] = 0.0;
	    for (k=0; k<4; k++ ) {
		ATA[i][j] += A[k][i]*A[k][j];
	    }
	}
    }
     
    return 0;
}

/************************************************/
int add_matrices  (double matrix1[4][4],double matrix2[4][4],
		   double result[4][4]){
    int i,j;

    for (i=0; i<4; i++ ) {
	for (j=0; j<4; j++ ) {
	    result[i][j] = matrix1[i][j] + matrix2[i][j];
	}
    }
    return 0;

}

    
