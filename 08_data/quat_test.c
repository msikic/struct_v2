# include <stdio.h>
# include <stdlib.h>
# include <math.h>
# include <time.h>
# include <string.h>

/******************************************/
/******************************************/

int main ( int argc, char * argv[]) {

    int test_number = 0;
    double alpha, beta;
    double delta = 0.3;
    double x[2][3] = {{0.0}};
    double y[2][3] = {{0.0}};
    double moment (int moment_n, double x[][3], double y[][3], double delta);
    double avg= 0.0, avg_sq = 0.0, avg_sq_sw = 0.0;
    double avg_exact= 0.0;

    if ( argc == 2 ) test_number = atoi(argv[1]); 

    switch (test_number) {
	
    case 0:

	x[0][2] = 1;
	y[0][2] = 1;


	avg = moment (1, x, y, delta);
	avg_exact = delta*delta/4*(1.0 - exp (-4.0/(delta*delta) ) );
	printf ("alpha: %8.3lf    beta: %8.3lf   avg:  %8.10lf    exact: %8.10lf  diff: %8.2le\n",
		alpha, beta, avg, avg_exact, fabs(avg-avg_exact)/avg_exact);
	break;
	
    case 1:
    
	alpha = M_PI/3;
	beta  = M_PI/6;

	x[0][2] = 1;
	x[1][2] = cos (beta);
	x[1][0] = sin (beta);
    
	y[0][2] = 1;
	y[1][2] = cos (alpha);
	y[1][0] = sin (alpha);


	avg_sq = moment (2, x, y, delta);
	avg_sq_sw = moment (2, y, x, delta);
	printf ("alpha: %8.3lf    beta: %8.3lf   avg_sq:  %8.10lf   switched: %8.10lf \n",
		alpha, beta, avg_sq, avg_sq_sw);
	break;

    case 2:
    
	alpha = beta  = M_PI/3;

	x[0][2] = 1;
	x[1][2] = cos (beta);
	x[1][0] = sin (beta);
    
	y[0][2] = 1;
	y[1][2] = cos (alpha);
	y[1][0] = sin (alpha);


	avg_sq = moment (2, x, y, delta);
	printf ("alpha: %8.3lf    beta: %8.3lf   avg_sq:  %8.10lf  \n",
		alpha, beta, avg_sq);
	break;

     case 3:
    
	alpha = beta  = M_PI/10;

	x[0][2] = 1;
	x[1][2] = cos (beta);
	x[1][0] = sin (beta);
    
	y[0][2] = 1;
	y[1][2] = -x[1][2];
	y[1][0] = -x[1][0];


	avg_sq = moment (2, x, y, delta);
	printf ("alpha: %8.3lf    beta: %8.3lf   avg_sq:  %8.10lf ( %8.3le)\n",
		alpha, beta, avg_sq, avg_sq);
	break;

    default:
	fprintf ( stderr, "Unrecognized test number: %s\n", argv[1]); 

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

    
