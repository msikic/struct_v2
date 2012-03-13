# include "struct.h"

/***************************************/
/***************************************/
/***************************************/

double F (double **X, int * x_type, int NX,
	  double **Y, int * y_type, int NY, double alpha) {
    
    int i,j,m;
    double aux, sum, f;
    double alpha_sq = alpha*alpha;
    double prefactor;
    double step = NR_POINTS/MAX_EXP_VALUE;
    f = 0;
    for (i=0; i < NX; i++ ) {	
	for (j=0; j < NY; j++) {

	    sum = 0;
	    for (m=0; m<3; m++) {
		aux = X[i][m] - Y[j][m];
		sum += aux*aux;
	    }
	   
	    if (  (aux = sum/alpha_sq) <= MAX_EXP_VALUE ) {
		/* type, defined in struct_geometry.h,
		   can be HELIX or PARALLEL (beta-strand, in this implementation)*/
		prefactor = (x_type[i]&y_type[j]) ? REWARD : PENALTY;
		f -= prefactor* exp_table[(int)(aux*step)];
	    }
	}
    }
    
    return f;
}

/***************************************/
/***************************************/
/***************************************/
int F_moments (double **x, int *x_type, int NX,
	       double **y, int *y_type, int NY, double alpha,
	       double *avg_ptr, double *avg_sq_ptr,
	       double *std_ptr) {

    double avg, avg_sq, stdev, aux;
    int factor;
    int x_ctr, y_ctr;
    double expected_sq_grid_lookup (double **x, int * x_type, int NX,
				    double **y, int * y_type, int NY,
				    double alpha);
    
    aux = alpha*alpha/4;
    avg = -aux*(1.0 - exp(-1/aux));
    factor = 0;
    for (x_ctr=0; x_ctr < NX; x_ctr++) {
	for (y_ctr=0; y_ctr < NY; y_ctr++) {
	    factor += (x_type[x_ctr]&y_type[y_ctr]) ? REWARD : PENALTY;
	}
    }
    avg *= factor;
    stdev = 0;
    avg_sq = expected_sq_grid_lookup (y, y_type, NY,
				      x, x_type, NX, alpha);

    if ( avg_sq < 0.0 ||   avg_sq - avg*avg < 0.0) {
	//printf ( " %8.3lf   %8.3lf   %8.3lf \n", avg, avg_sq, stdev);
	//fprintf (stderr, "Error in lookup table\n");
	//return 1;
	*avg_ptr    = 0.0;
	*avg_sq_ptr = 0.0;
	*std_ptr    = 0.0;
	return 0;
   }

    stdev = sqrt (avg_sq - avg*avg);
    
    *avg_ptr = avg;
    *avg_sq_ptr = avg_sq;
    *std_ptr = stdev;


    //exit (1);
    
    return 0;
}


/**********************************************************/
/**********************************************************/
int monte_carlo (double alpha,
		 double **x, int * x_type, int NX,
		 double **y, int * y_type, int NY,
		 double  *q_best, double *F_best_ptr) {

    int no_MC_steps = 50;
    int step;
    double q[4], exp_s[4], q_new[4];
    double **R;
    double **x_rotated;
    double F_old, F_min, F_new;
    double T = 1.0;
    
    srand48(time(0));

    if ( !(R=dmatrix (3,3))) return 1;
    if ( !(x_rotated=dmatrix (NX,3))) return 1;

    memcpy (q, q_best, 4*sizeof(double) );
    quat_to_R (q, R);
    rotate (x_rotated, NX, R, x);
    F_min = F_old = F( y, y_type, NY, x_rotated, x_type, NX, alpha);
    
    //printf ( "int  F_min %8.3lf   ", F_min);
    // vec_out (q_best, 4, "q_init");
	
    for (step=0; step < no_MC_steps; step++) {

	/* make  new q */
	random_q (exp_s, M_PI/300);
	multiply (exp_s, q, 0, q_new);
	
	
	/* evaluate new F */
	quat_to_R (q_new, R);
	rotate (x_rotated, NX, R, x);
	F_new = F( y, y_type, NY, x_rotated, x_type, NX, alpha);
	
	if ( F_new < F_old ||  exp ((F_old-F_new)/T) < drand48() ) {
	    /* accepted */
	    F_old = F_new;
	    memcpy (q, q_new, 4*sizeof(double) );
	    if ( F_min > F_new ) {
		F_min = F_new;
		memcpy (q_best, q, 4*sizeof(double) );
		//printf ( ">>>> step %4d  F_min %8.3lf  ", step, F_min);
		//vec_out (q_best, 4, "q_current");
	    }
	}
	
    }
    
    *F_best_ptr = F_min;

    free_dmatrix (R);
    free_dmatrix (x_rotated);

    
    
    return 0;
}

/**********************************************************/
/**********************************************************/
/**********************************************************/

int gradient(int first_call, double **x,
	     double **x2, int * x_type,  int NX,
	     double **y, int * y_type, int NY,
	     double alpha, double q[4], double grad[4]) {
    
    int i,j,k,m, n, p;
    double sum, f;
    double alpha_sq = alpha*alpha;
    double scalar2;
    static double **prefactor = NULL;
    static double ***cross = NULL, **scalar = NULL;
    double q_cross, qx, qy, blah;
    double exp_part;
    double dot;
    double aux, step = NR_POINTS/MAX_EXP_VALUE;
    
    if ( first_call ) {
	if (prefactor) free_dmatrix (prefactor);
	if (scalar)    free_dmatrix (scalar);
	if (cross)     free_d3matrix (cross);
	prefactor = NULL;
	scalar    = NULL;
	cross     = NULL;
	if  ( ! NX ) return 0; /* shutdown */
	    
	if ( ! (prefactor = dmatrix (NX, NY) )) exit (1);
	if ( ! (scalar    = dmatrix (NX, NY) )) exit (1);
	if ( ! (cross     = d3matrix (NX, NY, 3) )) exit (1);
	
    }
    for (i=0; i < NX; i++ ) {
	for (j=0; j < NY; j++) {
	    prefactor[i][j] = (x_type[i]&y_type[j]) ? REWARD : PENALTY;

	    scalar[i][j] = 0.0;
	    for (m=0; m<3; m++ ) {
		n = (m+1)%3;
		p = (m+2)%3;
		cross[i][j][m]  = x[i][n]*y[j][p];
		cross[i][j][m] -= x[i][p]*y[j][n] ;
		scalar[i][j]   += x[i][m]*y[j][m];
	    }
	}
    }
    
    blah = q[0]*q[0];
    for (m=1; m<=3; m++ ) blah -= q[m]*q[m];
	
    for (k=0; k<4; k++ ) grad[k] = 0.0;
    
    
    f = 0;
    for (i=0; i < NX; i++ ) {
	for (j=0; j < NY; j++) {

	    q_cross = 0;
	    qx = 0;
	    qy = 0;

	    for (m=0; m<3; m++ ) {
		q_cross  += q[m+1]*cross[i][j][m];
		qx       += q[m+1]*x[i][m];
		qy       += q[m+1]*y[j][m];
	    }
	    double q_norm = 0;
	    for (m=0; m<4; m++ ) q_norm += q[m]*q[m];

	    
	    /** if the exp big enough do not calculate  exp et all
		-- continue bcs the contribution to the gradient is 0 */
	    dot = 0;
	    for (k=0; k<3; k++ ) dot += x2[i][k]*y[j][k];
	    if (dot > 1) dot = 1.0;
	    
	    sum = 2*(1-dot);
	    if (  (aux = sum/alpha_sq) > MAX_EXP_VALUE ) continue;

	    
	    exp_part =  2*prefactor[i][j]/alpha_sq* exp_table[(int)(aux*step)];

	    for (k=0; k<4; k++ ) {
		if ( k== 0 ) {
		    dot =  q[0]*scalar[i][j] + q_cross;
		} else {
		    dot = -q[k]*scalar[i][j] + y[j][k-1]*qx + x[i][k-1]*qy
			+ q[0]*cross[i][j][k-1];
		}
		grad[k] += exp_part*dot ;
	    }
	    
	    
	}
    }
   
    /* trick to keep q normalized */
    scalar2  = 0;
    for (k=0; k<4; k++ ) scalar2 += q[k]*grad[k];
    /* Lagrange multiplier*/
    for (k=0; k<4; k++ ) grad[k] -= scalar2*q[k]; 
    
    
    return 0;

}

/**********************************************************/
/**********************************************************/
/**********************************************************/

int gradient_descent (int first_call, double alpha,
		      double **x, int * x_type, int NX,
		      double **y, int * y_type, int NY,
		      double  *q_best, double *F_best_ptr) {

    double q[4] = {0.0},  q_old[4] = {0.0}, grad[4] = {0.0};
    double F_old, F_new =  -1;
    double **x2, norm;
    double step_size = options.grad_step_size;
    double **R;
    int i, done;
    int step = 0, max_step = options.grad_max_step;
    int gradient(int first_call, double **X,
		 double **X2, int * x_type, int NX,
		 double **Y, int * y_type, int NY,
		 double alpha, double q[4], double grad[4]);

    if ( first_call && !NX ) { /* shutdown */
	gradient (  first_call, NULL, NULL, NULL, 0,
		    NULL, NULL, 0, 0.0, NULL, NULL);
	return 0;
    }
    
    if ( ! (R=dmatrix(3,3) )) return 1; /* compiler is bugging me otherwise */
    if ( ! (x2 = dmatrix (NX,3) )) return 1;
    
    memcpy (q, q_best, 4*sizeof(double));
    quat_to_R (q, R);
    rotate (x2, NX, R, x);
    F_old =  F( y, y_type, NY, 
		   x2, x_type, NX, alpha);
    
    if ( F_old > 0 ) exit(1);
	    
    
    for (i=0; i<4; i++) q_old[i] = 2;
    done = 0;
    for (step=1; step<=max_step && !done; step++) {

    
	gradient (  first_call && (step==1),
		    x,  x2,  x_type, NX, y, y_type, NY, 
		    alpha, q, grad);

	if ( isnan (grad[0]) ) { 
	    fprintf ( stderr, "in gradient descent: \n");
	    fprintf ( stderr, "step %3d   \n", step);
	    fprintf ( stderr," grad: %6.3lf  %6.3lf   %6.3lf  %6.3lf \n",
		    grad[0], grad[1], grad[2], grad[3] );
	    fprintf ( stderr," q    : %6.3lf  %6.3lf   %6.3lf  %6.3lf \n",
		    q[0], q[1], q[2], q[3] );
	    fprintf ( stderr, "F : %8.3lf   %8.3lf \n\n",
		      F_new, F_old );
	    return 17;
	}
	/* add gradient to the old q */
	norm = 0.0;
	for (i=0; i<4; i++) {
	    q[i] += step_size*grad[i];
	    norm += q[i]*q[i];
	}
	/* q-->R */
	quat_to_R (q, R);
		
	/* rotate X */
	rotate (x2, NX, R, x);
		
	/* calculate F */
	F_new =  F( y, y_type, NY, x2, x_type, NX, alpha);
 

	/* STOPPING CRITERION */
# if 1
	if (F_new < F_old){
	    double diff = fabs((F_old-F_new)/F_old );
	    if ( diff < options.grad_stop_tol) {
		done = 1;
		//printf ("reached precision\n");
	    } 
	    F_old = F_new; 
	    memcpy (q_old, q, 4*sizeof(double));

	} else {
	    step_size /= 2;
	    //printf ("changing step size: %8.2le \n", step_size);
	    if ( step_size < 1.e-6 ) {
	      done = 1;
	    }
	}
# else
	if ( F_old < 0.0 && F_new < F_old   ){
	    diff = fabs((F_old-F_new)/F_old );
	    if ( diff < options.grad_stop_tol) {
		done = 1;
		//printf ("reached precision\n");
	    } 
	}
	F_old = F_new; 
	memcpy (q_old, q, 4*sizeof(double));

# endif



    }

# if 0
    if ( step== max_step+1) {
	printf ("stopped by maxstep\n");
	exit (0);
    } else {
	printf ("stopped by tolerance: F_old: %8.3lf  F_new: %8.3lf  \n", F_old, F_new);
    }
# endif    
    memcpy (q_best, q, 4*sizeof(double));
    *F_best_ptr = F_new;
    
    free_dmatrix (R);
    free_dmatrix (x2);

    return 0;
    
}
    
/**********************************************************/
/**********************************************************/
/**********************************************************/

double expected_sq_grid_lookup (double **x,  int * x_type, int NX,
				double **y, int * y_type, int NY, double delta) {

    int i, j, k, l;
    //int s_i, s_j, s_k, s_l;
    int prefactor;
    int pos_ct, neg_ct, skipped, x_hel, y_hel;
    double cos_beta;
    double integral, total;
    double **cos_alpha, tmp_alpha;
    double total_pos, total_neg;

    if ( ! (cos_alpha= dmatrix(NY, NY) ) ) return 1;

    for (j=0; j<NY; j++) {
	cos_alpha[j][j] = 1.0;
	for (l=j+1; l<NY; l++) {
	    unnorm_dot (y[j], y[l], &tmp_alpha);
	    cos_alpha[j][l] = cos_alpha[l][j] = tmp_alpha;
	}
    }

    x_hel = 0;
    for (i=0; i<NX; i++) 
	if (x_type[i]==HELIX)  x_hel++;
    y_hel = 0;
    for (i=0; i<NY; i++) 
	if (y_type[i]==HELIX)  y_hel++;
   
    total = 0.0;
    total_pos = 0.0;
    total_neg = 0.0;
    pos_ct = neg_ct = skipped = 0;
    for (i=0; i<NX; i++) {
      //s_i = x_type[i]==HELIX ? 1:-1;
	
	for (k=0; k<NX; k++) {
	  //s_k = x_type[k]==HELIX ? 1:-1;
	    
	    unnorm_dot (x[i], x[k], &cos_beta);
	    
	    for (j=0; j<NY; j++) {
	      //s_j = y_type[j]==HELIX ? 1:-1;
		
		for (l=0; l<NY; l++) {
		  //s_l = y_type[l]==HELIX ? 1:-1;

		    integral  = lookup(cos_alpha[j][l], cos_beta);
		    
		 
		    if ( fabs(integral) < 1.e-7 ) {
		    	skipped ++;
		    	continue;
		    }
		    /***/
		    //prefactor  = s_i*s_j*s_k*s_l;
		    prefactor =
		      ((x_type[i]&y_type[j]) ? REWARD : PENALTY) *
		      ((x_type[k]&y_type[l]) ? REWARD : PENALTY);
		    total += prefactor*integral;
		    if (0) {
			if ( prefactor > 0) {
			    total_pos += integral;
			    pos_ct ++;
			} else {
			    total_neg += integral;
			    neg_ct ++;
			}

			// printf ( " %8.3lf  %8.3lf  %2d %2d %2d %2d   %8.2le   %8.2le \n",
			//	     cos_alpha[j][l], cos_beta,
			//	     s_i, s_j, s_k, s_l, integral, total);
		    } 

		}
	    }
	}
    }

    if (0) {
	/* there is some numerical problem here, cf  2k7jA.db 1us0A.db  */
	printf ( " NX: %3d  NY: %3d  loop size: %8d \n", NX, NY, NX*NX*NY*NY);
	printf ( " X helices: %3d     Y helices: %3d \n", x_hel, y_hel);
	printf ( "    pos: %8.2lf   neg: %8.2lf \n", total_pos, total_neg);
	printf ( " pos ct: %8d   neg ct: %8d   skipped: %8d \n", pos_ct, neg_ct, skipped);
	printf ( " %12.7lf \n", total);

	exit (0);
    }
    
    free_dmatrix (cos_alpha);
    
    return total;
}
