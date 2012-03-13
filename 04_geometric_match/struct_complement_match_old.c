# include "struct.h"

# define  TOP_RMSD 30
# define  BAD_RMSD 10.0
# define JACKFRUIT 8

/****************************************/
int complement_match (Representation* X_rep, Representation* Y_rep,
		      Map * map, int map_max,
		      int * map_ctr, int * map_best, int best_max, int parent_map){
			
    Penalty_parametrization penalty_params; /* for SW */
    double **x    = X_rep->full;
    int * x_type  = X_rep->full_type;
    int NX        = X_rep->N_full;
    double **y    = Y_rep->full;
    int * y_type  = Y_rep->full_type;
    int NY        = Y_rep->N_full;
    
    double F_effective = 0.0;
    double F_current;
    double q[4] = {0.0}, q_init[4] = {0.0};
    double **x_rotated = NULL;
    double **tr_x_rotated = NULL;
    double **R;
    double z_scr = 0.0, *z_best;
    double avg, avg_sq, stdev;
    double alpha = options.alpha;
    double rmsd, best_rmsd[TOP_RMSD];
    double **best_quat;
    double cutoff_rmsd = 3.0; /* <<<<<<<<<<<<<<<<< hardcoded */
    int *x_type_fudg, *y_type_fudg;
    int *anchor_x, *anchor_y, no_anchors;
    int no_top_rmsd = TOP_RMSD, chunk;
    int x_ctr, y_ctr, top_ctr;
    int **best_triple_x;
    int **best_triple_y;
    int x_triple[3], y_triple[3];
    int retval, done = 0;
    int best_ctr;
    int i, j;
    int t;
    int smaller;
    int my_map_ctr;
    int stored_new;
    int * x2y, map_unstable;
    //time_t  time_now, time_start;
    
    int cull_by_dna (Representation * X_rep, int *set_of_directions_x,
		 Representation * Y_rep, int *set_of_directions_y,
		     int set_size, Map *map, double cutoff_rmsd);
    int distance_of_nearest_approach (Representation * X_rep, int *set_of_directions_x,
				      Representation * Y_rep, int *set_of_directions_y,
				      int set_size,  double * rmsd_ptr);
    int same_hand_triple (Representation * X_rep, int *set_of_directions_x,
			  Representation * Y_rep, int *set_of_directions_y, int set_size);
    
    int find_map (Penalty_parametrization * params, Representation *X_rep,  Representation *Y_rep,
		  double ** R, double alpha, double * F_effective,  Map *map, 
		  int *anchor_x, int * anchor_y, int anchor_size );
    int find_next_triple (double **X, double **Y, 
			  int *x_type, int *y_type, int NX, int NY,
			  int *x_triple, int *y_triple);
    int gradient_descent (int first_call, double alpha,
			  double **x, int * x_type, int NX,
			  double **y, int * y_type, int NY,
			  double *q_best, double *F_best_ptr) ;
    int map_quality_metrics (Representation *X_rep, Representation *Y_rep,
			     double ** tr_x_rotated, Map * map, int *reasonable_angle_ct);
    int monte_carlo (double alpha,
		 double **x, int * x_type, int NX,
		 double **y, int * y_type, int NY,
		 double  *q_best, double *F_best_ptr);
    int opt_quat (double ** x, int NX, int *set_of_directions_x,
		  double ** y, int NY, int *set_of_directions_y,
		  int set_size, double * q, double * rmsd);
    int qmap (double *x0, double *x1, double *y0, double *y1, double * quat);
    int store_sorted (Map * map, int NX, int NY, int *map_best, int map_max,
		      double * z_best, int best_ctr,
		      double z_scr, int  my_map_ctr, int *stored);
    
	    
    
    map_best[0] = -1; /* it is the end-of-array flag */
    if ( *map_ctr >= map_max ) {
	fprintf (stderr, "Map array undersized.\n");
	exit (1);
    }

    smaller = (NX <= NY) ? NX : NY;
 
    /***********************/
    /* memory allocation   */
    /***********************/
    if ( ! (R=dmatrix(3,3) ) ) return 1; /* compiler is bugging me otherwise */
    if ( ! (x_rotated    = dmatrix (NX,3)) ) return 1;
    if ( ! (tr_x_rotated = dmatrix (NX,3)) ) return 1;
    if ( ! (best_quat    = dmatrix (no_top_rmsd,4)) ) return 1;
    if ( ! (best_triple_x    = intmatrix (no_top_rmsd,3)) ) return 1;
    if ( ! (best_triple_y    = intmatrix (no_top_rmsd,3)) ) return 1;
    if ( ! (z_best = emalloc(NX*NY*sizeof(double) )) ) return 1;
    if ( ! (x_type_fudg = emalloc(NX*sizeof(int) )) ) return 1;
    if ( ! (y_type_fudg = emalloc(NY*sizeof(int) )) ) return 1;
    if ( ! (anchor_x = emalloc(NX*sizeof(int) )) ) return 1;
    if ( ! (anchor_y = emalloc(NY*sizeof(int) )) ) return 1;

    penalty_params.custom_gap_penalty_x = NULL;
    penalty_params.custom_gap_penalty_y = NULL;
    //if ( ! (penalty_params.custom_gap_penalty_x = emalloc(NX*sizeof(double) )) ) return 1; 
    //if ( ! (penalty_params.custom_gap_penalty_y = emalloc(NY*sizeof(double) )) ) return 1; 
    /***********************/
    
    /***********************/
    /* expected quantities */
    /***********************/
    avg = avg_sq = stdev = 0.0;
    //if (options.postprocess) {
    if (0) {
	if (F_moments (x, x_type, NX, y, y_type, NY, alpha, &avg, &avg_sq, &stdev)) return 1;
    }
    /***********************/
    

    /***********************/
    /* initialization      */
    /***********************/
    best_ctr   = 0;
    penalty_params.gap_opening   = options.gap_open;
    penalty_params.gap_extension = options.gap_extend;
    penalty_params.endgap        = options.endgap;
    penalty_params.endgap_special_treatment = options.use_endgap;
    /***********************/

    /***************************************/
    /* find reasonble triples of SSEs      */
    /* that correspond in type             */
    /*  and can be mapped onto each other  */
    /***************************************/
    for (top_ctr=0; top_ctr<no_top_rmsd; top_ctr++) {
	best_rmsd[top_ctr] = BAD_RMSD+1;
	best_triple_x[top_ctr][0] = -1;
    }
	
    for (x_ctr=0; x_ctr < NX-2 && !done; x_ctr++) {
	
	for (y_ctr=0; y_ctr < NY-2 && !done; y_ctr++) {

	    if ( y_type[y_ctr] != x_type[x_ctr] ) continue;
	    
	    x_triple[0] = x_ctr;
	    y_triple[0] = y_ctr;
	    
	    if (find_next_triple (x, y,  x_type, y_type,
				  NX, NY,  x_triple, y_triple) ){
		continue;
	    }

	    if ( x_triple[1] < 0 ||  x_triple[2] < 0 ) continue;
	    if ( y_triple[1] < 0 ||  y_triple[2] < 0 ) continue;

	    /* do these three have  kind-of similar layout in space?*/
	    /* is handedness the same? */
	    if ( ! same_hand_triple ( X_rep, x_triple, Y_rep, y_triple, 3)) continue;
	    
	    /* are distances comparab;e? */
	    if (distance_of_nearest_approach ( X_rep, x_triple,
					       Y_rep, y_triple, 3, &rmsd)) continue;
	    if ( rmsd > cutoff_rmsd) continue;
	    
	    /* find q_init that maps the two triples as well as possible*/
	    if ( opt_quat ( x,  NX, x_triple, y, NY, y_triple, 3, q_init, &rmsd)) continue;


	    for (top_ctr=0; top_ctr<no_top_rmsd; top_ctr++) {
		
		if (  rmsd <= best_rmsd[top_ctr] ) {
			    
		    chunk = no_top_rmsd - top_ctr -1;

		    if (chunk) {
			memmove (best_rmsd+top_ctr+1, best_rmsd+top_ctr, chunk*sizeof(double)); 
			memmove (best_quat[top_ctr+1],
				 best_quat[top_ctr], chunk*4*sizeof(double)); 
			memmove (best_triple_x[top_ctr+1],
				 best_triple_x[top_ctr], chunk*3*sizeof(int)); 
			memmove (best_triple_y[top_ctr+1],
				 best_triple_y[top_ctr], chunk*3*sizeof(int)); 
		    }
		    best_rmsd[top_ctr] = rmsd;
		    memcpy (best_quat[top_ctr], q_init, 4*sizeof(double)); 
		    memcpy (best_triple_x[top_ctr], x_triple, 3*sizeof(int)); 
		    memcpy (best_triple_y[top_ctr], y_triple, 3*sizeof(int)); 
		    break;
		}
	    }
	    
	}
    }

# if 0
    for (top_ctr=0; top_ctr<no_top_rmsd; top_ctr++) {
	if ( best_rmsd[top_ctr] > BAD_RMSD ) break;
	printf (" %3d %8.3lf   ", top_ctr,  best_rmsd[top_ctr]);
	vec_out ( best_quat[top_ctr], 4, "quat: ");
	for (t=0; t<3; t++ ) {
	    printf ("\t %3d  %3d \n", best_triple_x[top_ctr][t]+1, best_triple_y[top_ctr][t]+1 );
	}
    }
    exit (1);
# endif

    /*********************************************/
    /*   main loop                               */
    /*********************************************/
    for (top_ctr=0; top_ctr<no_top_rmsd; top_ctr++) {
	
	if ( best_rmsd[top_ctr] > BAD_RMSD ) break;

	quat_to_R (best_quat[top_ctr], R);
	rotate (x_rotated, NX, R, x);


	F_current = F( y, y_type, NY, x_rotated, x_type, NX, alpha);

# if 0	
	printf ("\n***********************************\n");
	printf (" %3d %8.3lf  %8.3lf   ", top_ctr,  best_rmsd[top_ctr], F_current);
	vec_out ( best_quat[top_ctr], 4, "quat: ");
	for (t=0; t<3; t++ ) {
	    printf ("\t %3d  %3d \n", best_triple_x[top_ctr][t]+1, best_triple_y[top_ctr][t]+1 );
	}
# endif
	/* find map which uses the 2 triples as anchors */
	no_anchors = 3;
	find_map (&penalty_params, X_rep, Y_rep, R, alpha, &F_effective, map + (*map_ctr),
		   best_triple_x[top_ctr], best_triple_y[top_ctr], no_anchors);

	x2y = ( map + (*map_ctr) ) ->x2y;
	map_unstable  = 0;
	for (t=0; t<3; t++ ) {
	    if ( x2y[best_triple_x[top_ctr][t]] != best_triple_y[top_ctr][t] ) {
		map_unstable = 1;
	    }
	}
	if ( map_unstable) continue;
	
	/* dna here is not DNA but "distance of nearest approach" */
	cull_by_dna ( X_rep, best_triple_x[top_ctr], 
		      Y_rep, best_triple_y[top_ctr],  3,  map + (*map_ctr), cutoff_rmsd );
	
	//printf ("map after culling by dna:\n");
	//print_map (stdout, map+ (*map_ctr), NULL, NULL,  NULL, NULL, 1);

	/* monte that optimizes the aligned vectors only */
	for (i=0; i<NX; i++) {
	     x_type_fudg[i] = JACKFRUIT;
	}
	for (j=0; j<NY; j++) {
	     y_type_fudg[j] = JACKFRUIT*2;
	}

	no_anchors = 0;

	for (i=0; i<NX; i++) {
	     j = (map+(*map_ctr))->x2y[i];
	     if (j < 0 ) continue;
	     x_type_fudg[i] = x_type[i];
	     y_type_fudg[j] = y_type[j];
	     anchor_x[no_anchors] = i;
	     anchor_y[no_anchors] = j;
	     no_anchors ++;
	}


	if ( opt_quat ( x,  NX, anchor_x, y, NY, anchor_y, no_anchors, q, &rmsd)) continue;
	
	retval = monte_carlo ( alpha,  x, x_type_fudg, NX,
			       y,  y_type_fudg, NY, q, &F_current);
	if (retval) return retval;

	
	if (options.postprocess) {
	    z_scr = stdev ? (F_current - avg)/stdev : 0.0;
	} else {
	    z_scr =  0.0;
	}
	quat_to_R (q, R);
	/* store_image() is waste of time, but perhaps not critical */
	store_image (X_rep, Y_rep, R,  alpha, map + (*map_ctr));
	map_assigned_score ( X_rep, map + (*map_ctr));

	//printf ("map  %2d  assigned score:   %8.3lf      z_score: %8.3lf \n\n",
	//	*map_ctr+1, (map + (*map_ctr)) -> assigned_score, z_scr);


        /*   store the map that passed all the filters down to here*/
	my_map_ctr = *map_ctr;


	map[my_map_ctr].F       = F_current;
	map[my_map_ctr].avg     = avg;
	map[my_map_ctr].avg_sq  = avg_sq;
	map[my_map_ctr].z_score = z_scr;
	memcpy ( map[my_map_ctr].q, q, 4*sizeof(double) );
		
	/* recalculate the assigned score*/

	
	//if (top_ctr==24) exit (1);

	/************************/
	/* store sorted         */
	/************************/
	/* find the place for the new z-score */
	store_sorted (map, NX, NY, map_best, map_max,
		      z_best, best_ctr, -map[my_map_ctr].assigned_score, my_map_ctr, &stored_new);

	if ( stored_new ) { /* we want to keep this map */
	    (*map_ctr) ++;
	    best_ctr++;
	} /* otherwise this map space is reusable */
	    

	/* is this pretty much as good as it can get ? */
	if ( fabs (map[my_map_ctr].assigned_score - smaller)
	     < options.tol )  done = 1;


    }
    map_best[best_ctr] = -1;
  
    
    /******************************************************/
    /* look for the sub-map of a couple of best hits      */
    /******************************************************/
    /* initialization:*/
    
    map_consistence ( NX, NY, NULL, NULL, NULL, NULL, NULL); 
    
    best_ctr = 0;
    while (  map_best[best_ctr] >  -1 ) {
	best_ctr ++;
    }

    //exit (1);
    
    if (best_ctr) {
	int nr_maps = (best_ctr<options.number_maps_cpl)?
	    best_ctr : options.number_maps_cpl;
	int best_i;
	int consistent;
	double z;
	double total_assigned_score, score, best_score = -100;
	double gap_score;

	for (i=0; i<nr_maps; i++) { /* look for the complement */
	    best_i =  map_best[i];
	    
	    /*intialize the (list of) submatch map(s) */
	    if ( !map[best_i].submatch_best) {
		/* for now look for a single map only */
		/* TODO - would it be worth any to look at more maps?*/ 
		int map_max = 1;
		map[best_i].submatch_best = emalloc (map_max*sizeof(int) );
		if (! map[best_i].submatch_best) return 1;
	    }
	    map[best_i].submatch_best[0]    = -1;
	    map[best_i].score_with_children =  0;
	    map[best_i].compl_z_score       =  0;
	    
	    for (j=0; j<best_ctr; j++) {

		if (i==j) continue;
		
		map_complementarity ( map+best_i, map + map_best[j], &z);
			
		map_consistence ( NX, NY, map+best_i, map + map_best[j],
				  &total_assigned_score, &gap_score, NULL);
		consistent = ( (map+best_i)->assigned_score < total_assigned_score
			       && (map + map_best[j])->assigned_score
			       < total_assigned_score);
		if ( consistent ) {
		    score = total_assigned_score;
		    if (  score > best_score ) {
			best_score = score;
			map[best_i].submatch_best[0] = map_best[j];
			map[best_i].score_with_children = total_assigned_score;
			map[best_i].compl_z_score = z;
		    }
		}
	    }

	}
    }
    
    /**********************/
    /* garbage collection */
    gradient_descent (1, 0.0,  NULL, NULL, 0,
			       NULL, NULL, 0, NULL, NULL);
    free_dmatrix (R);
    free_dmatrix (x_rotated);
    free_dmatrix (tr_x_rotated);
    free_dmatrix (best_quat);
    free_imatrix (best_triple_x);
    free_imatrix (best_triple_y);
    free (z_best);
    free (x_type_fudg);
    free (y_type_fudg);
    free (anchor_x);
    free (anchor_y);
     
    if (penalty_params.custom_gap_penalty_x) free (penalty_params.custom_gap_penalty_x);
    if (penalty_params.custom_gap_penalty_y) free (penalty_params.custom_gap_penalty_y);
    /*********************/
    
    return 0;
}

/***************************************/
/***************************************/

int store_sorted (Map * map, int NX, int NY, int *map_best,
		  int map_max,  double * z_best, int best_ctr,
		  double z_scr, int  my_map_ctr, int *stored) {

    int ctr, chunk;
    int correlated;
    double  z;
    
    *stored = 1; /* unless correlated and smaller in z,
		    we will store the map */
    
    /* is this map correlated with something
       that we have already ?*/
    correlated = 0;
    for  (ctr = 0; ctr<best_ctr && ! correlated && ctr < map_max; ctr ++ ) {
	map_complementarity ( map+my_map_ctr, map+ctr,  &z);
	correlated = (z < options.z_max_corr);
	/* TODO store calculated correlations, so I
	   don't have to redo it */
    }
    if ( ctr == map_max ) *stored=0;
    /* correlated: if this score is bigger, replace the old map */
    if ( correlated ) {
	if ( z_best[ctr] <  z_scr ) {
	      /* move */
	    chunk = best_ctr - ctr;
	    if ( chunk ) {
		memmove (z_best+ctr, z_best+ctr+1, chunk*sizeof(double)); 
		memmove (map_best+ctr, map_best+ctr+1, chunk*sizeof(int)); 
	    }
	    best_ctr --;
	} else {
	    /* we'll discard this map entirely */
	    *stored = 0;
	}
    }

    if ( *stored ) { /* if this map is to be stored */
	/* find the place in the list*/
	ctr = 0;
	while (ctr<best_ctr && ( z_scr >= z_best[ctr]) && ctr<map_max) ctr++;

	if ( ctr==map_max ) {/* we alredy have enough maps
			       which are much better */
	    *stored = 0;
	} else {
	    /* move */
	    chunk = best_ctr - ctr;
	    if ( chunk ) {
		memmove (z_best+ctr+1, z_best+ctr, chunk*sizeof(double)); 
		memmove (map_best+ctr+1, map_best+ctr, chunk*sizeof(int)); 
	    }
	    /* store z_score*/
	    z_best[ctr]   = z_scr;
	    map_best[ctr] = my_map_ctr;
	}
    }   

    return 0;
}
/****************************************************************/
/****************************************************************/

int find_next_triple (double **X, double **Y, 
		    int *x_type, int *y_type,
		    int NX, int NY, int *x_triple, int *y_triple){

    int i, j, t, done = 0;
    int x_ctr;
    int y_ctr;
    
    /* next Y[y_ctr] cannot be numerically negative of X[x_ctr] */
    x_ctr = x_triple[0];
    y_ctr = y_triple[0];
    if ( x_type[x_ctr] != y_type[y_ctr] ) return 1;


    for (t=1; t<=2; t++ ) {
	x_ctr = x_triple[t-1];
	y_ctr = y_triple[t-1];

	if ( x_ctr == -1 || y_ctr==-1) break;
	
	x_triple[t] = -1;
	y_triple[t] = -1;
	done = 0;
	for (i=x_ctr+1; i<NX && !done; i++ ) {
	    for (j=y_ctr+1; j<NY && !done; j++ ) {
	
		if ( x_type[i] == y_type[j] ) {
		    x_triple[t] = i;
		    y_triple[t] = j;
		    done = 1;
		    
		}
	    }
	}
    }


    

    return 0;
}


/**********************************************************/
/**********************************************************/
/* check that the two are not mirror images - count on it being  a triple*/
int same_hand_triple (Representation * X_rep,  int *set_of_directions_x,
	       Representation * Y_rep, int *set_of_directions_y, int set_size) {
    
    double **x    = X_rep->full;
    double **x_cm = X_rep->cm;
    double **y    = Y_rep->full;
    double **y_cm = Y_rep->cm;
    double cm_vector[3], avg_cm[3], cross[3],  dx, dy, aux;
    int i, ray_a, ray_b,  ray_c, a, b, c;
    
    if (set_size !=3 ) return 0;

 
    a = 0;
    b = 1;
    c = 2;
    
    /*****************************/
    /*****************************/
    ray_a = set_of_directions_x[a];
    ray_b = set_of_directions_x[b];
    ray_c = set_of_directions_x[c];
    /* I better not use the cross prod here: a small diff
       int he angle makeing them pointing toward each other or
       away from each other changes the direction of cross prod;
       rather, use one vector, and distance between the cm's as the other */
    for (i=0; i<3; i++ ) {
	cm_vector[i] = x_cm[ray_b][i] - x_cm[ray_a][i];
    }
    normalized_cross (x[ray_a], x[ray_b], cross, &aux); 
    /* note I am makning another cm vector here */
    for (i=0; i<3; i++ ) {
	avg_cm[i] = (x_cm[ray_b][i] + x_cm[ray_a][i])/2;
	cm_vector[i] = x_cm[ray_c][i] - avg_cm[i];
    }
    unnorm_dot (cm_vector, cross, &dx);

    /*****************************/
    /*****************************/
    ray_a = set_of_directions_y[a];
    ray_b = set_of_directions_y[b];
    ray_c = set_of_directions_y[c];
    for (i=0; i<3; i++ ) {
	cm_vector[i] = y_cm[ray_b][i] - y_cm[ray_a][i];
    }
    normalized_cross (y[ray_a], y[ray_b], cross, &aux); 
    /* note I am makning another cm vector here */
    for (i=0; i<3; i++ ) {
	avg_cm[i] = (y_cm[ray_b][i] + y_cm[ray_a][i])/2;
	cm_vector[i] = y_cm[ray_c][i] - avg_cm[i];
    }
    /*note: unnorm_dot thinks it is getting unit vectors,
      and evrything that is >1 will be "rounded" to 1
      (similarly for -1) - it doesn't do the normalization itself*/
    unnorm_dot (cm_vector, cross, &dy);

     
    if ( dx*dy < 0 ) return 0;

    
    return 1;   /* this isn't err value - the handedness is the same */
}

/**********************************************************/
int distance_of_nearest_approach ( Representation * X_rep,  int *set_of_directions_x,
	       Representation * Y_rep, int *set_of_directions_y,
	       int set_size,  double * rmsd_ptr) {
    
    double **x    = X_rep->full;
    double **x_cm = X_rep->cm;
    double **y    = Y_rep->full;
    double **y_cm = Y_rep->cm;
    double cm_vector[3], cross[3],  distance_x, distance_y;
    double aux, rmsd;
    int i, ray_a, ray_b, a, b, prev_next, norm;
    
    if (set_size <=1 ) return 1;

    rmsd = 0.0;
    norm = 0;
    /* the rmsd for the remaining vectors is ... */
    for (a=0; a<set_size; a++) {

	for ( prev_next = -1;  prev_next <= +1; prev_next +=2 ) {
	    b = (set_size+a+prev_next)%set_size;
	
	    ray_a = set_of_directions_x[a];
	    ray_b = set_of_directions_x[b];
	    /* distance of nearest approach of ray b
	       to the cm of a, in the set of directions x */
	    for (i=0; i<3; i++ ) {
		cm_vector[i] = x_cm[ray_b][i] - x_cm[ray_a][i];
	    }
	    normalized_cross (cm_vector, x[ray_a], cross, &distance_x);

	    ray_a = set_of_directions_y[a];
	    ray_b = set_of_directions_y[b];
	    /* distance of nearest approach of ray b
	       to the cm of a, in the set of directions y */
	    for (i=0; i<3; i++ ) {
		cm_vector[i] = y_cm[ray_b][i] - y_cm[ray_a][i];
	    }
	    normalized_cross (cm_vector, y[ray_a], cross, &distance_y);


	    aux   = distance_x-distance_y;
	    rmsd += aux*aux;
	    norm ++;
# if 0
	    printf ("%d, %d   x:  %2d  %2d  y:  %2d  %2d  \n", a, b,
		    set_of_directions_x[a], set_of_directions_x[b], 
		    set_of_directions_y[a], set_of_directions_y[b]); 
	    printf (" distance x:  %8.3lf  distance y:  %8.3lf   difference:   %8.3lf \n",
		    distance_x,  distance_y, fabs (distance_x-distance_y));
# endif
	}
	
    }

    rmsd /= norm;
    rmsd  = sqrt(rmsd);
    
    //printf ("***** rmsd:  %8.3lf  \n\n", rmsd);

    *rmsd_ptr = rmsd;
    
    return 0;
}

/**********************************************************/
/**********************************************************/
int cull_by_dna (Representation * X_rep, int *set_of_directions_x,
		 Representation * Y_rep, int *set_of_directions_y,
		 int set_size, Map *map, double cutoff_rmsd) {
    
    double **x    = X_rep->full;
    double **x_cm = X_rep->cm;
    double **y    = Y_rep->full;
    double **y_cm = Y_rep->cm;
    double cm_vector[3], cross[3],  distance_x, distance_y;
    double aux, rmsd;
    int NX=X_rep->N_full;
    int i, j, ray_a, ray_b, a;
    int in_triple, coord_ctr, norm;
    
    if (set_size < 1 ) return 1;


    for (i=0; i<NX; i++) {
	j = map->x2y[i];
	if ( j<0) continue;

	/* check that i is not in the triple here */
	in_triple = 0;
	for (a=0; a<set_size; a++) {
	    if ( i==set_of_directions_x[a]) {
		in_triple =1;
		break;
	    }
	}
	if (in_triple) continue;

	rmsd = 0.0;
	norm = 0;
	/* the rmsd for the dna of this element to the triple */
	for (a=0; a<set_size; a++) {

	    /**************************************/
	    ray_a = set_of_directions_x[a];
	    ray_b = i;
	    /* distance of nearest approach of ray b
	       to the cm of a, in the set of directions x */
	    for (coord_ctr=0; coord_ctr<3; coord_ctr++ ) {
		cm_vector[coord_ctr] = x_cm[ray_b][coord_ctr] - x_cm[ray_a][coord_ctr];
	    }
	    normalized_cross (cm_vector, x[ray_a], cross, &distance_x);

	    ray_a = set_of_directions_y[a];
	    ray_b = j;
	    /* distance of nearest approach of ray b
	       to the cm of a, in the set of directions y */
	    for (coord_ctr=0; coord_ctr<3; coord_ctr++ ) {
		cm_vector[coord_ctr] = y_cm[ray_b][coord_ctr] - y_cm[ray_a][coord_ctr];
	    }
	    normalized_cross (cm_vector, y[ray_a], cross, &distance_y);

	    aux   = distance_x-distance_y;
	    rmsd += aux*aux;
	    norm ++;

	    /**************************************/
	    ray_b = set_of_directions_x[a];
	    ray_a = i;
	    /* distance of nearest approach of ray b
	       to the cm of a, in the set of directions x */
	    for (coord_ctr=0; coord_ctr<3; coord_ctr++ ) {
		cm_vector[coord_ctr] = x_cm[ray_b][coord_ctr] - x_cm[ray_a][coord_ctr];
	    }
	    normalized_cross (cm_vector, x[ray_a], cross, &distance_x);

	    ray_b = set_of_directions_y[a];
	    ray_a = j;
	    /* distance of nearest approach of ray b
	       to the cm of a, in the set of directions y */
	    for (coord_ctr=0; coord_ctr<3; coord_ctr++ ) {
		cm_vector[coord_ctr] = y_cm[ray_b][coord_ctr] - y_cm[ray_a][coord_ctr];
	    }
	    normalized_cross (cm_vector, y[ray_a], cross, &distance_y);

	    aux   = distance_x-distance_y;
	    rmsd += aux*aux;
	    norm ++;

	}

	rmsd /= norm;
	rmsd = sqrt(rmsd);

# if 0
	printf ("%2d  %2d : rmsd %8.3lf\n", i+1, j+1, rmsd);
# endif
	

	/* if rmsd is bigger than the cutoff, lose this from the mapping */
	if ( rmsd > cutoff_rmsd) {
	    map->x2y[i] = -1;
	    map->y2x[j] = -1;
	}
    }
    
     
    return 0;
}

/**********************************************************/
/**********************************************************/
int opt_quat ( double ** x, int NX, int *set_of_directions_x,
	       double ** y, int NY, int *set_of_directions_y,
	       int set_size, double * q, double * rmsd) {

    
    double * x_sub[set_size], * y_sub[set_size];
    int  ctr;
    int  i, j;
 
    double ATA     [4][4] = {{0.0}};
    double prev_ATA[4][4] = {{0.0}};
    double ATA_sum [4][4] = {{0.0}};
    double a[3] = {0.0}, b[3] = {0.0};
    
    int add_matrices  (double matrix1[4][4],double matrix2[4][4],
		       double result[4][4]);
    int construct_ATA (double ATA[4][4], double a[3], double  b[3]);

    /* note how we pass the matrix: pointer to the first element in the block */
    void dsyev_ (char * jobz, char *uplo,  int *n,
		  double *A, int * lda, double * w, double * work, int * lwork, int *info);

    if (!set_size) {
	*rmsd = -1;
	return 1;
    }

    memset ( &(q[0]), 0, 4*sizeof(double) );

    
    /* find the subset */
    ctr = 0;
    for ( ctr =0; ctr < set_size; ctr++ ) {
	x_sub[ctr] =  x[set_of_directions_x[ctr]];
	y_sub[ctr] =  y[set_of_directions_y[ctr]];
    }

    /* check: */
    if (0) {
	printf (" Number of vectors to match: %d. \n", set_size);
	for ( ctr =0; ctr < set_size; ctr++ ) {
	    printf ("\t x%1d   %10.4lf  %10.4lf  %10.4lf   ",
		    ctr, x_sub[ctr][0], x_sub[ctr][1], x_sub[ctr][2]);
	    printf ("\t y%1d   %10.4lf  %10.4lf  %10.4lf \n",
		    ctr, y_sub[ctr][0], y_sub[ctr][1], y_sub[ctr][2]);
	}
	exit (1);
    }

     
    /* B = ATA_sum matrix to diagonalize in order to get the quaternion */
    for ( ctr =0; ctr < set_size; ctr++ ) {
   	for (i=0; i<3; i++ ) {
	    a[i] = y_sub[ctr][i] + x_sub[ctr][i];
	    b[i] = y_sub[ctr][i] - x_sub[ctr][i];
	}
 	construct_ATA (ATA, a, b);
	add_matrices (prev_ATA, ATA, ATA_sum);
	memcpy (prev_ATA[0], ATA_sum[0], 4*4*sizeof(double));
    }
    for (i=0; i<4; i++ ) {
	for (j=0; j<4; j++ ) {
	    ATA_sum[i][j] /= set_size;
	}
    }
    /* diagonalize ATA_sum - the eigenvector corresponsing to the
       smallest lambda is the quaternion we are looking for; the
       eigenvalue is the rmsd*/
    /* use the nomenclature from dsyev*/
    char jobz= 'V'; /*Compute eigenvalues and eigenvectors.*/
    char uplo= 'U'; /* Upper triangle of A (the matrix we are diagonalizing) is stored; */
    int  n = 4;     /* order and the leading dimension of A */
    int  lda = 4;
    double ** A;
    int  info;
    int  lwork = 200;
    double w [4];
    double work[200];
    
    if ( !( A=dmatrix(4,4) ) ) exit (1);
    memcpy (A[0], ATA_sum[0], 4*4*sizeof(double));


   /* note how we pass the matrix: */
    dsyev_ ( &jobz, &uplo,  &n, A[0], &lda, w, work, &lwork, &info);
    if (  ! info) {
	*rmsd = sqrt (w[0]);
	for (i=0; i<4; i++ ) q[i] = A[0][i];
	if (0) {
	    /* w contains the eigenvalues */
	    printf ("\n");
	    for (i=0; i<4; i++ ) printf ("%8.3lf ", w[i]);
	    printf ("\nrmsd: %8.3lf \n", *rmsd);
	    printf ("quat:\n");
	    for (i=0; i<4; i++ ) printf ("%8.3lf ", q[i]);
	    printf ("\n");
	    /* printf (" opt lwork: %d\n", (int) work[0]); */
	}
    } else {
	fprintf (stderr, "Error in dsyev().\n");
	exit (1);
    }
    
   
    
    free_dmatrix(A);

    return 0;
    
}

/************************************/
/************************************/
