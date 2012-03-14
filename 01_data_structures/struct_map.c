# include "struct.h"


/********************************/
int initialize_map (Map *map, int NX, int NY ) {

    memset (map, 0, sizeof(Map));
    map->x2y = emalloc (NX*sizeof(int));
    map->y2x = emalloc (NY*sizeof(int));
    map->x2y_size = NX;
    map->y2x_size = NY;
    map->cosine  = dmatrix (NX, NY);
    map->image   = dmatrix (NX, NY);
    if ( !map->x2y ||  !map->y2x ||
	 !map->cosine || !map->image ) {
	return 1;
    }
    map->submatch_best = NULL;

    return 0;
}
/********************************/
int free_map (Map *map) {

    if ( !map) return 0;
    
    free (map->x2y);
    free (map->y2x);
    map->x2y_size = 0;
    map->y2x_size = 0;
    free_dmatrix (map->cosine);
    free_dmatrix (map->image); 
    if (map->submatch_best) {
	free (map->submatch_best );
    }


    return 0;
}

/********************************/
/********************************/
/********************************/

int store_image ( Representation *X_rep,  Representation *Y_rep, 
		  double **R, double alpha,  Map *map ){
    
    int i,j, i_is_helix;
    int NX=X_rep->N_full, NY=Y_rep->N_full;
    double cosine;
    double alpha_sq = alpha*alpha;
    double ** x_full_rotated = NULL, **y_full = Y_rep->full;
    double length_dif, exponent;

    if ( ! (x_full_rotated=dmatrix (NX,3) ) ) return 1;
    
    rotate (x_full_rotated, X_rep->N_full, R, X_rep->full);
    
    for (i=0; i < NX; i++ ) {
	i_is_helix =  ( X_rep->full_type[i] == HELIX);
	for (j=0; j < NY; j++) {

	    if ( X_rep->full_type[i] == Y_rep->full_type[j] ) {
		unnorm_dot (x_full_rotated[i], y_full[j], &cosine);
		exponent = 2*(1.0-cosine)/alpha_sq;
		if ( exponent < 10 ) {
		    map->image[i][j]  = exp (-exponent);
		    /* are we using the info about the length? */
		    if (options.use_length && cosine > 0.0) {
			length_dif = abs (X_rep->length[i] - Y_rep->length[j]);
			if ( i_is_helix ) {
			    map->image[i][j]  *=  1/(1.0+length_dif/options.H_length_mismatch_tol);
			} else {
			    map->image[i][j]  *=  1/(1.0+length_dif/options.S_length_mismatch_tol);
			}
		    }
		} else {
		    map->image[i][j] = 0.0;
		}
		map->cosine[i][j] = cosine;
	    } else {
		map->image[i][j]  = options.far_far_away;
		map->cosine[i][j] = -0.9999;
	    }
	}
    }
    free_dmatrix (x_full_rotated);
    
    return 0;
}


/***************************************/
/***************************************/
/***************************************/
int construct_translation_vecs ( Representation *X_rep,  Representation *Y_rep,
				 Map *map ){

    int i, j, k, map_size;
    int NX=X_rep->N_full, NY=Y_rep->N_full;
    double norm_x, norm_y, aux;
 
    /********************************************/
    /* find origin (geom center of mapped SSEs) */
    for (k=0; k<3; k++) {
	X_rep->origin[k] = 0;
	Y_rep->origin[k] = 0;
    }

    map_size = 0;
    for (i=0; i<NX; i++ ) {
	
	j =  map->x2y[i];
	if ( j < 0 ) continue;

	for (k=0; k<3; k++) {
	    X_rep->origin[k] += X_rep->cm[i][k];
	    Y_rep->origin[k] += Y_rep->cm[j][k];
	}
	map_size++;
    }

    for (k=0; k<3; k++) {
	X_rep->origin[k] /= map_size;
	Y_rep->origin[k] /= map_size;
    }

    
    /*******************************************/
    /* find translation vectors                */
    memset (X_rep->translation[0], 0, 3*NX*sizeof(double));
    memset (Y_rep->translation[0], 0, 3*NY*sizeof(double));
    
    for (i=0; i<NX; i++ ) {
	j =  map->x2y[i];
	if ( j < 0 ) continue;
	
	norm_x = norm_y = 0.0;
	for (k=0; k<3; k++) {
	    X_rep->translation[i][k] = aux = X_rep->cm[i][k] -  X_rep->origin[k];
	    norm_x += aux*aux;
	    Y_rep->translation[j][k] = aux = Y_rep->cm[j][k] -  Y_rep->origin[k];
	    norm_y += aux*aux;
	}
	norm_x = sqrt (norm_x);
	norm_y = sqrt (norm_y);

	X_rep->transl_norm[i] = norm_x;
	Y_rep->transl_norm[j] = norm_y;
    }
    
     
  
   
    
    return 0;
    
}
/***************************************/
/***************************************/
/***************************************/

int map_quality_metrics (Representation *X_rep, Representation *Y_rep,
			 double ** tr_x_rotated, Map * map, int *reasonable_match_ct) {
    
    int i,j,k;
    int NX = X_rep->N_full;
    int map_size, old_map_size;
    int done, loop_ct;
    double rmsd, dist_sq, aux;
    double avg_length_mismatch;

    /***********************************************/
    /* number of reasonable translation directions */
    map_size = 0;
    for (i=0; i < NX; i++ ) {	
	j = map->x2y[i];
	if ( j < 0) continue;
	map_size++;
    }



    /***********************************************/
    /* too distant .....                           */
    done = 0;
    loop_ct = 0;
    while ( !done ) {
	
	old_map_size = map_size;
	construct_translation_vecs (X_rep, Y_rep, map);
	rmsd     = 0;
	map_size = 0;
	
	for (i=0; i<NX; i++ ) {
	    j =  map->x2y[i];
	    if ( j < 0 ) continue;


	    /* note since not normalized this won't really be the cosine*/
	    /* I couldn't get this to be really useful */
	    /* unnorm_dot ( tr_x_rotated[i], Y_rep->translation[j], &cosine); */
	    /* now, this cosine refers to translation vectors */
	    /* if (cosine < 0.0 ) { */
/* 		map->x2y[i] = FAR_FAR_AWAY; */
/* 		map->y2x[j] = FAR_FAR_AWAY; */
/* 		continue; */
/* 	    } */
 	    
	
	    dist_sq = 0.0;
	    for (k=0; k<3; k++) {
		aux      = tr_x_rotated[i][k]-Y_rep->translation[j][k];
		dist_sq += aux*aux;
	    }
	    printf ("\t\t %3d --> %3d  dist:  %8.3lf\n", 
		    i+1, j+1,  sqrt(dist_sq) );

	    if (dist_sq < 100) { /**** hard coded parameter <<<<<<<<<<<< */
		rmsd += dist_sq;
		map_size ++;
	    } else {
		map->x2y[i] = FAR_FAR_AWAY;
		map->y2x[j] = FAR_FAR_AWAY;
	    }
	}
	loop_ct ++;
	done = (map_size==old_map_size || loop_ct >2);
	    
    }
    
    if (map_size) {
	rmsd /= map_size;
	map->rmsd = sqrt (rmsd);
    } else {
	map->rmsd = 100.0;
    }
    
    *reasonable_match_ct = map_size;

    
    /***********************************************/
    /* length mismatch of mapped SSEs              */
    avg_length_mismatch = 0;
    for (i=0; i < NX; i++ ) {
	j =  map->x2y[i];
	if ( j <  0 )  continue;
	avg_length_mismatch += fabs(X_rep->length[i] - Y_rep->length[j]);
    }
    map->size    =  map_size;
    avg_length_mismatch /= map_size;
    map->avg_length_mismatch  = avg_length_mismatch;

    
    
    return 0;
}
/***************************************/
/***************************************/
/***************************************/
 int find_map ( Penalty_parametrization * penalty_params,
	       Representation *X_rep,  Representation *Y_rep,
	       double ** R,
	       double alpha, double * F_effective,  Map *map,
	       int *anchor_x, int *anchor_y, int anchor_size) {
    
    int i, j;
    int NX=X_rep->N_full, NY=Y_rep->N_full;
    double F_eff, aln_score;
    
    /* "image" is the matrix of values (exp weight)*(-1 for SSE mismatch) */
    store_image (X_rep, Y_rep, R,  alpha, map);

    /* if the anchors are given, somebody insists
       they should be considered as already aligned */

    if (penalty_params->custom_gap_penalty_x)
	memset (penalty_params->custom_gap_penalty_x, 0, NX*sizeof(double) );
    if (penalty_params->custom_gap_penalty_y)
	memset (penalty_params->custom_gap_penalty_y, 0, NY*sizeof(double) );
    
    if (anchor_x && anchor_y) {
	 int a;
         for (a=0; a<anchor_size; a++){
	     
	      i = anchor_x[a];
	      if (penalty_params->custom_gap_penalty_x)
		  penalty_params->custom_gap_penalty_x[i] = options.far_far_away;
	      for (j=0; j<NY; j++) {
		   if ( j ==  anchor_y[a] )  continue;
		   map->image[i][j]  =  options.far_far_away;
		   map->cosine[i][j] = -0.9999;
	      }
	      
	      j = anchor_y[a];
	      if (penalty_params->custom_gap_penalty_y)
		  penalty_params->custom_gap_penalty_y[j] = options.far_far_away;
	      for (i=0; i<NX; i++) {
		   if ( i ==  anchor_x[a] )  continue;
		   map->image [i][j] =  options.far_far_away;
		   map->cosine[i][j] = -0.9999;
	      }

	 }
    }
    
    /* dynamic programming using the "image" or the hungarian algorithm */
    
    if (options.sw_score) {
        smith_waterman (penalty_params, NX, NY, map->image, map->x2y, map->y2x, &aln_score);
    } else {
        hungarian_alignment (NX, NY, map->image, map->x2y, map->y2x, &aln_score);
    }
    
    
    map_assigned_score (X_rep, map);

    F_eff = 0;
    for (i=0; i<NX; i++) {
	 j = map->x2y[i];
	 if ( j<0) continue;
	 if ( map->cosine[i][j] < options.far_away_cosine) {
	      map->x2y[i] = -1; 
	      map->y2x[j] = -1; /*everything crashes; investigate later */
	 } else {
	      F_eff --;
	 }    
    }

# if 0	
    /* for F effective, go back to the compact rep */
    F_eff = 0;
    for (i_rep=0; i_rep<X_rep->N_compact; i_rep++ ) {
	for (ctr=1; ctr <= X_rep->represents[i_rep][0]; ctr++) {
	    i = X_rep->represents[i_rep][ctr];
	    j = map->x2y[i];
	    
	    if ( j<0) continue;
	    
	    /* get rid of small cosines here (x2y is coming from SW or NW)
	       it can call something mapped even if the cosine is poor */
	    if ( map->cosine[i][j] < options.far_away_cosine) {
		map->x2y[i] = -1; 
		//map->y2x[j] = -1; /*everything crashes; investigate later */
	    } else {
		F_eff --;
		break; /*count only once */
	    }
	}
    }
# endif

    
    *F_effective = F_eff;
    
    return 0;
}

/************************************/
/************************************/
int map_assigned_score ( Representation *X_rep,  Map* map) {

    int i,j;
    int NX = X_rep->N_full;
    
    map->size = 0;
    map->assigned_score = 0;
    for (i=0; i < NX; i++ ) {
	j =  map->x2y[i];
	if ( j <  0 )  continue;
	
	map->assigned_score += map->image[i][j];
    }
    
    /* what's the difference btw this and map->size?*/
    map->matches = match_length (NX, map->x2y);
    return 0;

    
}
/************************************/
/************************************/
int smith_waterman (Penalty_parametrization *params, int max_i, int max_j, double **similarity,
		    int *map_i2j, int * map_j2i, double * aln_score) {

    double **F; /*alignment_scoring table*/
    char ** direction;
    double *custom_gap_penalty_x =  params->custom_gap_penalty_x;
    double *custom_gap_penalty_y =  params->custom_gap_penalty_y;
    double gap_opening   = params->gap_opening;
    double gap_extension = params->gap_extension;
    double endgap        = params->endgap;
    double penalty;
    double i_sim = 0.0, j_sim = 0.0, diag_sim = 0.0, max_sim = 0.0;
    double F_max;

    int use_endgap = params->endgap_special_treatment;
    int F_max_i, F_max_j;
    int i,j;

     /* allocate F */
    if ( ! (F = dmatrix( max_i+1, max_j+1)) ) return 1;
    if ( ! (direction = chmatrix ( max_i+1, max_j+1)) ) return 1;

    
    
    /* fill the table */
    F_max = options.far_far_away;
    F_max_i = 0;
    F_max_j = 0;
   
    for (i=0; i<= max_i; i++) {
	for (j=0; j<=max_j; j++) {

	    if ( !i && !j ) {
		F[0][0] = 0;
		direction[i][j] = 'd';
		continue;
	    }
	    
	    if ( i && j ){
		
		if ( direction[i-1][j] == 'i' ) {
		    /*  gap extension  */
		    if ( custom_gap_penalty_x && custom_gap_penalty_x[i-1] < 0 ) {
			penalty = custom_gap_penalty_x[i-1];
		    } else {
			penalty = (use_endgap&&j==max_j) ? endgap : gap_extension;
		    }
		} else {
		    /*  gap opening  */
		    if ( custom_gap_penalty_x && custom_gap_penalty_x[i-1] < 0 ) {
			penalty = custom_gap_penalty_x[i-1];
		    } else {
			penalty = (use_endgap&&j==max_j) ? endgap : gap_opening;
		    }
		}
		i_sim =  F[i-1][j] + penalty;
		
		if ( direction[i][j-1] == 'j' ) {
		    if ( custom_gap_penalty_y && custom_gap_penalty_y[j-1] < 0 ) {
			penalty = custom_gap_penalty_y[j-1];
		    } else {
			penalty = (use_endgap&&i==max_i) ? endgap : gap_extension;
		    }
		} else {
		    if ( custom_gap_penalty_y && custom_gap_penalty_y[j-1] < 0 ) {
			penalty = custom_gap_penalty_y[j-1];
		    } else {
			penalty = (use_endgap&&i==max_i) ? endgap : gap_opening;
		    }
		}
		j_sim = F[i][j-1] +  penalty;
       	
		
		diag_sim =  F[i-1][j-1] + similarity [i-1][j-1] ;
		

		max_sim = diag_sim;
		direction[i][j] = 'd';
		if ( i_sim > max_sim ){
		    max_sim = i_sim;
		    direction[i][j] = 'i';
		}
		if ( j_sim > max_sim ) {
		    max_sim = j_sim;
		    direction[i][j] = 'j';
		}
		
	    } else if ( j ) {
		
		if ( custom_gap_penalty_y && custom_gap_penalty_y[j-1] < 0 ) {
		    penalty = custom_gap_penalty_y[j-1];
		} else {
		    if ( use_endgap) {
			penalty = endgap;
		    } else {
			if ( direction[i][j-1] =='j' ) {
			    penalty = gap_extension;
			} else {
			    penalty = gap_opening;
			}
		    }
		}
		j_sim = F[i][j-1] + penalty;
		max_sim = j_sim;
		direction[i][j] = 'j';

		
	    } else if ( i ) {
		
		if ( custom_gap_penalty_x && custom_gap_penalty_x[i-1] < 0 ) {
		    penalty = custom_gap_penalty_x[i-1];
		} else {
		    if ( use_endgap) {
			penalty = endgap;
		    } else {
			if ( direction[i-1][j] == 'i' ) {
			    penalty =  gap_extension;
			} else {
			    penalty =  gap_opening;
			}
		    }
		}
		i_sim = F[i-1][j] + penalty;
		max_sim = i_sim;
		direction[i][j] = 'i';
		

	    } 

	    if ( max_sim < 0.0 ) max_sim = 0.0;
	    
	    F[i][j] = max_sim;
	    if ( F_max < max_sim ) {
		/* TODO: tie break here */
		F_max = max_sim;
		F_max_i = i;
		F_max_j = j;
		
	    }
	    
	}
    }
    
    /*retrace from the maximum element*/
    i= max_i;
    for( i= max_i; i>F_max_i; i--) map_i2j[i-1] = options.far_far_away;;
    for( j= max_j; j>F_max_j; j--) map_j2i[j-1] = options.far_far_away;;
    i = F_max_i;
    j = F_max_j;
    *aln_score = F[i][j]; 
    while ( i>0 ||  j >0 ) {
	//printf (" %4d  %4d  %8.3f  \n", i, j, F[i][j]);
	if ( i<0 || j<0 ) {
	    fprintf ( stderr, "Retracing error.\n");
	}
	switch ( direction[i][j] ) {
	case 'd':
	    //printf ( " %4d  %4d \n",  i, j);
	    map_i2j [i-1] = j-1;
	    map_j2i [j-1] = i-1;
	    i--;
	    j--; 
	    break;
	case 'i':
	    //printf ( " %4d  %4d \n",  i, -1);
	    map_i2j [i-1] = options.far_far_away;
	    i--; 
	    break; 
	case 'j':
	    //printf ( " %4d  %4d \n",  -1, j);
	    map_j2i [j-1] = options.far_far_away;
	    j--; 
	    break; 
	default: 
	    fprintf ( stderr, "Retracing error.\n");
		
	} 
    }

    /* free */ 
    free_dmatrix (F);
    free_cmatrix (direction);
    
    return 0; 
   
    
}

/**
 * Function that calculate alignment, and maps structures using Hungarian algorithm
 * @param NX
 * @param NY
 * @param similarity
 * @param x2y
 * @param y2x
 * @param alignment
 * @return 
 */

int hungarian_alignment (int NX, int NY, double **similarity, int * x2y, int * y2x, double * alignment ) {
    int i,j;
    *alignment  = 0;
    for (i =0; i < NX; ++i) x2y[i] = -10;
    for (i =0; i < NY; ++i) y2x[i] = -10;
     
    
    int multiplier = 1000; // precision level for conversion of double to int
    int **scoring_matrix;
    scoring_matrix = intmatrix(NX, NY);
    
    similarity_to_scoring(NX, NY, multiplier, similarity, scoring_matrix );
    
    
    hungarian_problem_t p;
    
    int matrix_size = hungarian_init(&p, scoring_matrix, NX,NY, HUNGARIAN_MODE_MAXIMIZE_UTIL);
    if (matrix_size < 0) printf("wrong matrix size for Hungarian algorithm\n");

    hungarian_solve(&p);

    // checking if the number of rows is greater than number of columns. Note: Hungarian 
    if (NX >= NY) {
        for (i = 0; i < NX; ++i) {
            for (j = 0; j < NY; ++j) {
                if (p.assignment[i][j] && similarity[i][j] > 0) {
                    x2y[i] = j;
                    y2x[j] = i;
                    *alignment += similarity[i][j];
                }
            }
        }
    }
    
  /* free used memory */
    hungarian_free(&p);
    free_imatrix(scoring_matrix);
    
    return 0;
    
}

/**
 * Function that convert similarity matrix of doubles to scoring matrix of integers
 * Scoring matrix is necessary to the Hungarian algorithm 
 * @param m input double matrix
 * @param rows
 * @param cols
 * @return matrix of integers
 */

void similarity_to_scoring(int NX, int NY, int multiplier, double** similarity, int ** hungarian_alignment ) {
  int i,j;
  
  for(i=0;i<NX;i++) {
      for(j=0;j<NY;j++) hungarian_alignment[i][j] = similarity[i][j] * multiplier;
  }
}




/************************************/
/************************************/
int needleman_wunsch (int max_i, int max_j, double **similarity,
		      int *map_i2j, int *map_j2i, double * aln_score) {

    double **F; /*alignment_scoring table*/
    char ** direction;
    double gap_opening   = options.gap_open;
    double gap_extension = options.gap_extend;
    double endgap        = options.endgap;
    double penalty;
    double i_sim = 0.0, j_sim = 0.0, diag_sim = 0.0, max_sim = 0.0;
    int use_endgap = options.use_endgap;
    int i,j;

    /* allocate F */
    if ( ! (F = dmatrix( max_i+1, max_j+1)) ) return 1;
    if ( ! (direction = chmatrix ( max_i+1, max_j+1)) ) return 1;


    /* fill the table */
    for (i=0; i<= max_i; i++) {
	for (j=0; j<=max_j; j++) {

	    if ( !i && !j ) { /* upper left corner */
		F[0][0] = 0;
		direction[i][j] = 'd';
		continue;
	    }
	    
	    if ( i && j ){ 
		if ( direction[i-1][j] == 'i' ) {
		    /*  gap extension  */
		    penalty = (use_endgap&&j==max_j) ? endgap : gap_extension;		    
		} else {
		    /*  gap opening  */
		    penalty = (use_endgap&&j==max_j) ? endgap : gap_opening;
		}
		i_sim =  F[i-1][j] + penalty;

		
		if ( direction[i][j-1] =='j' ) {
		    penalty = (use_endgap&&i==max_i) ? endgap : gap_extension;		    
		} else {
		    penalty = (use_endgap&&i==max_i) ? endgap : gap_opening;		    
		}
		j_sim = F[i][j-1] +  penalty;
       	
		
		diag_sim =  F[i-1][j-1] + similarity [i-1][j-1] ;
		
	    } else if ( j ) {
		
		if ( use_endgap) {
		    penalty = endgap;
		} else {
		    if ( direction[i][j-1] =='j' ) {
			penalty = gap_extension;
		    } else {
			penalty = gap_opening;
		    }
		}
		j_sim = F[i][j-1] + penalty;
		
		i_sim = diag_sim = options.far_far_away;

	    } else if ( i ) {
		if ( use_endgap) {
		    penalty = endgap;
		} else {
		    if ( direction[i-1][j] == 'i' ) {
			penalty =  gap_extension;
		    } else {
		        penalty =  gap_opening;
		    }
		}
		i_sim = F[i-1][j] + penalty;
		
		j_sim = diag_sim = options.far_far_away;
		
	    } 

	    max_sim = diag_sim;
	    direction[i][j] = 'd';
	    if ( i_sim > max_sim ){
		max_sim = i_sim;
		direction[i][j] = 'i';
	    }
	    if ( j_sim > max_sim ) {
		max_sim = j_sim;
		direction[i][j] = 'j';
	    }

	    F[i][j] = max_sim;
	    
	}
    }
    *aln_score = F[max_i][max_j];
    
    /*retrace*/
    i = max_i;
    j = max_j;
    while ( i>0 ||  j >0 ) {
	//printf (" %4d  %4d  %8.3f  \n", i, j, F[i][j]);
	switch ( direction[i][j] ) {
	case 'd':
	    //printf ( " %4d  %4d \n",  i, j);
	    map_i2j [i-1] = j-1;
	    map_j2i [j-1] = i-1;
	    i--;
	    j--; 
	    break;
	case 'i':
	    //printf ( " %4d  %4d \n",  i, -1);
	    map_i2j [i-1] = options.far_far_away;
	    i--; 
	    break; 
	case 'j':
	    //printf ( " %4d  %4d \n",  -1, j);
	    map_j2i [j-1] = options.far_far_away;
	    j--; 
	    break; 
	default: 
	    fprintf ( stderr, "Retracing error.\n");
		
	} 
    }

    /* free */ 
    free_dmatrix (F);
    free_cmatrix (direction);
    
    return 0; 
   
    
}

/************************************************************/
/************************************************************/

int match_length (int N, int * x2y) {
    int i, length = 0;
    for (i=0; i<N; i++ ) {
	if (x2y[i] >= 0 ) length ++;
    }
    return length;
}
/********************************/
/********************************/
/********************************/
int map_complementarity (Map *map1, Map *map2,  double *z) {
    
    int i, j, size = 0;
    int * map[2];
    double **array;
    double **image[2];
    double complementarity = 0.0, z_c,avg, avg_sq, stdev;
    double upper, lower;
    int ctr, length;

    size =  map1->x2y_size +  map2->y2x_size;
    array = dmatrix(2, size);
    if (! array ) return 1;
    
    upper = -100;
    lower = 100;
    ctr = 0;

    map[0]  =  map1->x2y;
    map[1]  =  map2->x2y;
    image[0] = map1->image;
    image[1] = map2->image;
    
	
    length = 0;
    for (i=0; i < map1->x2y_size ; i++ ) {
	if ( map[0][i] < 0 &&  map[1][i] < 0 ) continue;
	
	for (ctr=0; ctr<2; ctr++ ) {
	    if ( (j=map[ctr][i]) >= 0 ) {
		array[ctr][length] = image[ctr][i][j];
		if (upper < array[ctr][length] ) {
		    upper =  array[ctr][length];
		} else if ( lower > array[ctr][length] ) {
		    lower =  array[ctr][length];
		}
	    } else {
		array[ctr][length] = -10; 
	    }
	}
	length++;	
    }
    for (i=0; i<length; i++ ) {
	for (ctr=0; ctr<2; ctr++ ) {
	    if ( array[ctr][i] < -1 ) {
		array[ctr][i] = lower;
	    }
	}
    }

    /* printf (">>>>>>  l, u %6.2lf %6.2lf\n", lower, upper); */
   
    if ( length > 3) {
	double diff;
	double avg_1, avg_2, avg_3, avg_4;
	double u = upper, l = lower;
	double interval = u-l;
	int N = length;
	
	complementarity = 0.0;
	for (i=0; i < length; i++ ) {
	    /* printf (">>>   %2d    %6.2lf %6.2lf\n", i, */
	    /* 	    array[0][i], array[1][i]); */
	    diff = array[0][i] - array[1][i] ;
	    complementarity += diff*diff;
	}
	
	avg_1 = (u*u - l*l)/2/interval;
	avg_2 = (u*u*u - l*l*l)/3/interval;
	avg_3 = (u*u*u*u - l*l*l*l)/4/interval;
	avg_4 = (u*u*u*u*u - l*l*l*l*l)/5/interval;
   
	avg = 2*N*( avg_2 - avg_1*avg_1);
	avg_sq = N*(N-1)*( 4*avg_2*avg_2 - 8*avg_2*avg_1*avg_1 +
			   4*avg_1*avg_1*avg_1*avg_1)
	    + N*(2*avg_4 - 8*avg_3*avg_1 + 6*avg_2*avg_2);
	stdev = sqrt (avg_sq - avg*avg);
	z_c = (complementarity-avg)/stdev;
	/* printf (" complementarity =   %6.2lf   z =   %6.2lf \n", */
        	  /*   complementarity/length, z_c); */
	*z = z_c;
    } else {
	*z = -10;
    }

    
    //if ( length > 3) kendl1(array1, array2, length, tau, z, prob);

    free_dmatrix(array);
    
    return  0;
    
    
}
/********************************/
/********************************/
/********************************/
/*
int map_consistence (  int NX, int NY,  Map *map1, Map *map2,
		       double *total_ptr,  double *gap_score, FILE *fptr) {
*/
    
int map_consistence (  int NX, int NY,  Map * combined_map, Map *map1, Map *map2,
		       double *total_ptr,  double *gap_score, FILE *fptr) {
    int i,j;
    double val1, val2;
    double total = 0;
    double aln_score;
/*
    static Map combined_map = {0};
*/

/*

    if ( !map1  || !map2) { // intialize and return
	if ( combined_map.x2y ) free_map (&combined_map);
	if (NX && NY) if ( initialize_map (&combined_map, NX, NY) ) exit (1);
	return 0;
    }
*/

    
    if (!NX) NX = map1->x2y_size; /* TODO: rename */
    if (!NY) NY = map1->y2x_size; /* TODO: rename */
    
    for (i=0; i<NX; i++) {
	for (j=0; j<NY; j++) {
	    val1 =  map1->image[i][j];
	    val2 =  map2->image[i][j];

	    if ( val1 > val2) {
		//combined_map.image[i][j]  = val1;
                combined_map->image[i][j]  = val1;
		//combined_map.cosine[i][j] = map1->cosine[i][j];
                combined_map->cosine[i][j] = map1->cosine[i][j];
	    } else {
		//combined_map.image[i][j]  = val2;
                combined_map->image[i][j]  = val2;
		//combined_map.cosine[i][j] = map2->cosine[i][j];
                combined_map->cosine[i][j] = map2->cosine[i][j];
	    }
	}
    }
      

    /* Needleman on combined image */
/*
    needleman_wunsch (NX, NY, combined_map.image,
		      combined_map.x2y,  combined_map.y2x, &aln_score );
*/
    needleman_wunsch (NX, NY, combined_map->image,
		      combined_map->x2y,  combined_map->y2x, &aln_score );
    
    /* how do lengths compare? */
    //combined_map.matches = match_length (NX, combined_map.x2y);
    combined_map->matches = match_length (NX, combined_map->x2y);
    

/*
    if ( ! combined_map.matches ) {
	combined_map.assigned_score = 0.0;
    } else {
	for (i=0; i<NX; i++) {
	    j = combined_map.x2y[i];
	    if (j<0) continue;
	    total += combined_map.image[i][j];
	}
	combined_map.assigned_score = total;
    }
*/
    
        if ( ! combined_map->matches ) {
	combined_map->assigned_score = 0.0;
    } else {
	for (i=0; i<NX; i++) {
	    //j = combined_map.x2y[i];
            j = combined_map->x2y[i];
            
	    if (j<0) continue;
	    total += combined_map->image[i][j];
	}
	combined_map->assigned_score = total;
    }
    //if (gap_score) check_gap_lengths (&combined_map, gap_score);
    /* find z-score for each map separately */
    
     
    if (fptr) {
	fprintf (fptr," map 1 matches %d \n", map1->matches );
	fprintf (fptr," map 2 matches %d \n", map2->matches );
	//fprintf (fptr," combined  matches %d \n", combined_map.matches);
        fprintf (fptr," combined  matches %d \n", combined_map->matches);
	if (gap_score) 
	fprintf (fptr," gap score %8.3lf\n", *gap_score);
	//print_map (fptr, &combined_map, NULL, NULL,  NULL, NULL, 0);
        print_map (fptr, combined_map, NULL, NULL,  NULL, NULL, 0);
	
	if (0) {
	    fprintf (fptr,"image:\n");
	    for (i=0; i<NX; i++) {
		for (j=0; j<NY; j++) {
		    //if (  combined_map.image[i][j] < -9 ) {
                    if (  combined_map->image[i][j] < -9 ) {
			fprintf (fptr,"%6s ", "  X  " );
		    } else {
			//fprintf (fptr,"%6.2lf ", combined_map.image[i][j]);
                        fprintf (fptr,"%6.2lf ", combined_map->image[i][j]);
		    }
		}
		fprintf (fptr,"\n");
	    }
	    fprintf (fptr,"aln score: %6.2lf\n", aln_score);
	}
	
    }
    /* pass NULL for total_ptr  if I don't need it */
    //if (total_ptr) *total_ptr =   combined_map.assigned_score ; 
    if (total_ptr) *total_ptr =   combined_map->assigned_score ; 
    
    return  0;
}
/***************************************/
/***************************************/
