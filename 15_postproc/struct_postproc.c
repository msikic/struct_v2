# include "struct.h"


//# define DEBUG

int smith_waterman_2 (int max_i, int max_j, double **similarity,
			  int *map_i2j, int * map_j2i, double * aln_score);

# define MAX_DIST_TO_CONSIDER 9.0

/********************************/
int postprocess (Descr *descr1, Protein * protein1, Representation *rep1, 
		 Descr *descr2, Protein * protein2, Representation *rep2, 
                 Map *map, Score * score ){

    int no_res_1 = protein1->length, no_res_2= protein2->length;
    int resctr1, resctr2;
    int element_ctr_1, element_ctr_2;
    int *element_1_begin, *element_1_end; /* "element" here means SSE */
    int *element_2_begin, *element_2_end;
    int *element_1_begin_pdb, *element_1_end_pdb; /* the labels that the beginnings and ends of an element have in pdb*/
    int *element_2_begin_pdb, *element_2_end_pdb;
    int num_pdb_id;
    int map_size;
    int *residue_map_i2j, *residue_map_j2i;
    
    int last_res_prev_loop_1, last_res_prev_loop_2; 
    int first_res_prev_loop_1, first_res_prev_loop_2;
    int last_res_1, last_res_2; 
    int first_res_1, first_res_2;
    int last_res_next_loop_1, last_res_next_loop_2; 
    int first_res_next_loop_1, first_res_next_loop_2;
    //int last_element_1, last_element_2;
    int type, *type_1, *type_2;
    
    char tmp[PDB_ATOM_RES_NO_LEN+1] = {'\0'};
    double d, d0 = 10.0;
    double aln_score, rmsd;
    double ca1[3], ca2[3], rotated_ca1[3];
    double ** similarity;
    double **x, **y;
    double **R, T[3], q[4];
    double total_score = 0;
  
    Element * element;
    /* for the MC: */
    int max_no_steps = 20, no_steps = 0;
    int done = 0;
    // int toggle = 0;
    double *current_q,  *old_q, *current_T, *old_T;
    double *best_q, *best_T;
    double old_score = total_score, current_score = 0.0, d_mc = 0.5;
    double max_score;
    /// double t_mc, 
    double d_init;
    double aux; // step = NR_POINTS/MAX_EXP_VALUE;
   
     
    double alignment_score (Protein * protein1, Protein * protein2, int * residue_map_i2j,
			    double **R, double *T,  double d0 );
    int  alignment_size  (int * residue_map_i2j, int no_res_1 );
    int  closeness_score (Descr *descr1, Representation *rep1, Representation *rep2, Map * map,
			  int *element_1_begin, int *element_1_end,
			  int *element_2_begin, int *element_2_end,
			  Protein *protein1, Protein *protein2,
			  double **R, double *T, double d0, double ** similarity, double * score_ptr);
    int following_loop (int *element_begin, int *element_end,
			int no_of_elements, int no_of_res, 
			int element_ctr, int * first_res, int * last_res);  
    int map2rotation (Protein *protein1, Protein *protein2, int *residue_map_i2j,
		       double **x, double **y, double *q, double *T, double *rmsd);
    int preceding_loop (int *element_begin, int *element_end,
			int element_ctr, int * first_res, int * last_res);
    
    if ( ! (R=dmatrix(3,3) ) ) return 1; /* compiler is bugging me otherwise */
    
    construct_translation_vecs (rep1, rep2, map);
    
    /* make sure that we have all the info we might need */
  
    /* define matrix, the size of nr of residues in set of SSEs 
       x nr of  residues in the other set of SSEs, and fill it with -1 */
    similarity = dmatrix (no_res_1, no_res_2);
    if ( !similarity ) return 1;
    
    for (resctr1=0; resctr1<no_res_1; resctr1++) {
	for (resctr2=0; resctr2<no_res_2; resctr2++) {
	    similarity[resctr1][resctr2] = -1;
	}
    }
    /* alloc */
    if ( ! (element_1_begin = emalloc (no_res_1*sizeof(int) )) ) return 1;
    if ( ! (element_1_end   = emalloc (no_res_1*sizeof(int) )) ) return 1;
    
    if ( ! (element_2_begin = emalloc (no_res_2*sizeof(int) )) ) return 2;
    if ( ! (element_2_end   = emalloc (no_res_2*sizeof(int) )) ) return 2;
    
    if ( ! (element_1_begin_pdb = emalloc (no_res_1*sizeof(int) )) ) return 1;
    if ( ! (element_1_end_pdb   = emalloc (no_res_1*sizeof(int) )) ) return 1;
    
    if ( ! (element_2_begin_pdb = emalloc (no_res_2*sizeof(int) )) ) return 2;
    if ( ! (element_2_end_pdb   = emalloc (no_res_2*sizeof(int) )) ) return 2;
    
    if ( ! (type_1   = emalloc (no_res_1*sizeof(int) )) ) return 1;
    if ( ! (type_2   = emalloc (no_res_2*sizeof(int) )) ) return 2;
    

    if ( ! (residue_map_i2j = emalloc (no_res_1*sizeof(int) )) ) return 1;
    if ( ! (residue_map_j2i = emalloc (no_res_2*sizeof(int) )) )  return 2;

    if ( ! (x = dmatrix (3, no_res_1+no_res_2)))  exit(1);
    if ( ! (y = dmatrix (3, no_res_1+no_res_2)))  exit(1);
    
    /*********************************************************************/
    /*********************************************************************/
    /* some bookkeeping */
    /* find the beginning and the end of SSEs according to the pdb_tag   */
    for (element_ctr_1=0; element_ctr_1 < descr1->no_of_elements; element_ctr_1++) {
	element = descr1->element+element_ctr_1;
	/* get rid of the insertion tag, in case we
	   haven't shaken it off along the way */
	memcpy (tmp,  element->begin_id, PDB_ATOM_RES_NO_LEN*sizeof(char) );
	element_1_begin_pdb[element_ctr_1] = atoi (tmp);
	
	memcpy (tmp,  element->end_id, PDB_ATOM_RES_NO_LEN*sizeof(char) );
	element_1_end_pdb[element_ctr_1] = atoi (tmp);
    }
   for (element_ctr_1=0; element_ctr_1 < descr1->no_of_elements; element_ctr_1++) {
	printf ("orig  %3d  %3d  %3d \n", element_ctr_1 , element_1_begin_pdb[element_ctr_1], element_1_end_pdb[element_ctr_1]);
    }


    for (element_ctr_2=0; element_ctr_2 < descr2->no_of_elements; element_ctr_2++) {
	element = descr2->element+element_ctr_2;
	/* get rid of the insertion tag, in case we
	   haven't shaken if off along the way */
	memcpy (tmp,  element->begin_id, PDB_ATOM_RES_NO_LEN*sizeof(char) );
	element_2_begin_pdb[element_ctr_2] = atoi (tmp);
	
	memcpy (tmp,  element->end_id, PDB_ATOM_RES_NO_LEN*sizeof(char) );
	element_2_end_pdb[element_ctr_2] = atoi (tmp);
    }
    
    /* translate the beginning and the end of SSEs from pdb_tag to the position on the sequence*/
    element_ctr_1 = 0;
    for (resctr1=0; resctr1<no_res_1; resctr1++) {
	
	memcpy (tmp, protein1->sequence[resctr1].pdb_id, PDB_ATOM_RES_NO_LEN*sizeof(char) );
	num_pdb_id = atoi (tmp);
	
	/* break if we are outside of the last element */
	if ( num_pdb_id > element_1_end_pdb [descr1->no_of_elements-1]) {
	    /*see for example 1k4j - the first residue has the number 399 */
	    continue;
	}
	
	if (  num_pdb_id == element_1_end_pdb[element_ctr_1] ) {
	    element_1_end [element_ctr_1] = resctr1;
	    element_ctr_1++;
	}
	
	if ( num_pdb_id == element_1_begin_pdb[element_ctr_1] ) {
	    element_1_begin [element_ctr_1] = resctr1;
	}
	/* also if there is some disagreement about where the first element starts: */
	if ( ! resctr1 && num_pdb_id >= element_1_begin_pdb[element_ctr_1] ) {
	    element_1_begin [element_ctr_1] = resctr1;
	}
    }	    

  
    element_ctr_2 = 0;
    for (resctr2=0; resctr2<no_res_2; resctr2++) {
	
	memcpy (tmp, protein2->sequence[resctr2].pdb_id, PDB_ATOM_RES_NO_LEN*sizeof(char) );
	num_pdb_id = atoi (tmp);

	/* break if we are outside of the last element */
	if ( num_pdb_id > element_2_end_pdb [descr2->no_of_elements-1]) break;
	
	if (  num_pdb_id==element_2_end_pdb  [element_ctr_2] ) {
	    element_2_end [element_ctr_2] = resctr2;
	    element_ctr_2++;
	}
	
	if ( num_pdb_id == element_2_begin_pdb[element_ctr_2] ) {
	    element_2_begin [element_ctr_2] = resctr2;
	}
	/* also if there is some disagreement about where the first element starts: */
	if ( ! resctr2 && num_pdb_id >= element_2_begin_pdb[element_ctr_2] ) {
	    element_2_begin [element_ctr_2] = resctr2;
	}
    }
   
    /*****************************************************************/
    /* we'll allow for a bit of sloppiness here to accomodaite cath's
       mis-definition of domains:  */
    if ( element_1_end[descr1->no_of_elements-1] >= protein1->length) {
	 element_1_end[descr1->no_of_elements-1]  = protein1->length-1;
    }
    if ( element_2_end[descr2->no_of_elements-1] >= protein2->length) {
	 element_2_end[descr2->no_of_elements-1]  = protein2->length-1;
    }
    
  
   
    /*****************************************************************/
    /* make the information about SSE type available by residue number
       - we'll use it below */
    memset (type_1, 0, no_res_1*sizeof(int) );
    for (element_ctr_1=0; element_ctr_1<descr1->no_of_elements; element_ctr_1++) {
	type = descr1->type[element_ctr_1];
	for (resctr1=element_1_begin[element_ctr_1]; resctr1<= element_1_end[element_ctr_1]; resctr1++) {
	    type_1[resctr1] = type;
	}
    }

    memset (type_2, 0, no_res_2*sizeof(int) );
    for (element_ctr_2=0; element_ctr_2<descr2->no_of_elements; element_ctr_2++) {
	type = descr2->type[element_ctr_2];
	for (resctr2=element_2_begin[element_ctr_2]; resctr2 <= element_2_end[element_ctr_2]; resctr2++) {
	    type_2[resctr2] = type;
	}
    }

    /************************************************************/
    /************************************************************/
    /* ALIGNMENT, round 1                                       */
    /************************************************************/
    /* for all mapped blocks calculate similarity as exp (-d/d0) */
    total_score = 0.0;
    quat_to_R (map->q, R);

    closeness_score (descr1, rep1, rep2, map, element_1_begin, element_1_end,
		     element_2_begin, element_2_end, protein1, protein2,
		     R, NULL, d0, similarity, &total_score);
	    

    /* run Smith-Waterman and use the mapped CA to find the transformation        */
    /* I have another copy of SW here (the first one is in struct_map.c)          */
    /* so I wouldn't fumble with parameters - the two should be joined eventually */
    smith_waterman_2 (no_res_1, no_res_2, similarity,
		      residue_map_i2j,  residue_map_j2i, &aln_score);

    map2rotation (protein1, protein2, residue_map_i2j, x, y, q, T, &rmsd);
    
    quat_to_R (q, R);
    current_score = alignment_score (protein1, protein2, residue_map_i2j, R, T, d0);

# ifdef DEBUG
    int i,j;
    
    printf (" before MC:  total: %8.3lf   aln: %8.3lf   aln fixed: %8.3lf  \n\n",
    		    total_score, aln_score, current_score);
    for (i=0; i<3; i++) {
	for (j=0; j<3; j++) fprintf (stdout,"%8.3lf ", R[i][j]);
	fprintf (stdout,"%8.3lf  \n", T[i]);
    }

# endif

    /*********************************************************/
    /* fiddle iteratively with the transformation            */
    if ( ! (current_q = emalloc (4*sizeof(double)) )) return 1;
    if ( ! (old_q     = emalloc (4*sizeof(double)) )) return 1;
    if ( ! (best_q    = emalloc (4*sizeof(double)) )) return 1;
    if ( ! (current_T = emalloc (3*sizeof(double)) )) return 1;
    if ( ! (old_T     = emalloc (3*sizeof(double)) )) return 1;
    if ( ! (best_T    = emalloc (3*sizeof(double)) )) return 1;
    
    srand48 (time (0));
    
    memcpy (current_q, q, 4*sizeof(double));
    memcpy (    old_q, q, 4*sizeof(double));
    memcpy (   best_q, q, 4*sizeof(double));
    
    memcpy (    old_T, T, 3*sizeof(double));
    memcpy (   best_T, T, 3*sizeof(double));
    memcpy (current_T, T, 3*sizeof(double));

    quat_to_R ( current_q, R);

    d_init = d0;


    /* t_mc = exp ( (1.0- (double)anneal_round)/10.0); */
    d_mc = d_init;
//    t_mc = 5;
	
    memcpy (current_q, best_q, 4*sizeof(double));
    memcpy (current_T, best_T, 3*sizeof(double));
    memcpy (old_q, best_q, 4*sizeof(double));
    memcpy (old_T, best_T, 3*sizeof(double));
    old_score = 0;
    max_score = 0;
    no_steps  = 0;
    
//    toggle = 1;
    done = 0;
    
    while (no_steps < max_no_steps && !done ) {

	    
	closeness_score (descr1, NULL, NULL, map, element_1_begin, element_1_end,
			 element_2_begin, element_2_end, protein1, protein2,
			 R, current_T, d_mc, similarity, &total_score);
	    
	
	smith_waterman_2 (no_res_1, no_res_2, similarity, residue_map_i2j, residue_map_j2i, &aln_score);

	map2rotation (protein1, protein2, residue_map_i2j, x, y, current_q,  current_T, &rmsd);
	
	quat_to_R ( current_q, R);
	current_score = alignment_score (protein1, protein2, residue_map_i2j, R, current_T, d_mc);


	//printf (" step %3d:  total: %8.3lf   aln: %8.3lf   aln fixed: %8.3lf  \n",
	//	no_steps+1,  total_score, aln_score, current_score);
	
	if ( current_score >  max_score )  {
	    max_score = current_score;
	    memcpy (best_q, current_q, 4*sizeof(double));
	    memcpy (best_T, current_T, 3*sizeof(double));
	}

	if (old_score) done = ( fabs(old_score-current_score)/old_score < 0.01);
	old_score = current_score;
	no_steps++;
    
    }

    
    memcpy (q, best_q, 4*sizeof(double));
    memcpy (T, best_T, 3*sizeof(double));
	 
    

    free (current_q);
    free (old_q);
    free (best_q);
    free (current_T);
    free (old_T);
    free (best_T);
   


   
    /************************************************************/
    /************************************************************/
    /* ALIGNMENT, round 2                                       */
    /************************************************************/
    /************************************************************/
    /* find the similarity matrix for this new rotation
       -- this time extending to neighboring elements */
    for (resctr1=0; resctr1<no_res_1; resctr1++) {
	for (resctr2=0; resctr2<no_res_2; resctr2++) {
	    similarity[resctr1][resctr2] = -1;
	}
    }
   // last_element_1 = 0;
    // last_element_2 = 0;
    for (element_ctr_1=0; element_ctr_1 < descr1->no_of_elements; element_ctr_1++) {

	
	element_ctr_2 = map->x2y[element_ctr_1];
	if (element_ctr_2 <  0) continue;

	// last_element_1 = element_ctr_1;
	// last_element_2 = element_ctr_2;
	


	/****************************************************/
	/* try to extend the match to the previous loop     */

	/* the first and the last residues of the preceding loops */
	if ( element_1_begin[element_ctr_1] == 0 ) {
	     first_res_prev_loop_1 = 0;
	} else {
	     preceding_loop (element_1_begin, element_1_end, element_ctr_1,
			&first_res_prev_loop_1, &last_res_prev_loop_1);
	}
	
	if ( element_2_begin[element_ctr_2] == 0 ) {
	     first_res_prev_loop_2 = 0;
	} else {
	     preceding_loop (element_2_begin, element_2_end, element_ctr_2,
			&first_res_prev_loop_2, &last_res_prev_loop_2);
	}

	/****************************************************/
	/* now try to extend the match to the following loop */
	/* the first and the last residues of the following loops */
	/* there is some redundancy here (prev/next) but we can */
        /* get rid of it when it becomes a bottleneck */
	if ( element_1_end[element_ctr_1] >= no_res_1 -1) {
	     last_res_next_loop_1 = no_res_1 -1;
	} else {
	     following_loop (element_1_begin, element_1_end, descr1->no_of_elements, no_res_1,
			element_ctr_1, &first_res_next_loop_1, &last_res_next_loop_1);
	}
	
	if ( element_2_end[element_ctr_2] >= no_res_2 -1) {
	     last_res_next_loop_2 = no_res_2 -1;
	} else {
	     following_loop (element_2_begin, element_2_end, descr2->no_of_elements, no_res_2,
			element_ctr_2, &first_res_next_loop_2, &last_res_next_loop_2);
	}

	/* this is an exact map point  -- these two we won't try aligning to anything else */
	for (resctr1 = first_res_prev_loop_1; resctr1<= last_res_next_loop_1; resctr1++) {
	    /* find the representative atom (CA) for residue resctr1*/
	    if ( find_Calpha ( protein1, resctr1, ca1 ) ) {
		/* find_Calpha returns 1 on failure */
		for (resctr2= first_res_prev_loop_2; resctr2<= last_res_next_loop_2; resctr2++)
		    similarity[resctr1][resctr2] = 0.0;
		continue;
	    }
	    /* rotate and  translate*/
	    point_rot_tr (ca1, R, T, rotated_ca1);
		
	    for (resctr2= first_res_prev_loop_2; resctr2<= last_res_next_loop_2; resctr2++) {
		/* skip if we are sure the two are not of the same type */
		if ( (type_1[resctr1]|type_2[resctr2]) == (HELIX|PARALLEL) )  continue;
	

		/* find the representative atom (CA) for residue resctr2*/
		if ( find_Calpha (protein2, resctr2,  ca2 ) ) {
		    similarity[resctr1][resctr2] = 0.0;
		    continue;
		}
		d = two_point_distance (rotated_ca1, ca2);
		
		if ( d<MAX_DIST_TO_CONSIDER ) { 
		    aux = d/d0;
		    if (  aux <  MAX_EXP_VALUE ) {
			//similarity[resctr1][resctr2] = exp_table[(int)(aux*step)];
			similarity[resctr1][resctr2] = exp (-aux);
		    } else {
			similarity[resctr1][resctr2] = 0.0;
		    }
		}
	    }
	}



    }
    
    /* another ad-hoc fix: did we skip an element of the
       same type in both sequences? */
    int gap = 0;
    int gap_begin_1 = 0, gap_end_1 = 0;
    int gap_begin_2 = 0, gap_end_2 = 0;
    int all_healthy;
    for (element_ctr_1=0; element_ctr_1 < descr1->no_of_elements; element_ctr_1++) {
     	element_ctr_2 = map->x2y[element_ctr_1];

	if (element_ctr_2>=0) {
	  if (gap && gap <=2) {
	    gap_end_1 = element_ctr_1-1;
	    gap_end_2 = element_ctr_2-1;

	    all_healthy  = (gap_begin_1 >=0);
	    all_healthy &= (gap_end_1 >=0);
	    all_healthy &= (gap_begin_2 >=0);
	    all_healthy &= (gap_end_2 >=0);
	    all_healthy &= ( gap_end_1 - gap_begin_1 <= 2);
	    all_healthy &= ( gap_end_2 - gap_begin_2 <= 2);

	    if ( all_healthy ) {
		first_res_1  = element_1_begin[gap_begin_1];
		last_res_1   = element_1_end  [gap_end_1];
		first_res_2  = element_2_begin[gap_begin_2];
		last_res_2   = element_2_end  [gap_end_2];
		for (resctr1=first_res_1; resctr1<=last_res_1; resctr1++) {
		    if ( find_Calpha ( protein1, resctr1, ca1 ) ) {
			/* find_Calpha returns 1 on failure */
			for (resctr2= first_res_prev_loop_2; resctr2<= last_res_next_loop_2; resctr2++)
			    similarity[resctr1][resctr2] = 0.0;
			continue;
		    }
		    /* rotate and  translate*/
		    point_rot_tr (ca1, R, T, rotated_ca1);
		    /* the two, however, cannot be matched if they are not
		       of the same type (helix, strand or unassigned)*/
		    for (resctr2=first_res_2; resctr2<=last_res_2; resctr2++) {
		
			if ( (type_1[resctr1]|type_2[resctr2]) == (HELIX|PARALLEL) )  continue;
			
			if ( find_Calpha (protein2, resctr2,  ca2 ) ) {
			    similarity[resctr1][resctr2] = 0.0;
			    continue;
			}

			d = two_point_distance (rotated_ca1, ca2);
			if ( d<MAX_DIST_TO_CONSIDER ) { 
			    aux = d/d0;
			    if (  aux <  MAX_EXP_VALUE ) {
				//similarity[resctr1][resctr2] = exp_table[(int)(aux*step)];
				similarity[resctr1][resctr2] = exp (-aux);
			    } else {
				similarity[resctr1][resctr2] = 0.0;
			    }
			}
		    }
		}
	    }
	    gap = 0;
	    gap_begin_1 = element_ctr_1+1;
	    gap_begin_2 = element_ctr_2+1;
	    continue;
	  }
	  gap++;
	}
    }
    
    /***********************************************/
    /***********************************************/
    /***********************************************/
    /***********************************************/
    
    memset (residue_map_i2j, 0, no_res_1*sizeof(int)); 
    memset (residue_map_j2i, 0, no_res_2*sizeof(int)); 
    smith_waterman_2 (no_res_1, no_res_2, similarity,
		      residue_map_i2j,  residue_map_j2i, &aln_score);
    map2rotation (protein1, protein2, residue_map_i2j, x, y, q, T, &rmsd);
    
    quat_to_R (q, R);
 
    
# ifdef DEBUG
    printf ("in postp function: rmsd %8.4lf \n", rmsd );

    printf ("size  %4d  %4d\n", no_res_1, no_res_2);
    printf ("residue-level map: \n");
    for (resctr1=0; resctr1<no_res_1; resctr1++) {
	find_Calpha ( protein1, resctr1, ca1 );
	point_rot_tr (ca1, R, T, rotated_ca1);

	resctr2 = residue_map_i2j[resctr1];

	if ( resctr2 < 0 ) continue;
	
	if (find_Calpha (protein2, resctr2,  ca2)) continue;

	d = two_point_distance (rotated_ca1, ca2);
	printf ("%3d --> %3d   %8.4lf  %8.4lf    \n", resctr1, resctr2, d, similarity[resctr1][resctr2] );
   }
    exit (1);
# endif
    
    aln_score = alignment_score (protein1, protein2, residue_map_i2j, R, T, d0);
    map_size  = alignment_size (residue_map_i2j, protein1->length);		

    memcpy (&(map->q[0]), &q[0], 4*sizeof(double) );
    memcpy (&(map->T[0]), &T[0], 3*sizeof(double) );
    
    map->x2y_residue_level  = residue_map_i2j;
    map->y2x_residue_level  = residue_map_j2i;
    
    map->x2y_residue_l_size = no_res_1;
    map->y2x_residue_l_size = no_res_2;
    

    /*************************************************************************/
    map->res_almt_length     = map_size;
    map->res_rmsd            = rmsd;

    score->res_almt_score  = aln_score;
    score->res_almt_length = map_size;
    score->res_rmsd            = rmsd;

    free_dmatrix(R);
    free_dmatrix (similarity);
    free_dmatrix (x);
    free_dmatrix (y);
    free (element_1_begin);
    free (element_2_begin);
    free (element_1_end);
    free (element_2_end);
    free (element_1_begin_pdb);
    free (element_2_begin_pdb);
    free (element_1_end_pdb);
    free (element_2_end_pdb);
    free (type_1);
    free (type_2);

 
     
    return 1;
}


/*************************************************************************/
/*************************************************************************/
int  alignment_size ( int * residue_map_i2j, int no_res_1 ) {

    int resctr1, resctr2;
    int aln_size =0;
    
    for (resctr1=0; resctr1<no_res_1; resctr1++) {
	resctr2 = residue_map_i2j[resctr1];
	if (resctr2 < 0 ) continue;
	
	aln_size++;
    }
    return aln_size;
}

/*************************************************************************/
/*************************************************************************/
double alignment_score (Protein * protein1, Protein * protein2, int * residue_map_i2j,
			double **R, double *T,  double d0 ) {

    double d;
    double aln_score = 0.0;
    double ca1[3], ca2[3], rotated_ca1[3];
    int resctr1, resctr2;
    int no_res_1 = protein1->length;
    int aln_size =0;
    
    for (resctr1=0; resctr1<no_res_1; resctr1++) {
	resctr2 = residue_map_i2j[resctr1];
	if (resctr2 < 0 ) continue;
	
	find_Calpha ( protein1, resctr1, ca1 );
	find_Calpha ( protein2, resctr2, ca2 );
	
	/* rotate and  translate*/
	point_rot_tr (ca1, R, T, rotated_ca1);
		
	d = two_point_distance (rotated_ca1, ca2);
	aln_score +=  exp (-d/d0);

	//printf ("\t  %3d --> %3d    %8.2lf   %8.2lf   \n",
	//	resctr1, resctr2, d,  exp (-d/d0));
	aln_size++;
    }
    return aln_score;
}


/*************************************************************************/
/*************************************************************************/
int  closeness_score (Descr *descr1, Representation *rep1,  Representation *rep2, Map *map,
		      int *element_1_begin, int *element_1_end,
		      int *element_2_begin, int *element_2_end,
		      Protein *protein1, Protein *protein2,
		      double **R, double *fixed_T, double d0, double ** similarity, double * score_ptr) {

    int element_ctr_1, element_ctr_2;
    int resctr1, resctr2, i, j;
    double d, score;
    double T[3];
    double ca1[3], ca2[3], rotated_ca1[3];
    double aux; // step = NR_POINTS/MAX_EXP_VALUE;
    /***********/
    int find_Calpha ( Protein *protein, int  resctr, double ca[3] );
    
    if ( !rep1 )  {
	if ( !fixed_T) {
	    fprintf (stderr, "Error in closeness_score().\n");
	    exit (1);
	}
	memcpy (T, fixed_T, 3*sizeof(double));
    } else if (!rep2 ) {
	fprintf (stderr, "Error in closeness_score() -- both reps or none ...\n");
	exit (1);
    }
  
    /* for all mapped blocks calculate similarity as exp (-d/d0) */
    score = 0.0;
    
    for (element_ctr_1=0; element_ctr_1 < descr1->no_of_elements; element_ctr_1++) {
	
	 element_ctr_2 = map->x2y[element_ctr_1];
	 if (element_ctr_2 < 0) continue;

	 if ( rep1 ) { /* "exploding" from the origin */
	     for (i=0; i<3; i++) {
		 T[i] = 0.0;
		 for (j=0; j<3; j++) T[i] += R[i][j]*rep1->translation[element_ctr_1][j];
	     }
	 }

# ifdef DEBUG
	 printf ( "elmt %3d  T:  %8.3lf  %8.3lf  %8.3lf \n", element_ctr_1, T[0], T[1], T[2]);
	           printf ( "\t   elmt2  %3d  T:  %8.3lf  %8.3lf  %8.3lf \n", element_ctr_2,
			    rep2->translation[element_ctr_2][0],
			    rep2->translation[element_ctr_2][1],
			    rep2->translation[element_ctr_2][2]);
# endif
	 
	 for (resctr1=element_1_begin[element_ctr_1]; 
	      resctr1<= element_1_end[element_ctr_1]; resctr1++) {
	    
	      /* find the representative atom (CA) for residue resctr1*/
	      find_Calpha ( protein1, resctr1, ca1 );
	    

              if ( rep1) {/* translate to the origin */
		  for (i=0; i<3; i++) ca1[i] -= rep1->origin[i];
	      }
	      
	      /* rotate and  translate out*/
	      point_rot_tr (ca1, R, T, rotated_ca1);
	    
	      for (resctr2=element_2_begin[element_ctr_2]; 
		   resctr2<= element_2_end[element_ctr_2]; resctr2++) {
			/* find the representative atom (CA) for residue resctr2*/
		   find_Calpha (protein2, resctr2,  ca2 );
                   /* translate to cm & out */
		   if ( rep2 ) {
		       for (i=0; i<3; i++) ca2[i] += -rep2->origin[i] + rep2->translation[element_ctr_2][i];
		   }
#		   /* finally, find the distance & assign the "similarity" score */
		   d = two_point_distance (rotated_ca1, ca2);
		   if ( d<MAX_DIST_TO_CONSIDER ) { 
		       aux = d/d0;
		       if (  aux <  MAX_EXP_VALUE ) {
			   //similarity[resctr1][resctr2] = exp_table[(int)(aux*step)];
			   similarity[resctr1][resctr2] = exp (-aux);
		       } else {
			   similarity[resctr1][resctr2] = 0.0;
		       }
		   }

		   score += similarity[resctr1][resctr2];
# ifdef DEBUG
		   printf ( "\t %3d  %3d     %8.3lf \n",  resctr1, resctr2, d);
# endif
	      }
	 }

    }

    *score_ptr = score;
    
    return 0;
}


/*************************************************************************/
/*************************************************************************/
int preceding_loop (int *element_begin, int *element_end,
		    int element_ctr, int * first_res, int * last_res) {
    
    if ( element_ctr == 0 ) {
	/* there is no previous element */
	*first_res = 0;
    } else {
	*first_res = element_end[element_ctr-1]+1;
    }
    *last_res = element_begin[element_ctr]-1;

    return 0;
}

/*************************************************************************/
/*************************************************************************/
int following_loop (int *element_begin, int *element_end,
		    int no_of_elements, int no_of_res, 
		    int element_ctr, int * first_res, int * last_res) {
    
    *first_res = element_end[element_ctr]+1;
    if (element_ctr >= no_of_elements-1) {
	*last_res = no_of_res - 1;
    } else {
	*last_res = element_begin[element_ctr+1]-1;
    }

    return 0;
}

/*************************************************************************/
/*************************************************************************/
double  two_point_distance (double point1[3], double point2[3] ) {
    int i;
    double aux;
    double d = 0;
    for (i=0; i<3; i++) {
	aux = point1[i] - point2[i];
	d  += aux*aux;
    }
    return  sqrt (d);
}
/*************************************************************************/
/*************************************************************************/
int point_rot_tr (double point_in[3], double **R, double T[3],double point_out[3]) {

    int i, j;
    /* rotate */
    for (i=0; i<3; i++) {
	point_out[i] = 0.0;
	for (j=0; j<3; j++) point_out[i] += R[i][j]*point_in[j];
    }
    /* translate out */
    for (i=0; i<3; i++) point_out[i] += T[i];
    return 0;
}
/*************************************************************************/
/*************************************************************************/
int smith_waterman_2 (int max_i, int max_j, double **similarity,
		      int *map_i2j, int * map_j2i, double * aln_score) {

    double **F; /*alignment_scoring table*/
    char   **direction;
    double gap_opening   = -0.2;
    double gap_extension = -0.1;
    double endgap        =  0.0;
    double penalty;
    double i_sim = 0.0, j_sim = 0.0, diag_sim = 0.0, max_sim = 0.0;
    double F_max;
    int use_endgap = 1;
    int F_max_i, F_max_j;
    int i,j;

     /* allocate F */
    if ( ! (F = dmatrix( max_i+1, max_j+1)) ) return 1;
    if ( ! (direction = chmatrix ( max_i+1, max_j+1)) ) return 1;


    
    /* fill the table */
    F_max = FAR_FAR_AWAY;
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
    for( i= max_i; i>F_max_i; i--) map_i2j[i-1] = FAR_FAR_AWAY;;
    for( j= max_j; j>F_max_j; j--) map_j2i[j-1] = FAR_FAR_AWAY;;
    i = F_max_i;
    j = F_max_j;
    *aln_score = F[i][j]; 
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
	    map_i2j [i-1] = FAR_FAR_AWAY;
	    i--; 
	    break; 
	case 'j':
	    //printf ( " %4d  %4d \n",  -1, j);
	    map_j2i [j-1] = FAR_FAR_AWAY;
	    j--; 
	    break; 
	default: 
	    fprintf (stderr, "Retracing error.\n");
		
	} 
    }

    /* free */ 
    free_dmatrix (F);
    free_cmatrix (direction);
    
    return 0; 
   
    
}
/************************************************************/
/************************************************************/
int find_Calpha ( Protein *protein, int  resctr, double ca[3] ){

    Residue * res = protein->sequence + resctr;

    if ( ! res ) return 1;
    if ( ! res->Ca ) return 1;
    
    ca[0]= res->Ca->x;
    ca[1]= res->Ca->y;
    ca[2]= res->Ca->z;


    return 0;

}

/************************************************************/
/************************************************************/



# if 0
    printf ("protein1 length: %d \n",no_res_1);
    descr_out (stdout, descr1);
    for (element_ctr_1=0; element_ctr_1 < descr1->no_of_elements; element_ctr_1++) {
	printf ("SSE %2d begins at %3d(%s)   ends at %3d (%s) \n",
		element_ctr_1,
		element_1_begin [element_ctr_1], protein1->sequence[element_1_begin [element_ctr_1]].pdb_id,
		element_1_end [element_ctr_1], protein1->sequence[element_1_end [element_ctr_1]].pdb_id);
	
    }
    printf ("protein2 length: %d \n",no_res_2);
    descr_out (stdout, descr2);
    for (element_ctr_2=0; element_ctr_2 < descr2->no_of_elements; element_ctr_2++) {
	printf ("SSE %2d begins at %3d(%s)   ends at %3d (%s) \n",
		element_ctr_2,
		element_2_begin [element_ctr_2], protein2->sequence[element_2_begin [element_ctr_2]].pdb_id,
		element_2_end [element_ctr_2], protein2->sequence[element_2_end [element_ctr_2]].pdb_id);
	
    }

    return 0;

    for (resctr1=0; resctr1< no_res_1; resctr1++) {
	printf (" %d  type %d \n", resctr1+1, type_1[resctr1]);
    }
    for (resctr2=0; resctr2< no_res_2; resctr2++) {
	printf (" %d  type %d \n", resctr2+1, type_2[resctr2] );
    }
    exit (1);
# endif
 
	
# if 0
	/* 	printf (" %4d: %5.2lf %5.2lf %5.2lf   %4d: %5.2lf %5.2lf %5.2lf  %8.2lf \n", */
/* 			resctr1, rotated_ca1[0], rotated_ca1[1],  rotated_ca1[2], */
/* 			resctr2, ca2[0], ca2[1], ca2[2], */
/* 			d); */


	d = 0;
	for (i=0; i<3; i++) {
	    aux = rotated_tr[i] - rep2->translation[element_ctr_2][i];
	    d += aux*aux;
	}


	printf ("elmts %2d--> %2d, cm distance= %6.2lf \n",
		element_ctr_1+1, element_ctr_2+1, sqrt(d));
# if 0
    for ( i =0; i < 3; i++ ) {
	for ( j =0; j < 3; j++ ) {
	    printf ("%8.3lf ",  R[i][j]);
	}
	printf ("%8.3lf \n", T[i]);
    }
    printf ("alignment score extended: %8.2f,   map size: %d,  rmsd: %8.2f\n\n", aln_score, map_size, rmsd);
# endif

# endif


