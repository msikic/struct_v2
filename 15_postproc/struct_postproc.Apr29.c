# include "struct.h"

int smith_waterman_2 (int max_i, int max_j, double **similarity,
			  int *map_i2j, int * map_j2i, double * aln_score);

/********************************/
int postprocess (Descr *descr1, Protein * protein1, Representation *rep1, 
		 Descr *descr2, Protein * protein2, Representation *rep2, 
                 Map *map, Score * score ){

    int no_res_1 = protein1->length, no_res_2= protein2->length;
    int resctr1, resctr2;
    int element_ctr_1, element_ctr_2;
    int *element_1_begin, *element_1_end; /* "element" here means SSE */
    int *element_2_begin, *element_2_end;
    int num_pdb_id;
    int i, j;
    int map_count, map_size;
    int *map_i2j, *map_j2i;
    
    int last_res_1, last_res_2; 
    int first_res_1, first_res_2;
    int last_element_1, last_element_2;
    int type, *type_1, *type_2;
    
    char tmp[PDB_ATOM_RES_NO_LEN+1] = {'\0'};
    double d, d0 = 10.0;
    double aln_score, rmsd;
    double ca1[3], ca2[3], rotated_ca1[3];
    double rotated_tr[3];
    double ** similarity;
    double **x, **y;
    double **R, T[3], q[4];
   
    Element * element;
    
    int following_loop (int *element_begin, int *element_end,
			int no_of_elements, int no_of_res, 
			int element_ctr, int * first_res, int * last_res);  
    int map2rotation (double **x, double **y, int map_size,
		      double q[4], double T[3], double *rmsd);
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
    
    if ( ! (type_1   = emalloc (no_res_1*sizeof(int) )) ) return 1;
    if ( ! (type_2   = emalloc (no_res_2*sizeof(int) )) ) return 2;
    

    if ( ! (map_i2j = emalloc (no_res_1*sizeof(int) )) ) return 1;
    if ( ! (map_j2i = emalloc (no_res_2*sizeof(int) )) )  return 2;

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
	element_1_begin[element_ctr_1] = atoi (tmp);
	
	memcpy (tmp,  element->end_id, PDB_ATOM_RES_NO_LEN*sizeof(char) );
	element_1_end[element_ctr_1] = atoi (tmp);
    }

    for (element_ctr_2=0; element_ctr_2 < descr2->no_of_elements; element_ctr_2++) {
	element = descr2->element+element_ctr_2;
	/* get rid of the insertion tag, in case we
	   haven't shaken if off along the way */
	memcpy (tmp,  element->begin_id, PDB_ATOM_RES_NO_LEN*sizeof(char) );
	element_2_begin[element_ctr_2] = atoi (tmp);
	
	memcpy (tmp,  element->end_id, PDB_ATOM_RES_NO_LEN*sizeof(char) );
	element_2_end[element_ctr_2] = atoi (tmp);
    }
    
    /* translate the beginning and the end of SSEs from pdb_tag to the position on the sequence*/
    element_ctr_1 = 0;
    for (resctr1=0; resctr1<no_res_1; resctr1++) {
	
	memcpy (tmp, protein1->sequence[resctr1].pdb_id, PDB_ATOM_RES_NO_LEN*sizeof(char) );
	num_pdb_id = atoi (tmp);

	/* break if we are outside of the last element */
	if ( num_pdb_id > element_1_end [descr1->no_of_elements-1]) {
	    /*see for example 1k4j - the first residue has the number 399 */
	    continue;
	}
	
	if (  num_pdb_id == element_1_end[element_ctr_1] ) {
	    element_1_end [element_ctr_1] = resctr1;
	    element_ctr_1++;
	}
	
	if ( num_pdb_id == element_1_begin[element_ctr_1] ) {
	    element_1_begin [element_ctr_1] = resctr1;
	}
	/* also if there is some diagreement about where the first element starts: */
	if ( ! resctr1 && num_pdb_id >= element_1_begin[element_ctr_1] ) {
	    element_1_begin [element_ctr_1] = resctr1;
	}
    }	    


    element_ctr_2 = 0;
    for (resctr2=0; resctr2<no_res_2; resctr2++) {
	
	memcpy (tmp, protein2->sequence[resctr2].pdb_id, PDB_ATOM_RES_NO_LEN*sizeof(char) );
	num_pdb_id = atoi (tmp);

	/* break if we are outside of the last element */
	if ( num_pdb_id > element_2_end [descr2->no_of_elements-1]) break;
	
	if (  num_pdb_id==element_2_end  [element_ctr_2] ) {
	    element_2_end [element_ctr_2] = resctr2;
	    element_ctr_2++;
	}
	
	if ( num_pdb_id == element_2_begin[element_ctr_2] ) {
	    element_2_begin [element_ctr_2] = resctr2;
	}
	/* also if there is some disagreement about where the first element starts: */
	if ( ! resctr2 && num_pdb_id >= element_2_begin[element_ctr_2] ) {
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
	    if ( resctr1 >= no_res_1) {
		printf ("%s   %d :  %d  %d \n",descr1->name, protein1->length,  resctr1,  no_res_1);
	    } 
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
    /* for all mapped blocks calculate similarity as exp (-d/5) */
    for (element_ctr_1=0; element_ctr_1 < descr1->no_of_elements; element_ctr_1++) {

	
	element_ctr_2 = map->x2y[element_ctr_1];
	if (element_ctr_2 < 0) continue;

	/* if there is an SSE in the other structure that
	   this one maps to, find the rotated translation vector */
	quat_to_R (map->q, R);
	for (i=0; i<3; i++) {
	    rotated_tr[i] = 0.0;
	    for (j=0; j<3; j++) rotated_tr[i] += R[i][j]*rep1->translation[element_ctr_1][j];
	}

	

	for (resctr1=element_1_begin[element_ctr_1]; resctr1<= element_1_end[element_ctr_1]; resctr1++) {
	    
	    /* find the representative atom (CA) for residue resctr1*/
	    find_Calpha ( protein1, resctr1, ca1 );
	    
	    /* translate to the origin */
	    for (i=0; i<3; i++) ca1[i] -= rep1->origin[i];
	    
	    /* rotate and  translate out*/
	    point_rot_tr (ca1, R, rotated_tr, rotated_ca1);
	    
	    for (resctr2=element_2_begin[element_ctr_2]; resctr2<= element_2_end[element_ctr_2]; resctr2++) {
		/* find the representative atom (CA) for residue resctr2*/
		find_Calpha (protein2, resctr2,  ca2 );
		/* translate to cm & out */
		for (i=0; i<3; i++) ca2[i] += -rep2->origin[i] + rep2->translation[element_ctr_2][i];
		/* finally, find the distance & assign the "similarity" score */
		d = two_point_distance (rotated_ca1, ca2);
		similarity[resctr1][resctr2] = exp (-d/d0);
	    }
	}
    }
    
    
    /* run Smith-Waterman and use the mapped CA to find the transformation        */
    /* I have another copy of SW here (the first one is in struct_map.c)          */
    /* so I wouldn't fumble with parameters - the two should be joined eventually */
    smith_waterman_2 (no_res_1, no_res_2, similarity,
		      map_i2j,  map_j2i, &aln_score);

    map_count = 0;
    for (resctr1=0; resctr1<no_res_1; resctr1++) {
	resctr2 = map_i2j[resctr1];
	if (resctr2 < 0 ) continue;
	map_count ++;
    }
    map_size = map_count;
     

    /*  find the transformation mapping this particular set
	of atoms one onto another */

    map_count = 0;
    for (resctr1=0; resctr1<no_res_1; resctr1++) {
	resctr2 = map_i2j[resctr1];
	if (resctr2 < 0 ) continue;
	
	find_Calpha ( protein1, resctr1, ca1 );
	find_Calpha ( protein2, resctr2, ca2 );

	if (map_count >= no_res_1+no_res_2 ) {
	    fprintf (stderr, "bleep\n"); exit (1);
	}
	for (i=0; i<3; i++) {
	    x[i][map_count] = ca1[i];
	    y[i][map_count] = ca2[i];
	}
	map_count ++;
    }

    map2rotation (x, y, map_size, q, T, &rmsd);
    quat_to_R (q, R);
    
    /* calculate some sort of the match score */
    aln_score = 0.0;
    for (resctr1=0; resctr1<no_res_1; resctr1++) {
	resctr2 = map_i2j[resctr1];
	if (resctr2 < 0 ) continue;
	
	find_Calpha ( protein1, resctr1, ca1 );
	find_Calpha ( protein2, resctr2, ca2 );
	
	/* rotate and  translate*/
	point_rot_tr (ca1, R, T, rotated_ca1);
		
	d = two_point_distance (rotated_ca1, ca2);
	aln_score +=  exp (-d/d0);
    }
    
# if 0
    printf ("\n");
    for ( i =0; i < 3; i++ ) {
	for ( j =0; j < 3; j++ ) {
	    printf ("%8.3lf ",  R[i][j]);
	}
	printf ("%8.3lf \n", T[i]);
    }
    printf ("alignment score: %8.2f,   map size: %d,  rmsd: %8.2f\n\n", aln_score, map_size, rmsd);
# endif
    
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
    last_element_1 = 0;
    last_element_2 = 0;
    for (element_ctr_1=0; element_ctr_1 < descr1->no_of_elements; element_ctr_1++) {

	
	element_ctr_2 = map->x2y[element_ctr_1];
	if (element_ctr_2 <  0) continue;

	last_element_1 = element_ctr_1;
	last_element_2 = element_ctr_2;
	

	/* this is an exact map point  -- these two we won't try aligning to anything else */
	for (resctr1=element_1_begin[element_ctr_1]; resctr1<= element_1_end[element_ctr_1]; resctr1++) {
	    /* find the representative atom (CA) for residue resctr1*/
	    find_Calpha ( protein1, resctr1, ca1 );
	    /* rotate and  translate*/
	    point_rot_tr (ca1, R, T, rotated_ca1);
		
	    for (resctr2=element_2_begin[element_ctr_2]; resctr2<= element_2_end[element_ctr_2]; resctr2++) {
		/* find the representative atom (CA) for residue resctr2*/
		find_Calpha (protein2, resctr2,  ca2 );

		d = two_point_distance (rotated_ca1, ca2);
		similarity[resctr1][resctr2] = exp (-d/d0);
	    }
	}

	/****************************************************/
	/* try to extend the match to the previous loop     */
	if ( element_1_begin[element_ctr_1] == 0 ) continue;
	if ( element_2_begin[element_ctr_2] == 0 ) continue;

	/* the first and the last residues of the preceding loop */
	preceding_loop (element_1_begin, element_1_end, element_ctr_1,
			&first_res_1, &last_res_1);
	
	preceding_loop (element_2_begin, element_2_end, element_ctr_2,
			&first_res_2, &last_res_2);

	
	/******************/
	for (resctr1=first_res_1; resctr1<=last_res_1; resctr1++) {
	    find_Calpha ( protein1, resctr1, ca1 );
	    /* rotate and  translate*/
	    point_rot_tr (ca1, R, T, rotated_ca1);
	    /* the two, however, cannot be matched if they are not
	       of the same type (helix, strand or unassigned)*/
	    for (resctr2=first_res_2; resctr2<=last_res_2; resctr2++) {
		if ( type_1[resctr1] != type_2[resctr2] ) continue;
		
		find_Calpha (protein2, resctr2,  ca2);

		d = two_point_distance (rotated_ca1, ca2);
		if ( d<9.0 ) { 
		    similarity[resctr1][resctr2] = exp (-d/d0);
		}
	    }
	}


	/****************************************************/
	/* now try to extend the match to the following loop */
	if ( element_1_end[element_ctr_1] >= no_res_1 -1) continue;
	if ( element_2_end[element_ctr_2] >= no_res_2 -1) continue;

	/* the first and the last residues of the following loop */
	following_loop (element_1_begin, element_1_end, descr1->no_of_elements, no_res_1,
			element_ctr_1, &first_res_1, &last_res_1);
	
	following_loop (element_2_begin, element_2_end, descr2->no_of_elements, no_res_2,
			element_ctr_2, &first_res_2, &last_res_2);


	
	/******************/
	for (resctr1=first_res_1; resctr1<=last_res_1; resctr1++) {
	    find_Calpha ( protein1, resctr1, ca1 );
	    /* rotate and  translate*/
	    point_rot_tr (ca1, R, T, rotated_ca1);
	    /* the two, however, cannot be matched if they are not
	       of the same type (helix, strand or unassigned)*/
	    for (resctr2=first_res_2; resctr2<=last_res_2; resctr2++) {
		if ( type_1[resctr1] != type_2[resctr2] ) continue;
		
		find_Calpha (protein2, resctr2,  ca2);

		d = two_point_distance (rotated_ca1, ca2);
		if ( d<9.0 ) { 
		    similarity[resctr1][resctr2] = exp (-d/d0);
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
		find_Calpha ( protein1, resctr1, ca1 );
		/* rotate and  translate*/
		point_rot_tr (ca1, R, T, rotated_ca1);
		/* the two, however, cannot be matched if they are not
		   of the same type (helix, strand or unassigned)*/
		for (resctr2=first_res_2; resctr2<=last_res_2; resctr2++) {
		  if ( type_1[resctr1] != type_2[resctr2] ) continue;
		
		  find_Calpha (protein2, resctr2,  ca2);

		  d = two_point_distance (rotated_ca1, ca2);
		  if ( d<9.0 ) { 
		    similarity[resctr1][resctr2] = exp (-d/d0);
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
/* 	prev_2  = map->x2y[element_ctr_1-1]; */
/* 	if (prev_2 < 0) continue; */
/* 	next_2  = map->x2y[element_ctr_1+1]; */
/* 	if (next_2 < 0) continue; */
/* 	if (next_2 - prev_2 != 2 )continue; */
/* 	if (descr1->type[element_ctr_1+1] != descr2->type[prev_2+1] ) continue; */
/* 	printf ("looks like we've skipped %d --> %d \n",  */
/* 		element_ctr_1+1, prev_2+1); */
    }
    /***********************************************/
    /***********************************************/
    /***********************************************/
    /***********************************************/
    
    memset (map_i2j, 0, no_res_1*sizeof(int)); 
    memset (map_j2i, 0, no_res_2*sizeof(int)); 
    smith_waterman_2 (no_res_1, no_res_2, similarity,
		      map_i2j,  map_j2i, &aln_score);


    map_count = 0;
    for (resctr1=0; resctr1<no_res_1; resctr1++) {
	resctr2 = map_i2j[resctr1];
	if (resctr2 < 0 ) continue;
	map_count ++;
    }
    map_size = map_count;
   
    for (i=0; i<3; i++)  memset (x[i], 0, (no_res_1+no_res_2)*sizeof(double) );
    for (i=0; i<3; i++)  memset (y[i], 0, (no_res_1+no_res_2)*sizeof(double) );
 
    map_count = 0;
    for (resctr1=0; resctr1<no_res_1; resctr1++) {
	resctr2 = map_i2j[resctr1];
	if (resctr2 < 0 ) continue;
	
	find_Calpha ( protein1, resctr1, ca1 );
	find_Calpha ( protein2, resctr2, ca2 );

	if (map_count >= no_res_1+no_res_2 ) {
	    fprintf (stderr, "bleep\n"); exit (1);
	}
	for (i=0; i<3; i++) {
	    x[i][map_count] = ca1[i];
	    y[i][map_count] = ca2[i];
	}
	map_count ++;
    }
    
    
    /******************************************/
    /* find the rotation                      */
    map2rotation (x, y, map_size, q, T, &rmsd);
    quat_to_R (q, R);
 
    /* calculate some sort of the match score */
    aln_score = 0.0;
    for (resctr1=0; resctr1<no_res_1; resctr1++) {
	resctr2 = map_i2j[resctr1];
	if (resctr2 < 0 ) continue;
	
	find_Calpha ( protein1, resctr1, ca1 );
	find_Calpha ( protein2, resctr2, ca2 );
	
	/* rotate and  translate*/
	point_rot_tr (ca1, R, T, rotated_ca1);
		
	d = two_point_distance (rotated_ca1, ca2);
	aln_score +=  exp (-d/d0);

    }

    memcpy (&(map->q[0]), &q[0], 4*sizeof(double) );
    memcpy (&(map->T[0]), &T[0], 3*sizeof(double) );
    
    map->x2y_residue_level  = map_i2j;
    map->y2x_residue_level  = map_j2i;
    
    map->x2y_residue_l_size = no_res_1;
    map->y2x_residue_l_size = no_res_2;
    
# if 0
    for ( i =0; i < 3; i++ ) {
	for ( j =0; j < 3; j++ ) {
	    printf ("%8.3lf ",  R[i][j]);
	}
	printf ("%8.3lf \n", T[i]);
    }
    printf ("alignment score extended: %8.2f,   map size: %d,  rmsd: %8.2f\n\n", aln_score, map_size, rmsd);
# endif
    /*************************************************************************/
    map->res_almt_score  = aln_score;
    map->res_almt_length = map_size;
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
    free (type_1);
    free (type_2);

 
     
    return 1;
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
    char ** direction;
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
    int i;
    
    for (i=0; i<res->no_atoms; i++) {
	if ( ! res->atom[i].backbone) continue;
	if ( ! strcmp(res->atom[i].type, "CA")) {
	    ca[0]= res->atom[i].x;
	    ca[1]= res->atom[i].y;
	    ca[2]= res->atom[i].z;
	    break;
	}
    }


    return 0;

}



# if 0
	    /* printf (" res %3d  (%s) belongs to element %2d\n", */
	    /* 	    resctr1, protein1->sequence[resctr1].pdb_id, element_ctr_1); */
# endif



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

# endif

	
# if 0
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
# endif


