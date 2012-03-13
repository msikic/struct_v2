# include "struct.h"

/* output for further processing by another program*/
int rec_map_out_for_postproc (Map * map, int map_ctr, 
			      Representation *X_rep, Representation *Y_rep, int depth) {
    
    static FILE * fptr = NULL;
    int print_map_for_postproc (FILE *fptr,  Map * map,
				Representation *X_rep, Representation *Y_rep);
    int print_transformation (FILE* fptr, Map * map, int depth);

    if ( ! depth) {
	char outname[BUFFLEN] = {'\0'};
	/* as a temp measure: all maps have the  same name */
	sprintf (outname, "%s.for_postp",
		 options.outname);
	fptr  = efopen (outname, "w");
	if ( !fptr) exit (1);

	//fprintf (fptr, "map");
	if (map[map_ctr].submatch_best && map[map_ctr].submatch_best[0]) {
	    //fprintf (fptr, " with submap");
	}
	//fprintf (fptr, "\n");
    } else {
	//fprintf (fptr, "submap");
    }
    //print_map_for_postproc (fptr,  map+map_ctr, X_rep, Y_rep);
    print_transformation   (fptr,  map+map_ctr, depth);

			  
    if (map[map_ctr].submatch_best && map[map_ctr].submatch_best[0] > -1 ) {
	
	int sub_best_ctr, sub_map_ctr;
	sub_best_ctr = 0;
	while (  (sub_map_ctr = map[map_ctr].submatch_best[sub_best_ctr] ) > -1 ) {
	    if ( depth < 1 )  {/* TODO this is only two submaps*/
		rec_map_out_for_postproc ( map, sub_map_ctr, X_rep, Y_rep, depth+1); 
	    }
	    sub_best_ctr ++;
	    if (sub_best_ctr==1) break; /*TODO check this cutoff out --output only the best */
	} 
    }
 
    if ( ! depth) {
	fclose (fptr);
    }
    return 0;
    
}
/************************************************/
/************************************************/
/************************************************/

/* output for a human reader (in verbose option)*/

int recursive_map_out (Map * map, int map_ctr, 
		       Descr * descr1, Descr * descr2,
		       Protein *protein1, Protein *protein2,
		       int depth) {
    
    int ctr;
 
    if ( ! depth) {
	printf ("####################################\n");
	printf ("**  %s   %s ", descr1->name, descr2->name);
	if (map[map_ctr].submatch_best && map[map_ctr].submatch_best[0] > -1) {
	    printf ("  with submap (%d)", map[map_ctr].submatch_best[0]+1);
	}
	printf ("\n");
    } else {
	for (ctr=0; ctr < depth; ctr++) printf ("\t");
    }
    printf ("map: %3d      geometric z_score: %6.3lf \n", map_ctr+1,
	    map[map_ctr].z_score);
    print_map (stdout, map+map_ctr, descr1, descr2,   protein1, protein2,  depth);
    
    if (map[map_ctr].submatch_best && map[map_ctr].submatch_best[0] > -1 ) {
	
	int sub_best_ctr, sub_map_ctr;
	sub_best_ctr = 0;
	while ( sub_best_ctr < options.number_maps_cpl &&
		(sub_map_ctr = map[map_ctr].submatch_best[sub_best_ctr] ) > -1 ) {
	    for (ctr=0; ctr < depth; ctr++) printf ("\t");
	    printf ("map %3d   submatch  complementarity z-score:  %6.3lf \n",
		    map_ctr+1, map[map_ctr].compl_z_score); 
	    recursive_map_out (map, sub_map_ctr, descr1, descr2,
			       protein1, protein2, depth+1); 
	    sub_best_ctr ++;
	    if (sub_best_ctr==1) break; /*output only the best */
	} 
    }
    return 0; 
}
/************************************************/
/************************************************/
/************************************************/
int print_transformation (FILE* fptr, Map * map, int depth){

     
    double **R;
    int i, j, ctr;
    if ( ! (R=dmatrix(3,3) ) ) return 1; /* compiler is bugging me otherwise */
    quat_to_R (map->q, R);
    //fprintf (fptr,"transformation\n");
    for (i=0; i<3; i++) {
	if ( depth) for ( ctr=0; ctr < depth; ctr++) fprintf (fptr, "\t");
	for (j=0; j<3; j++) fprintf (fptr,"%8.3lf ", R[i][j]);
	fprintf (fptr,"%8.3lf  \n", map->T[i]);
    }

    free_dmatrix (R);
   
  
    return 0;
    
}

/************************************************/
/************************************************/
/************************************************/
int print_transformation_0 (FILE* fptr, Map * map, int depth){

    double **R;
    int i, j, ctr;
    if ( ! (R=dmatrix(3,3) ) ) return 1; /* compiler is bugging me otherwise */
    quat_to_R (map->q, R);
    fprintf (fptr,"rotation\n");
    for (i=0; i<3; i++) {
	if ( depth) for ( ctr=0; ctr < depth; ctr++) fprintf (fptr, "\t");
	for (j=0; j<3; j++) fprintf (fptr,"%8.3lf ", R[i][j]);
	fprintf (fptr,"\n");
    }

    free_dmatrix (R);
   
    return 0;
    
}

/************************************************/
/************************************************/
/************************************************/
int print_map_for_postproc (FILE *fptr, Map * map,
			    Representation *X_rep,  Representation *Y_rep) {
    
    int i, j, k;

    fprintf ( fptr, "translations\n");
    for (i=0; i<map->x2y_size; i++) {
	    
	j = map->x2y[i];
	if ( j < 0 ) continue;
	fprintf ( fptr, " %3d ", i+1);
	for (k=0; k<3; k++) {
	    fprintf ( fptr," %8.3lf ",  X_rep->translation[i][k]);
	}
	fprintf ( fptr, " %3d ", j+1);
	for (k=0; k<3; k++) {
	    fprintf ( fptr," %8.3lf ",  Y_rep->translation[j][k]);
	}
	fprintf ( fptr, "\n");
    }

    fprintf ( fptr, "cms\n");
    for (i=0; i<map->x2y_size; i++) {
	    
	j = map->x2y[i];
	if ( j < 0 ) continue;
	fprintf ( fptr, " %3d ", i+1);
	for (k=0; k<3; k++) {
	    fprintf ( fptr," %8.3lf ",  X_rep->cm[i][k]);
	}
	fprintf ( fptr, " %3d ", j+1);
	for (k=0; k<3; k++) {
	    fprintf ( fptr," %8.3lf ",  Y_rep->cm[j][k]);
	}
	fprintf ( fptr, "\n");
    }
    fprintf (fptr,"origin1 ");
    for (k=0; k<3; k++) {
	fprintf (fptr,"%8.3lf ", X_rep->origin[k]);
    }
    fprintf ( fptr, "\n");
    fprintf (fptr,"origin2 ");
    for (k=0; k<3; k++) {
	fprintf (fptr,"%8.3lf ", Y_rep->origin[k]);
    }
    fprintf ( fptr, "\n");
    
    return 0;
}
/************************************************/
/************************************************/
/************************************************/
int print_map (FILE *fptr, Map * map, Descr * descr1, Descr * descr2,
	       Protein * protein1, Protein *protein2, int tab) {
    
    int i, j, index_x, index_y;
    char type;
    
    /* header */
    if (options.print_header) {
	if (tab) fprintf ( fptr, "\t");
	fprintf ( fptr, " %9s   %4s %6s   %8s  %9s  %9s\n", "SSE map", "type",
		  " cosine", "Gauss wt.", "res range", "res range");
	if ( descr1 && descr2) {
	    if (tab) fprintf ( fptr, "\t");
	    fprintf ( fptr, "%5s-->%5s   %4s   %6s %8s in %5s   in %5s\n", descr1->name, descr2->name,
		      "", "", "",  descr1->name, descr2->name);
	}
    }
    /* the actual map */
    for (i=0; i<map->x2y_size; i++) {
	    
	index_x =  i;
	index_y =  map->x2y[i];
	if ( index_y <  0 ) continue;

	type = ' ';
	if ( descr1 ) {
	    if (descr1->type[index_x] == HELIX) {
		type = 'H';
	    } else {
		type = 'S';
	    }
	}
	
	if (tab) fprintf ( fptr, "\t");
	fprintf ( fptr, " %2d --> %2d  %4c  %6.2lf   %8.1le ",
		  index_x+1, index_y+1,  type, map->cosine[index_x][index_y],
		map->image[index_x][index_y]);
	if ( descr1 && descr2)  
	    fprintf ( fptr,  " %4s -%4s  %4s -%4s",
		     descr1->element[index_x].begin_id,
		     descr1->element[index_x].end_id,
		     descr2->element[index_y].begin_id,
		     descr2->element[index_y].end_id);
	fprintf ( fptr, "\n");
    }
    if (tab) fprintf ( fptr, "\t");
    fprintf ( fptr, "total assigned score %8.4lf \n", map->assigned_score);
    //print_transformation   (fptr,  map, 0);

    /* output the alignment */
    if (options.postprocess) {
	
	int pos, last_pos, chunk;
	int last_mapped_x, last_mapped_y;
	int first_mapped_x, first_mapped_y;
	int *aligned_seq_x, *aligned_seq_y;
	int chunk_size = 50;
	char formatted_almt_line[BUFFLEN]; /* BUFFLEN is 150 */
	int format_almt_line (Protein * protein, int * aligned_residues,
			      int first_pos, int last_pos, int chunk_size, char ret_str[]);
	
	if ( !protein1 || ! protein2 ) {
	    fprintf (stderr, "From print_map(): postprocessing without structure (?).\n");
	    return 1;
	}
	if ( !(aligned_seq_x= emalloc ( (protein1->length + protein2->length)*sizeof(int)))) return 1;
	if ( !(aligned_seq_y= emalloc ( (protein1->length + protein2->length)*sizeof(int)))) return 1;
	

	
	last_mapped_x = map->x2y_residue_l_size-1;
	while (last_mapped_x>=0 && map->x2y_residue_level[last_mapped_x] < 0 ) last_mapped_x--;
	last_mapped_x++;

	last_mapped_y = map->y2x_residue_l_size-1;
	while (last_mapped_y>=0 && map->y2x_residue_level[last_mapped_y] < 0 ) last_mapped_y--;
	last_mapped_y++;
    
# if 0 
	for (i=0; i<map->x2y_residue_l_size; i++) {
	    j = map->x2y_residue_level[i];
	    /* looks like I've switched the two around someplace */
	    printf ("**  %2d  %3s %2c ", i, protein2->sequence[i].pdb_id, 
		    protein2->sequence[i].res_type_short);
	    if ( j >= 0 ) {
		printf ("**  %2d  %3s  %2c ", j, protein2->sequence[j].pdb_id,
			protein1->sequence[j].res_type_short);
	    }

	    printf ("\n");
	}

	for (j=0; j<map->y2x_residue_l_size; j++) {
	    i = map->y2x_residue_level[j];
	    /* looks like I've switched the two around someplace */
	    printf ("**  %2d  %3s %2c ", j, protein2->sequence[j].pdb_id, 
		    protein2->sequence[j].res_type_short);
	    if ( i >= 0 ) {
		printf ("**  %2d  %3s  %2c ", i, protein2->sequence[i].pdb_id,
			protein1->sequence[i].res_type_short);
	    }

	    printf ("\n");
	}
	exit (1);
# endif
    
	
	first_mapped_x = -1;
	do {
	    first_mapped_x ++;
	} while (first_mapped_x <map->x2y_residue_l_size
		 && map->x2y_residue_level[first_mapped_x] < 0);
	
	last_mapped_x = map->x2y_residue_l_size;
	do {
	    last_mapped_x--;
	} while (last_mapped_x>=0 && map->x2y_residue_level[last_mapped_x] < 0 ); 



	first_mapped_y = -1;
	do {
	    first_mapped_y ++;
	} while (first_mapped_y <map->y2x_residue_l_size
		 && map->y2x_residue_level[first_mapped_y] < 0);
	
	last_mapped_y = map->y2x_residue_l_size;
	do {
	    last_mapped_y--;
	} while (last_mapped_y>=0 && map->y2x_residue_level[last_mapped_y] < 0 ); 




# if 0	    
	printf ("first mapped x: %d\n", first_mapped_x);
	printf ("first mapped y: %d\n", first_mapped_y);
	printf ("last  mapped x: %d\n", last_mapped_x);
	printf ("last  mapped y: %d\n", last_mapped_y);
# endif	
	
	i=first_mapped_x; j=first_mapped_y; pos = 0;
	while (i<map->x2y_residue_l_size ||  j<map->y2x_residue_l_size) {
	    
	    if ( i == map->x2y_residue_l_size ) {
		aligned_seq_x[pos] = -1;
		aligned_seq_y[pos] =  j;
		j++;
		pos ++;
		
	    } else if (j == map->y2x_residue_l_size) {
		aligned_seq_x[pos] =  i;
		aligned_seq_y[pos] = -1;
		i++;
		pos ++;
		
	    } else if ( map->x2y_residue_level[i] >= 0  && map->y2x_residue_level[j] >= 0 ) {
		aligned_seq_x[pos] =  i;
		i++;
		aligned_seq_y[pos] =  j;
		j++;
		pos ++;
		
	    } else if ( map->x2y_residue_level[i] < 0 ) {
		aligned_seq_x[pos] =  i;
		i++;
		aligned_seq_y[pos] = -1;
		pos ++;
		
		    
	    } else {
		aligned_seq_x[pos] = -1;
		aligned_seq_y[pos] =  j;
		j++;
		pos ++;
		
	    }
	    if ( i > last_mapped_x || j > last_mapped_y ) break;
	}

	last_pos = pos;
	fprintf ( fptr, "\n");
	if (tab) fprintf ( fptr, "\t");
	printf ("alignment score A: %8.2f\nalignment size (number of residues): %d\n",
		map->res_almt_score,  map->res_almt_length);
	printf ("residue rmsd (\\AA): %8.2f\n\n",  map->res_rmsd);
	fprintf ( fptr, "alignment:\n\n");

	pos = 0;
	for ( chunk=0; chunk <= last_pos/chunk_size; chunk++) {


	    format_almt_line (protein1, aligned_seq_y, pos, last_pos, chunk_size, formatted_almt_line);
	    fprintf (fptr,  "%s\n", formatted_almt_line);
	    format_almt_line (protein2, aligned_seq_x, pos, last_pos, chunk_size, formatted_almt_line);
	    fprintf (fptr,  "%s\n", formatted_almt_line);
	    fprintf (fptr, "\n\n");
	    
	    pos += chunk_size;
	}

	/* print also to a file, for the visualization */
	/* TODO: this should be cleaned up, this is not the
	   place to do this output */
	{
	    FILE *fptr2 = stdout;
	    double **R;
	    int resctr_2, resctr_1;
	    double ca1[3], ca2[3], rotated_ca2[3], d;
	    /* open file */
	    char outname[BUFFLEN] = {'\0'};
	    /* as a temp measure: all maps have the  same name */
	    sprintf (outname, "%s.alignment",
		 options.outname);
	    fptr2  = efopen (outname, "w");
	    if ( !fptr2) exit (1);
	    /* the plain old alignment */
	    fprintf ( fptr2, "%% alignment:\n%%\n");

	    pos = 0;
	    for ( chunk=0; chunk <= last_pos/chunk_size; chunk++) {

		 format_almt_line (protein1, aligned_seq_y, pos, last_pos, chunk_size, formatted_almt_line);
		 fprintf (fptr2,  "%% %s\n", formatted_almt_line);
		 format_almt_line (protein2, aligned_seq_x, pos, last_pos, chunk_size, formatted_almt_line);
		 fprintf (fptr2,  "%% %s\n", formatted_almt_line);
		 fprintf (fptr2, "%%\n%%\n");
	    
		 pos += chunk_size;
	    }

	    /* find explicit rotation */
	    if ( ! (R=dmatrix(3,3) ) ) return 1; /* compiler is bugging me otherwise */
 
	    quat_to_R (map->q, R);
	    fprintf ( fptr2, "%% distances of the aligned Ca's:\n%%\n" );
	    fprintf (fptr2, "%% pdb_id_1  type_1    distance  pdb_id_2  type_2\n"); 
	    for (pos=0; pos < last_pos; pos++) {
		 resctr_2 = aligned_seq_x[pos];
		 resctr_1 = aligned_seq_y[pos];

		 if ( resctr_1 < 0) {
		      fprintf (fptr2, "%7c  %7c ",  '-', '-');
		 } else {
		      fprintf (fptr2, "%7s  %7c ",
			      protein1->sequence[resctr_1].pdb_id, 
			      protein1->sequence[resctr_1].res_type_short);
		 }
		 if ( resctr_1 < 0 || resctr_2 < 0 ) {
		      fprintf (fptr2, " %10c ", '-');
		 } else {
		      /** evaluate the distance */
		      find_Calpha ( protein1, resctr_1, ca1 );
		      find_Calpha ( protein2, resctr_2, ca2 );
		      point_rot_tr (ca2, R, map->T, rotated_ca2);
		
		      d = two_point_distance (rotated_ca2, ca1);
		      fprintf (fptr2, " %10.2lf ", d);
		     
		 }
		 if ( resctr_2 < 0) {
		      fprintf (fptr2, "%7c  %7c ",  '-', '-');
		 } else {
		      fprintf (fptr2, "%7s  %7c ",
			      protein2->sequence[resctr_2].pdb_id, 
			      protein2->sequence[resctr_2].res_type_short);
		 }
		 fprintf (fptr2, "\n");
	    }
	    fclose (fptr2);
	    free_dmatrix (R);

	}
	free (aligned_seq_x);
	free (aligned_seq_y);
    }
   
    
    return 0;
}



/*********************************************/

int format_almt_line (Protein * protein, int * aligned_residues,
		      int first_pos, int last_pos, int chunk_size, char ret_str[]) {

    int pos, i, resctr;

    memset (ret_str, 0 , BUFFLEN*sizeof(char));

    
    /* look for the first label in this line*/
    pos = first_pos;
    resctr = -1;
    for (i=0; i< chunk_size && pos < last_pos; i++) {
	resctr = aligned_residues[pos];
	if ( resctr > -1) break;
	pos ++;
    }
    if (  resctr > -1 ) {
	sprintf (ret_str, " %3s ", protein->sequence[resctr].pdb_id);
    } else {
	sprintf (ret_str, " %3s ", "");
    }
    
    pos = first_pos;	    
    /* the next  chunk_size cahracters in the  protein */
    for (i=0; i< chunk_size && pos < last_pos; i++) {
	resctr = aligned_residues[pos];
	if (  resctr > -1 ) {
	    sprintf (ret_str,  "%s%c", ret_str, protein->sequence[resctr].res_type_short);
	} else {
	    sprintf (ret_str,  "%s%c", ret_str, '.');
	}
	pos++;

    }


    /* look for the last label*/
    pos -- ;
    resctr = -1;
    for (i=0; i< chunk_size && pos > 0; i++) {
	resctr = aligned_residues[pos];
	if ( resctr > -1) break;
	pos --;
    }
    if (  resctr > -1 ) {
	sprintf (ret_str, "%s %3s ", ret_str, protein->sequence[resctr].pdb_id);
    }
	    
 
    return 0;

}
