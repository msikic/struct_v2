# include "struct.h"

/**************************************************/
/***************************************************/
/**************************************************/
/***************************************************/
/**************************************************/
/***************************************************/
/**************************************************/
/***************************************************/
/**************************************************/
/***************************************************/
int rep_initialize (Representation * rep, Descr * descr  ){

    int N = descr->no_of_elements;
    int i;
    int construct_compact (Representation *rep);
    
    rep->N_full = N;
    
    rep->full_no_of_strands = descr->no_of_strands;
    if ( ! (rep->full  = emalloc(N*sizeof(double*))) ) return 1;
    for (i=0; i<N; i++ ) {
	/* note we are not allocating here - only storing the pointer*/
	rep->full[i] = descr->element[i].p[0];
    }
    if ( ! (rep->cm  = emalloc(N*sizeof(double*))) ) return 1;
    for (i=0; i<N; i++ ) {
	rep->cm[i]   = descr->element[i].cm[0];
    }
   
    if ( ! (rep->translation  = dmatrix(N,3)) ) return 1;
    if ( ! (rep->transl_norm  = emalloc(N*sizeof(double)))) return 1;
    
    rep->full_type = descr->type;
    rep->full_sheet_id  = descr->sheet_id;
    rep->length  = descr->length;

    if ( ! (rep->compact  = dmatrix (N,3) )) return 1;
    if ( ! (rep->compact_type  = emalloc (N*sizeof (int)) )) return 1;
    
	 
    if ( ! (rep->is_rep_by   = emalloc ( N*sizeof (int))) ) return 1;
    /* 1 for the counter, 1 to represent itself: */
    if ( ! (rep->represents  = intmatrix (N, N )) )return 1;
    if ( construct_compact (rep) ) return 1;

    return 0;
 }

/**************************************************/
/***************************************************/
 
 int rep_shutdown  (Representation * rep) {

     free (rep->full);
     free (rep->cm);
    
     free (rep->compact_type);
     free (rep->is_rep_by );
     free (rep->transl_norm);

     free_imatrix (rep->represents);
     free_dmatrix (rep->compact );
     free_dmatrix (rep->translation);
   

     memset (rep, 0, sizeof(Representation) );
     
     return 0;
 }

/******************************************/
int find_cluster_avg (int *cluster, int * strand_id,
		      Representation * rep, int ** nbrs,  double *avg ) {

    int i, j, s;
    int first_cluster_member, cluster_member;
    int cluster_size = cluster[0];
    double norm;

    if (! cluster_size) return 1;
    
    first_cluster_member =  cluster[1]; /* chooses direction */

    for (j=0; j<3; j++) avg[j] = 0.0;
    for ( s=1; s <= cluster_size; s++ ) {
	cluster_member = cluster[s];
	i = strand_id [cluster_member];
	for (j=0; j<3; j++) {
	    avg[j] +=
		nbrs[first_cluster_member][cluster_member]
		*rep->full[i][j];
	}
    }
    norm = 0.0;
    for (j=0; j<3; j++)  {
	avg[j] /=  cluster_size;
	norm   += avg[j]*avg[j];
    }
    norm = sqrt(norm);
    for (j=0; j<3; j++) avg[j] /= norm;


    return 0;

}

/**********************************************************/
/**********************************************************/
/**********************************************************/
/**********************************************************/

int construct_compact (Representation *rep) {

     int i, j, new_ctr = 0;

     if (options.sheet_cosine < 1.0 ) {
	 int count;
	 /* attempt some conservative (=strict) clustering
	    of the strands belonging to the same sheet,
	    according to the "closeness" of their direction vectors*/
	 int  **cluster;
	 int no_of_clusters, c, s;
	 int cluster_member;
	 int find_cluster_avg (int * cluster, int * strand_id,
		      Representation * rep,  int ** nbrs,  double *avg );
	 /* find sheets */
	 int ** element, no_of_elements;   /* strand id's for each sheet */
	 int ** nbrs;  /* neighboring matrix for cluster counting */
	 int  sse_ctr, sse_ctr_2, sse_type, type_ctr[3] = {0};
	 int *cluster_representative, *visited, cluster_id;
	 int *parent_cluster;
	 int total_nbrs, cluster_ctr;
	 double dot, **avg;

	 no_of_elements = rep->N_full;

	 if ( ! (element = intmatrix (3, rep->N_full)))
	     return 1;
	 if ( ! (nbrs    = intmatrix (no_of_elements, no_of_elements) ) )
	     return 1;
	 if ( ! (cluster = intmatrix (no_of_elements+1, no_of_elements+1)))
	     return 1;
	 
	 if ( !( visited  = emalloc ( no_of_elements*sizeof(int))) ) return 1;
	 if ( !( parent_cluster= emalloc (rep->N_full*sizeof(int)))) return 1;
	 if ( !( cluster_representative=emalloc(no_of_elements*sizeof(int))))
	     return 1;
	 
	 if ( ! ( avg = dmatrix ( no_of_elements, 3) ) ) return 1;
	 
	 /*******************************************/
	 /*   extract helices and strands           */
	 /*******************************************/
	 for (sse_ctr=0; sse_ctr < rep->N_full; sse_ctr++ ) {
	     sse_type = rep->full_type[sse_ctr]/2; /* I have the types defined as 1, 2 or 4*/
	     element[sse_type][ type_ctr[sse_type] ] = sse_ctr;
	     type_ctr[sse_type]++;
	     
	     parent_cluster[sse_ctr] = -1;
	 }


	 /***********************************************************/
	 /*   for each type, find nbr matrix and related clusters   */
	 /***********************************************************/

	 cluster_ctr = 0;
	 for ( sse_type = 0; sse_type < 3; sse_type++ ) {

	     /* neighbor matrix */
	     total_nbrs = 0;
	     for (sse_ctr=0; sse_ctr < type_ctr[sse_type]; sse_ctr++) {
		 
		 nbrs [sse_ctr][sse_ctr] = 1;
		 i = element[sse_type][sse_ctr];
		 
		 for (sse_ctr_2=sse_ctr+1; sse_ctr_2< type_ctr[sse_type]; sse_ctr_2++){

		     j = element[sse_type][sse_ctr_2];
		     unnorm_dot (rep->full[i], rep->full[j], &dot);
		     if  (  dot  > options.sheet_cosine ) {
			 nbrs [sse_ctr_2][sse_ctr]= nbrs [sse_ctr][sse_ctr_2] = 1;
			 total_nbrs ++;
		     }
		 }
	     }
	     if (!total_nbrs) continue;
	     
	     /* neighbor matrix */
	     if ( type_ctr[sse_type] == 2 && nbrs [0][1] < 0 ) {
		 continue;
	     } else if (type_ctr[sse_type] == 2 && nbrs [0][1] > 0){
		 no_of_clusters = 1;
		 memset (cluster[0], 0, (no_of_elements+1)*(no_of_elements+1)*sizeof(int) );
		 cluster[1][0] = 2; /* size */
		 cluster[1][1] = 0;
		 cluster[1][2] = 1;
	     } else {
		 /* TODO: does cluster_couter  have a bug for the cluster size == 2? */
		 /* or maybe if there nbrs matrix is all 0 ? */
		 /* maybe it doesn't -   just not sure */ 
		 cluster_counter (type_ctr[sse_type], nbrs, 
				  &no_of_clusters, cluster);
	     }
	     
	     /* for each cluster larger than 1, find a representative */
	     
	     for ( c=1; c <= no_of_elements ; c++) { /* c=0 stores isolated strands */
		 if ( cluster[c][0] < 2 ) continue;
		 
		 for ( s=1; s <= cluster[c][0]; s++ ) {
		     cluster_member = cluster[c][s];
		      /* all sse's in this cluster use the same representative */
		     parent_cluster[cluster_member] = cluster_ctr;
		 }
		 
		 /* find average direction to use as a representative : */
		 find_cluster_avg (cluster[c], element[sse_type],
				   rep, nbrs, avg[cluster_ctr]);
# if 0
		 printf ("type %1d   cluster %2d,  size: %2d\n ",
			sse_type, cluster_ctr, cluster[c][0]);
		 printf ("\tcluster members: " );
		 for ( s=1; s <= cluster[c][0]; s++ ) {
		     cluster_member = cluster[c][s];
		     i = element[sse_type][cluster_member];
		     printf (" %2d:%2d ", cluster_member, i );
		 }
		 printf ("\n");
		 vec_out ( avg[cluster_ctr], 3 , "avg");
		 printf ("\n\n");
# endif
		 cluster_ctr++;
		 
	     } /* end loop over all clusters  */
	     
	 }/* end loop over types (helix or strand) */

	 //exit(0);
	 /*******************************************************/
	 /* assign representatives for each element             */
	 /*******************************************************/
	 new_ctr = 0;
	 cluster_id = -1;
	 memset ( visited, 0, no_of_elements*sizeof(int) );	 
	 
	 for (i=0; i< rep->N_full; i++) {
	     rep->represents[i][0] = 0;
	 }
	 
	 for (i=0; i< rep->N_full; i++) {

	     cluster_id = parent_cluster[i];
	     
	     if (   cluster_id < 0) {
		 
		 memcpy (rep->compact[new_ctr], rep->full[i], 3*sizeof(double));
		 rep->compact_type  [new_ctr] = rep->full_type[i];
		 rep->is_rep_by[i] = new_ctr;
		 new_ctr++;
		 
	     } else if ( ! visited[cluster_id] ) {
		 
		 /* this  cluster we see  for the 1st time */
		 visited[cluster_id] = 1;
		 
		 for (j=0; j<3; j++) rep->compact[new_ctr][j] = avg[cluster_id][j];
		 rep->compact_type  [new_ctr] = rep->full_type[i];

		 rep->is_rep_by[i] = new_ctr;
		 cluster_representative[cluster_id] = rep->is_rep_by[i];
		 new_ctr ++;
		 
	     } else { /* this cluster already seen: keep track of the representative */
		 rep->is_rep_by[i] = cluster_representative[cluster_id];
	     }
	     
	     /* sort the representatives */
	     j = rep->is_rep_by[i];
	     rep->represents[j][0] += 1;
	     count =  rep->represents[j][0];
	     rep->represents[ j ] [count] = i;
	     
	 }

	 rep->N_compact = new_ctr;
	 free_imatrix (element);
	 free_imatrix (nbrs);
	 free_imatrix (cluster);
	 
	 free (visited);
	 free (parent_cluster);
	 free (cluster_representative);
	 
	 free_dmatrix(avg);

     /******************************************/
     } else { /* clustering not called for */
	 for (i=0; i< rep->N_full; i++) {
	     memcpy (rep->compact[i], rep->full[i], 3*sizeof(double));
	     rep->compact_type[i]   = rep->full_type[i];
	     rep->is_rep_by[i] = i;
	     rep->represents[i][0] = 1; /* this is the count */
	     rep->represents[i][1] = i;
	 }
	 rep->N_compact = rep->N_full;
     }

     
# if 0
     int i_rep, ctr;
     printf ("full set size: %d\n",  rep->N_full);
     for (i=0; i< rep->N_full; i++) {
	 printf ("%2d is represented by ", i);
	 j = rep->is_rep_by[i];
	 printf ("%2d ", j);
	 printf ("\n");
     }
     
     for (i_rep=0; i_rep<rep->N_compact; i_rep++ ) {
	 printf ("%2d  represents ", i_rep);
	 for (ctr=1; ctr <= rep->represents[i_rep][0]; ctr++) {
	     i = rep->represents[i_rep][ctr];
	     printf (" * %2d  *", i);
	 }
	 printf ("\n");
     }
     printf (" effective size: %d \n", rep->N_compact);
     exit (1);
# endif
	 
     return 0;
 }


/**********************************************************/
/**********************************************************/
/**********************************************************/
/**********************************************************/

int construct_compact_beta_sheet (Representation *rep) {

     int i, j, new_ctr = 0;
     int j1, j2;

     if (options.sheet_cosine < 1.0 && rep->full_no_of_strands > 1 ) {
	 int count, done;
	 /* attempt some conservative (=strict) clustering
	    of the strands belonging to the same sheet,
	    according to the "closeness" of their direction vectors*/
	 int  **cluster;
	 int no_of_clusters, c, c_array_size;
	 int s1, s2;
	 int cluster_member;
	 int find_cluster_avg (int * cluster, int * strand_id,
		      Representation * rep,  int ** nbrs,  double *avg );
	 /* find sheets */
	 int ** strand_per_sheet;   /* strand id's for each sheet */
	 int ** nbrs;  /* neighboring matrix for cluster counting */
	 int sheet_id, strand_ctr, str_ctr_2;
	 int * sheet_size, no_of_strands;
	 int * parallel, *cluster_representative, *visited, cluster_id;
	 int * parent_cluster;
	 int total_nbrs, cluster_ctr;
	 double dot, **avg;

	 no_of_strands = rep->full_no_of_strands;

	 if ( !(strand_per_sheet=intmatrix (rep->N_full, no_of_strands)))
	     return 1;
	 if ( !( nbrs=intmatrix (no_of_strands, no_of_strands) ) )
	     return 1;
	 if ( ! ( cluster=intmatrix (no_of_strands+1, no_of_strands+1)))
	     return 1;
	 
	 if (!(sheet_size =emalloc ( rep->N_full*sizeof(int) ) ) )
	     return 1;
	 c_array_size =  no_of_strands;
	 if ( !( parallel = emalloc ( no_of_strands*sizeof(int)))) return 1;
	 if ( !( visited  = emalloc ( no_of_strands*sizeof(int))) ) return 1;
	 if ( !( parent_cluster= emalloc (rep->N_full*sizeof(int)))) return 1;
	 if ( !( cluster_representative=emalloc(no_of_strands*sizeof(int))))
	     return 1;
	 
	 if ( ! ( avg = dmatrix ( no_of_strands, 3) ) ) return 1;
	 
	 /*******************************************************/
	 /* organize strands into sheets                        */
	 /*******************************************************/
	
	 for (i=0; i< rep->N_full; i++) {
	     if ( rep->full_type[i] == HELIX ) continue;
	     
	     sheet_id = rep->full_sheet_id[i]-1;
	     strand_ctr = sheet_size[sheet_id];
	     strand_per_sheet[sheet_id][strand_ctr] = i;
	     sheet_size[sheet_id]++;
	     
	     parent_cluster[i] = -1;
	    
	 }
	 /*******************************************************/
	 /* for each sheet find nbr matrix and related clusters */
	 /*******************************************************/

	 cluster_ctr = 0;
	 for ( sheet_id = 0; sheet_id <  rep->N_full; sheet_id ++ ) {

	     if (sheet_size[sheet_id] < 2 ) continue;
	     
	     /* neighbor matrix */
	     total_nbrs = 0;
	     for (strand_ctr=0; strand_ctr < sheet_size[sheet_id]; strand_ctr++) {
		 
		 nbrs [strand_ctr][strand_ctr] = 1;
		 i = strand_per_sheet[sheet_id][strand_ctr];

		 for (str_ctr_2=strand_ctr+1; str_ctr_2<sheet_size[sheet_id]; str_ctr_2++){

		     j = strand_per_sheet[sheet_id][str_ctr_2];
		     unnorm_dot (rep->full[i], rep->full[j], &dot);
		     if  ( (fabs ( dot ) > options.sheet_cosine) ) {
			 nbrs [str_ctr_2][strand_ctr]
			     = nbrs [strand_ctr][str_ctr_2] 
			     = ( dot > 0 ) ? 1 : -1;
			 total_nbrs ++;
		     }
		 }
	     }
	     if (!total_nbrs) continue;
	     
	     /* neighbor matrix */
	     if ( sheet_size[sheet_id] == 2 && nbrs [0][1] < 0 ) {
		 continue;
	     } else if (sheet_size[sheet_id] == 2 && nbrs [0][1] > 0){
		 no_of_clusters = 1;
		 memset (cluster[0], 0, (no_of_strands+1)*(no_of_strands+1)*sizeof(int) );
		 cluster[1][0] = 2; /* size */
		 cluster[1][1] = 0;
		 cluster[1][2] = 1;
	     } else {
		 /* TODO: does cluster_couter  have a bug for the cluster size == 2? */
		 /* or maybe if there nbrs matrix is all 0 ? */
		 /* maybe it doesn't -   just not sure */ 
		 cluster_counter (sheet_size[sheet_id], nbrs, 
				  &no_of_clusters, cluster);
	     }
	     
	     /* for each cluster of strands, larger than 1 if parallel,
		larger than 2 if antiparallel, find a representative */
	     
	     for ( c=1; c <= c_array_size ; c++) { /* c=0 stores isolated strands */
		 if ( cluster[c][0] < 2 ) continue;

		 /* are all the  members of this cluster parallel ?*/
		 /* s1 and s2 run over all strands in cluster c */
		 parallel[cluster_ctr] = 1; /* innocent until proven guilty */
		 done = 0;
		 for ( s1=1; s1 <= cluster[c][0] && !done; s1++ ) {
		     for ( s2=s1+1; s2 <= cluster[c][0] && !done; s2++ ) {
			 if ( nbrs[s1-1][s2-1] < 0 ) {
			     parallel[cluster_ctr] = 0;
			     done = 1;
			 }
		     }
		 }
		 /* if  ! all_parallel && if cluster < 3, continue;
		     this cannot be compactified */
		 if ( (! parallel[cluster_ctr]) && (cluster[c][0] < 3) )  continue;
		 		 
		 
		 for ( s1=1; s1 <= cluster[c][0]; s1++ ) {
		     cluster_member = cluster[c][s1];
		     i = strand_per_sheet[sheet_id][cluster_member];
		     parent_cluster[i] = cluster_ctr; /* all strands in this
						       cluster use the same
						       representative */
		 }
		 /* find average direction to use as a representative : */
		 find_cluster_avg (cluster[c], strand_per_sheet[sheet_id],
		 	   rep, nbrs, avg[cluster_ctr]);

# if 0
		 printf ("sheet: %2d,   sheet size: %2d,  cluster %2d,  size: %2d, parallel: %1d \n ",
			 sheet_id, sheet_size[sheet_id], cluster_ctr,
			 cluster[c][0], parallel[cluster_ctr]);
		 printf ("\tcluster members: " );
		 for ( s1=1; s1 <= cluster[c][0]; s1++ ) {
		     cluster_member = cluster[c][s1];
		     i = strand_per_sheet[sheet_id] [cluster_member];
		     printf (" %2d:%2d ", cluster_member, i );
		 }
		 printf ("\n");
		 vec_out ( avg[cluster_ctr], 3 , "avg");
		 printf ("\n\n");
# endif
		 
		 cluster_ctr++;
		 
	     } /* end loop over all cluters in this sheet */
	     
	 }/* end loop over beta-sheets */

	 
	 /*******************************************************/
	 /* assign representatives for each element             */
	 /*******************************************************/
	 new_ctr = 0;
	 cluster_id = -1;
	 memset ( visited, 0, no_of_strands*sizeof(int) );	 
	 for (i=0; i< rep->N_full; i++) {
	     rep->represents[i][0] = 0;
	 }
	 
	 for (i=0; i< rep->N_full; i++) {

	     if (  rep->full_type[i] != HELIX ) cluster_id = parent_cluster[i];
	     
	     if (  rep->full_type[i] == HELIX || cluster_id < 0) { 
		 /* this is helix, or no descriptor attached*/
		 memcpy (rep->compact[new_ctr], rep->full[i], 3*sizeof(double));
		 rep->compact_type  [new_ctr] = rep->full_type[i];
		 rep->is_rep_by[i] = new_ctr;
		 new_ctr++;
		 
	     } else if ( ! visited[cluster_id] ) {
		 /* this  cluster we see  for the 1st time */
		 visited[cluster_id] = 1;
		 
		 for (j=0; j<3; j++) rep->compact[new_ctr][j] = avg[cluster_id][j];
		 rep->compact_type  [new_ctr] = rep->full_type[i];

		 if ( parallel[cluster_id] ) {
		     rep->is_rep_by[i] = new_ctr;
		 } else {
		     /* signal two representatives with the - sign */
		     rep->is_rep_by[i] = -(new_ctr+1);
		     new_ctr++;
		     for (j=0; j<3; j++) rep->compact[new_ctr][j] = -avg[cluster_id][j];
		     rep->compact_type  [new_ctr] = rep->full_type[i];
		 }
		 cluster_representative[cluster_id] = rep->is_rep_by[i];
		 new_ctr ++;
	     } else { /* this cluster already seen: keep track of the representative */
		 rep->is_rep_by[i] = cluster_representative[cluster_id];
	     }
	     /* sort the representatives */
	     if ( (j=rep->is_rep_by[i]) < 0 ) { j1=-j-1; j2=-j;}
	     else {j1=j2=j;};
	     for (j=j1; j<=j2; j++) {
		 count = ( ++(rep->represents[j][0]) );
		 rep->represents[ j ] [count] = i;
	     } 
	 }
	 rep->N_compact = new_ctr;
	 free_imatrix (strand_per_sheet);
	 free_imatrix (nbrs);
	 free_imatrix (cluster);
	 
	 free (sheet_size);
	 free (parallel);
	 free (visited);
	 free (parent_cluster);
	 free (cluster_representative);
	 
	 free_dmatrix(avg);
	 
    } else { /* clustering not called for */
	 for (i=0; i< rep->N_full; i++) {
	     memcpy (rep->compact[i], rep->full[i], 3*sizeof(double));
	     rep->compact_type[i]   = rep->full_type[i];
	     rep->is_rep_by[i] = i;
	     rep->represents[i][0] = 1; /* this is the count */
	     rep->represents[i][1] = i;
	 }
	 rep->N_compact = rep->N_full;
     }

# if 0
	 int i_rep, ctr;
	 for (i=0; i< rep->N_full; i++) {
	     printf ("%2d is represented by ", i);
	     if ( (j=rep->is_rep_by[i]) < 0 ) { j1=-j-1; j2=-j;}
	     else {j1=j2=j;};
	     for (j=j1; j<=j2; j++) {
		 printf ("%2d ", j);
	     }
	     printf ("\n");
	 }
	 
	 for (i_rep=0; i_rep<rep->N_compact; i_rep++ ) {
	     printf ("%2d  represents ", i_rep);
	     for (ctr=1; ctr <= rep->represents[i_rep][0]; ctr++) {
		 i = rep->represents[i_rep][ctr];
		 printf (" * %2d  *", i);
	     }
	     printf ("\n");
	 }
	 printf (" effective size: %d \n", rep->N_compact);
	 exit (1);
# endif
     
     
     return 0;
 }


