# include "struct.h"

int input_by_strand (FILE * fptr, Descr * description);
int input_by_sheet (FILE * fptr, Descr * description);


int get_next_descr (FILE * fptr,  Descr * description) {

    int retval;
    
    if ( ! fptr ) {
	fprintf (stderr, "Error in get_next_decr: read on closed fptr (?)\n");
	return -1;
	
    } 
    
    retval = input_by_strand  (fptr,  description);
    return retval;

}
/*************************************************/
int input_by_strand (FILE * fptr, Descr * description) {

    int no_of_elements = 0;
    int no_of_helices = 0, no_of_strands = 0;
    int type, retval;
    int element_ctr = 0; 
    int descr_initialized =0;
    char line[BUFFLEN];
    Element * element;

    if ( feof(fptr) ) return -1;
    
    if (description->element) descr_shutdown(description);
    
    element_ctr = 0;
    memset (line,  0, BUFFLEN);
    
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( line[0] == '#' ) break; /* the end of record */ 
	if ( ! strncmp(line, "name", 4) ) {

	    sscanf (line, "%*s %s", description->name);
	    	    
	} else if ( ! strncmp(line, "number", 6) ) {

	    if ( strstr(line, "residues")) {
		sscanf ( line, "%*s %*s%*s  %d", &(description->no_of_residues));
	    } else if ( strstr(line, "helices")){
		sscanf ( line, "%*s %*s %*s  %d", &no_of_helices);
		no_of_elements += no_of_helices;
		description->no_of_helices   = no_of_helices;
		description->no_of_elements += no_of_helices;
	    } else if (strstr(line, "strands")){
		sscanf ( line, "%*s %*s %*s  %d", &no_of_strands);
		/* if ( ! options.use_perp) no_of_strands/= 2; */
		no_of_elements += no_of_strands;
		description->no_of_strands = no_of_strands;
		description->no_of_elements += no_of_strands;
	    } else {
		fprintf (stderr, "Unrecognized record in the input:\n%s\n",
			 line);
		return 1;
	    } 
		
	} else if (  (!strncmp(line, "HELIX", 5)) ||
		    (!strncmp(line, "STRAND", 6))  ) {
	    
	    if ( ! descr_initialized ) {
		descr_initialized = 1;
		if ( (retval = descr_init(description) ) ) return retval;
	    }

	    
	    /* type check */
	    sscanf ( line, "%*s %d ", &type);
	    
	    /* if (  options.use_perp || */
            /* 		  (!options.use_perp  && type != PERP) ) { */
		
		
	    /* store element info */
	    if ( element_ctr >= description->no_of_elements ) {
		fprintf (stderr, "number of elements (%d) in %s\n",
			 element_ctr+1, description->name) ;
		fprintf (stderr, "does not match the number declared (%d)\n",
			 description->no_of_elements);
		return 1;
	    }
	    element = description->element + element_ctr;
	    sscanf ( line, "%*s %d  %d   %d %s %s %lf %lf %lf %lf %lf %lf \n",
		     &(description->type[element_ctr]),
		     &(description->sheet_id[element_ctr]),
		     &(element->length), element->begin_id, element->end_id,
		     element->p[0], element->p[0]+1, element->p[0]+2, 
		     element->cm[0], element->cm[0]+1, element->cm[0]+2);
	    description->length[element_ctr] = element->length;

	    element_ctr++;
	    /*  } */
	    
	    
	} 
	memset (line,  0, BUFFLEN);
    }
    if ( !no_of_elements && feof(fptr) ) return -1;
    if ( description->no_of_elements !=  element_ctr ) {
	fprintf (stderr, "number of elements (%d) in %s\n",
		 element_ctr, description->name) ;
	fprintf (stderr, "does not match the number declared (%d)\n",
		 description->no_of_elements);
	return 1;
    }
    
    /* make sure that p is normalized */
    
    for ( element_ctr=description->no_of_elements-1; element_ctr >=0; element_ctr-- ) {
	double norm, aux;
	int i;
	norm = 0;
	for (i=0; i<3; i++) {
	    aux = description->element[element_ctr].p[0][i];
	    norm += aux*aux;
	}
	norm = sqrt (norm);
	for (i=0; i<3; i++) {
	    description->element[element_ctr].p[0][i] /= norm;
	}
    }
    
    return 0;
}


# if 0
/*************************************************/
int input_by_sheet (FILE * fptr, Descr * description) {

    int no_of_helices, no_of_sheets;
    int h, s, retval;
    int element_ctr = 0; 
    int descr_initialized =0;
    char line[BUFFLEN];
    Element *element, * prev_element;
    
    if ( feof(fptr) ) return -1;
    
    if (description->element) descr_shutdown(description);
    
    element_ctr = 0;
    memset (line,  0, BUFFLEN);
    
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	if ( line[0] == '#' ) break; /* the end of record */ 
	if ( ! strncmp(line, "name", 4) ) {
	    
	    sscanf (line, "%*s %s", description->name);
	    	    
	    
	} else if ( ! strncmp(line, "number", 6) ) {
	    
	    if ( strstr(line, "residues")) {
		sscanf ( line, "%*s %*s %*s %d", &(description->no_of_residues));
	    } else if ( strstr (line, "helices")){
		sscanf ( line, "%*s %*s %*s %d", &no_of_helices);
		description->no_of_elements += no_of_helices;
		description->no_of_helices = no_of_helices;
	    } else if ( strstr (line, "paral")){
		sscanf (line, "%*s %*s%*s %*s %d", &no_of_sheets);
		description->no_of_elements += no_of_sheets;
		description->no_of_p_sheets = no_of_sheets;
	    } else if ( strstr(line, "antip")){
		sscanf ( line,"%*s %*s%*s %*s %d", &no_of_sheets);
		description->no_of_elements += 2*no_of_sheets; /* put both directions
						     in the set of vectors */
		description->no_of_a_sheets = no_of_sheets;
	    } else {
		fprintf (stderr, "Unrecognized record in the input:\n%s\n",
			 line);
		return 1;
	    }
		
	} else if ( ! (h=strncmp(line, "HELIX", 5)) ||
		    ! (s=strncmp(line, "SHEET", 5))  ) {

	    if ( ! descr_initialized ) {
		descr_initialized = 1;
		if ( (retval = descr_init(description) ) ) return retval;
	    }
	    /* store element info */
	    if ( element_ctr >= description->no_of_elements ) {
		fprintf (stderr, "number of elements (%d)\n", element_ctr) ;
		fprintf (stderr, "does not match the number declared (%d)\n",
			 description->no_of_elements);
		return 1;
 	    }
	    element = description->element + element_ctr;
	    sscanf ( line, "%*s %d  %d %s %s %lf %lf %lf %lf %lf %lf \n",
		     &(description->type[element_ctr]),
		     &(element->length), element->begin_id, element->end_id,
		     element->p[0], element->p[0]+1, element->p[0]+2, 
		     element->cm[0], element->cm[0]+1, element->cm[0]+2);
	    description->length[element_ctr] = element->length;
	    
	    element_ctr++;

	    if ( description->type[element_ctr-1] == ANTIPARALLEL ) {
		prev_element = element;
		element = description->element + element_ctr;
		description->type[element_ctr] = description->type[element_ctr-1];
		element->length   =  prev_element->length;
		memcpy (element->begin_id, prev_element->begin_id,
			(PDB_SHEET_BEGIN_LEN+2)*sizeof(char) );
		memcpy (element->end_id, prev_element->end_id,
			(PDB_SHEET_END_LEN+2)*sizeof(char) );
		element->p[0][0]  = -prev_element->p[0][0]; 
		element->p[0][1]  = -prev_element->p[0][1]; 
		element->p[0][2]  = -prev_element->p[0][2]; 
		memcpy (element->cm[0],prev_element->cm[0], 3*sizeof(double) ) ;
		description->length[element_ctr] = element->length;
	    
		element_ctr++;

	    }
	    
	} 
	memset (line,  0, BUFFLEN);
    }
    if ( !description->no_of_elements && feof(fptr) ) return -1;
    if ( description->no_of_elements != element_ctr ) {
	fprintf (stderr, "number of elements (%d)\n", element_ctr) ;
	fprintf (stderr, "does not match the number declared (%d)\n",
		 description->no_of_elements);
	return 1;
    }

    /* make sure that p is normalized */
    
    for ( element_ctr=description->no_of_elements-1; element_ctr >=0; element_ctr-- ) {
	double norm, aux;
	int i;
	norm = 0;
	for (i=0; i<3; i++) {
	    aux = description->element[element_ctr].p[0][i];
	    norm += aux*aux;
	}
	norm = sqrt (norm);
	for (i=0; i<3; i++) {
	    description->element[element_ctr].p[0][i] /= norm;
	}
    }
    
    return 0;
}

# endif

int free_element (Descr * descr ) {

    int ctr;
    Element * element = descr->element;

    for (ctr=0; ctr < descr->no_of_elements;  ctr++ ) {
	free_dmatrix ( (element+ctr)->p );
	free_dmatrix ( (element+ctr)->cm );
    }
    free (element);
    descr-> no_of_elements = 0;

    
    return 0;
}


int alloc_element (Descr * descr) {

    int ctr, no_of_elements = descr->no_of_elements;
    Element * element;
    
    element = emalloc ( no_of_elements*sizeof(Element) );
    if ( ! element) return 1;
    
    for (ctr=0; ctr < no_of_elements; ctr++ ) {
	(element+ctr)->p = dmatrix(1,3);
	if ( ! ((element+ctr)->p ) ) return 1;
	(element+ctr)->cm = dmatrix(1,3);
	if ( ! ((element+ctr)->cm) ) return 1;
    }

    descr->element = element;
    descr->no_of_elements = no_of_elements;
    
    return 0;
}
/****************************************************************/

int descr_init ( Descr * description ) {
		
    int retval, no_of_elements = description->no_of_elements;
    
    retval = alloc_element (description);
    if (retval) return retval;
		
    description->type = emalloc (no_of_elements*sizeof(int));
    if (!description->type) return 1;

    description->length = emalloc (no_of_elements*sizeof(int));
    if (!description->length) return 1;

    description->sheet_id = emalloc (no_of_elements*sizeof(int));
    if (!description->sheet_id) return 1;


    return 0;
    
}


int descr_shutdown ( Descr * description ) {
    
    int retval;
    
    if ( description->element) {
	retval = free_element ( description );
	if (retval) return retval;
    }

    if ( description->type) free ( description->type );
    if ( description->length) free ( description->length );
    if ( description->sheet_id) free ( description->sheet_id );

    memset (description, 0, sizeof(Descr) );

    return 0;
}

