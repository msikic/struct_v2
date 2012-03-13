# include "struct.h"

int descr_out (FILE * fptr, Descr * descr) {

    int element_ctr, i;
    Element * element;
    
    fprintf (fptr, "name: %s \n", descr->name);
    fprintf (fptr, "no of elements: %d \n", descr->no_of_elements);
    
    for (element_ctr=0; element_ctr < descr->no_of_elements; element_ctr++) {
	element = descr->element+element_ctr;
	if ( descr->type[element_ctr] == HELIX ) {
	    fprintf (fptr, "HELIX  %3d %5d  %5s %5s ",
		     descr->type[element_ctr],
		     element->length,  element->begin_id, element->end_id);
	    for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", element->p[0][i]);
	    for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", element->cm[0][i]);
	    fprintf (fptr, "\n");
	} else {
	    fprintf (fptr, "STRAND %3d %5d  %5s %5s ",
		     PARALLEL,
		     element->length,  element->begin_id, element->end_id);
	    for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", element->p[0][i]);
	    for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", element->cm[0][i]);
	    fprintf (fptr, "\n");
	}
	
    }
    fprintf (fptr, "#\n"); /* end of record */
    
    return 0;
}


# if 0

/* junkyard */
	    
fprintf (fptr, "STRAND %3d %5d  %5s %5s ",
	 PERP,
	 element->length,  element->begin_id, element->end_id);
for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", element->perp[0][i]);
for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", element->foot[0][i]);
fprintf (fptr, "\n");

# endif
