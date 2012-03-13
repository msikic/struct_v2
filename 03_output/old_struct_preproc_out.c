# include "struct.h"

int output (FILE * fptr, Protein * protein) {

    int helix_ctr, strand_ctr, i, ctr;
    int is_helix, is_strand;
    Helix  * helix;
    Strand * strand;
    
    ctr = 0;
    while (ctr < protein->length) {
	is_helix = 0;
	for ( helix_ctr=0; helix_ctr<protein->no_helices; helix_ctr++ ) {
	    helix = protein->helix+helix_ctr;
	    if ( helix->begin<=ctr &&  helix->end>=ctr )  {
		is_helix = 1;
		ctr =  helix->end+1;
		fprintf (fptr, "HELIX  %5d  %5s %5s ",
			 helix->length, helix->begin_id, helix->end_id);
		for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", helix->p[i]);
		for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", helix->cm[i]);
		fprintf (fptr, "\n");
		break;
	    }
	}
	if ( is_helix) continue;

	is_strand = 0;
	for ( strand_ctr=0; strand_ctr<protein->no_strands; strand_ctr++ ) {
	    strand = protein->strand+strand_ctr;
	    if ( strand->begin<=ctr &&  strand->end>=ctr )  {
		is_strand = 0;
		ctr =  strand->end+1;
		fprintf (fptr, "STRAND %5d  %5s %5s ",
			 strand->length, strand->begin_id, strand->end_id);
		for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", strand->p[0][i]);
		for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", strand->cm[0][i]);
		fprintf (fptr, "\n");
		break;
	    }
	}
	if ( is_strand) continue;
	ctr++;
    }
    fprintf (fptr, "#\n");

    return 0;
}

/******************************************************/
/******************************************************/

int output_by_type (FILE * fptr, Protein * protein) {

    int helix_ctr, strand_ctr, i;
    Helix  * helix;
    Strand * strand;
    
    
    for ( helix_ctr=0; helix_ctr<protein->no_helices; helix_ctr++ ) {
	helix = protein->helix+helix_ctr;
	fprintf (fptr, "HELIX  %5d  %5s %5s ", helix->length, helix->begin_id, helix->end_id);
	for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", helix->p[i]);
	for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", helix->cm[i]);
	fprintf (fptr, "\n");
    }
    
    for ( strand_ctr=0; strand_ctr<protein->no_strands; strand_ctr++ ) {
	strand = protein->strand+strand_ctr;
	fprintf (fptr, "STRAND %5d  %5s %5s ",
		 strand->length, strand->begin_id, strand->end_id);
	for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", strand->p[0][i]);
	for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", strand->cm[0][i]);
	fprintf (fptr, "\n");
    }
    fprintf (fptr, "#\n");

    return 0;
}
