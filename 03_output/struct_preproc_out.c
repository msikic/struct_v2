# include "struct.h"

int output_by_strand (char * name,  char * outname,char chain, Protein * protein)  {

    int helix_ctr, strand_ctr, ctr, i;
    int is_helix, vect_ctr;
    int eff_no_hlx, eff_no_strnd;
    Helix  * helix;
    Strand * strand;
    FILE * fptr;

    if ( ! (fptr= efopen (outname, "w") ) ) return 1;
    
    eff_no_hlx = 0;
    for ( helix_ctr=0; helix_ctr<protein->no_helices; helix_ctr++ ) {
	helix = protein->helix+helix_ctr;
	eff_no_hlx += helix->no_of_vectors;
    }
    eff_no_strnd = 0;
    for ( strand_ctr=0; strand_ctr<protein->no_strands; strand_ctr++ ) {
	strand = protein->strand+strand_ctr;
	eff_no_strnd += strand->no_of_vectors;
	/* eff_no_strnd += strand->no_of_vectors*2; */
    }
    
    fprintf (fptr, "name: %4s", name);
 
    fprintf  (fptr, "\n");
    fprintf (fptr, "number of residues: %d\n", protein->length);
    fprintf (fptr, "number of helices:  %d\n", eff_no_hlx);
    fprintf (fptr, "number of strands:   %d\n",eff_no_strnd );
    
    
    ctr = 0;
    while (ctr < protein->length) {
	is_helix = 0;
	for ( helix_ctr=0; helix_ctr<protein->no_helices; helix_ctr++ ) {
	    helix = protein->helix+helix_ctr;
	    
	    if ( helix->begin<=ctr &&  helix->end>=ctr )  {
		is_helix = 1;
		ctr =  helix->end+1;
		for (vect_ctr=0; vect_ctr < helix->no_of_vectors; vect_ctr++) {
		    fprintf (fptr, "HELIX   %3d  %3d %5d  %5s %5s ",
			     HELIX, 0, helix->length, helix->begin_id, helix->end_id);
		    for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", helix->p[vect_ctr][i]);
		    for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", helix->cm[vect_ctr][i]);
		    fprintf (fptr, "\n");
		}
		break;
	    }
	}
	if (is_helix) continue;
	
	for ( strand_ctr=0; strand_ctr<protein->no_strands; strand_ctr++ ) {
	    strand = protein->strand+strand_ctr;
	    if ( strand->begin<=ctr   &&  strand->end>=ctr )  {
		ctr =  strand->end+1;
		for (vect_ctr=0; vect_ctr < strand->no_of_vectors; vect_ctr++) {
		    fprintf (fptr, "STRAND  %3d  %3d %5d  %5s %5s ",
			     PARALLEL, strand->sheet_number, /* not sure what's in this field */
			     strand->length, strand->begin_id, strand->end_id);
		    for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", strand->p[vect_ctr][i]);
		    for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", strand->cm[vect_ctr][i]);
		    fprintf (fptr, "\n");

	      
/* 		    fprintf (fptr, "STRAND  %3d  %3d %5d  %5s %5s ", */
/* 			     PERP, strand->sheet_number,  */
/* 			     strand->length, strand->begin_id, strand->end_id); */
/* 		    for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", strand->perp[vect_ctr][i]); */
/* 		    for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", strand->foot[vect_ctr][i]); */
/* 		    fprintf (fptr, "\n"); */
		    
		    
		}
		break;
	    }
	}

	ctr++;
    }
    fprintf (fptr, "#\n");

    fclose (fptr);
    return 0;
}

/******************************************************/
/******************************************************/
# if 0
int output_by_sheet (FILE * fptr, char * name, char chain, Protein * protein) {

    int helix_ctr, sheet_ctr, i, ctr;
    int is_helix;
    int type, no_of_parallel;
    Helix  * helix;
    Sheet * sheet;
    
    no_of_parallel = 0;
    for ( sheet_ctr=0; sheet_ctr<protein->no_sheets; sheet_ctr++ ) {
	sheet = protein->sheet+sheet_ctr;
	no_of_parallel +=  sheet->parallel;
    }
    
    printf ("name: %4s", name);
    if ( chain ) printf ("%c", chain);
    printf  ("\n");
    printf ("number of residues: %d\n", protein->length);
    printf ("number of helices:  %d\n", protein->no_helices);
    printf ("number of parallel sheets:   %d\n", no_of_parallel);
    printf ("number of antiparallel sheets:   %d\n",
	    protein->no_sheets - no_of_parallel );
 
    ctr = 0;
    while (ctr < protein->length) {
	is_helix = 0;
	for ( helix_ctr=0; helix_ctr<protein->no_helices; helix_ctr++ ) {
	    helix = protein->helix+helix_ctr;
	    
	    if ( helix->begin<=ctr &&  helix->end>=ctr )  {
		is_helix = 1;
		ctr =  helix->end+1;
		fprintf (fptr, "HELIX %3d %5d  %5s %5s ", ( type=HELIX) , 
			 helix->length, helix->begin_id, helix->end_id);
		for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", helix->p[i][0]);
		for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", helix->cm[i][0]);
		fprintf (fptr, "\n");
		break;
	    }
	}
	if (is_helix) continue;
	


	for ( sheet_ctr=0; sheet_ctr<protein->no_sheets; sheet_ctr++ ) {
	    sheet = protein->sheet+sheet_ctr;
	    if ( sheet->first_strand->begin<=ctr &&  sheet->first_strand->end>=ctr )  {
		ctr =  sheet->first_strand->end+1;
		type = sheet->parallel? PARALLEL : ANTIPARALLEL ;
		fprintf (fptr, "SHEET %3d %5d  %5s %5s ", type,
			 sheet->number_of_strands, sheet->first_strand->begin_id,
			 sheet->first_strand->end_id);
		for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", sheet->p[i]);
		for (i=0; i<3; i++) fprintf (fptr, "%8.2lf", sheet->first_strand->cm[0][i]);
		fprintf (fptr, "\n");
		
		break;
	    }
	}
	ctr++;
    }
    fprintf (fptr, "#\n");

    return 0;
}
# endif
/******************************************************/
/******************************************************/
