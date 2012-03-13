# include "struct.h"


int protein_shutdown (Protein * protein) {

    protein->length = 0;
    if (protein->sequence) free(protein->sequence);
    
    if (protein->no_helices) {
	int helix_ctr;
	for (helix_ctr=0; helix_ctr<protein->no_helices; helix_ctr++) {
	    if (protein->helix[helix_ctr].p ) free_dmatrix (protein->helix[helix_ctr].p );
	    if (protein->helix[helix_ctr].cm) free_dmatrix (protein->helix[helix_ctr].cm );
	}
	protein->no_helices = 0;
    }
    if (protein->helix)free (protein->helix);
    
    if (protein->no_strands) {
	int strand_ctr;
	for (strand_ctr=0; strand_ctr<protein->no_strands; strand_ctr++) {
	    if (protein->strand[strand_ctr].p) free_dmatrix (protein->strand[strand_ctr].p);
	    if (protein->strand[strand_ctr].cm) free_dmatrix (protein->strand[strand_ctr].cm);
	}
	protein->no_strands = 0;
    } 
    if (protein->strand) free (protein->strand);
    protein->no_sheets= 0;
    if (protein->sheet) free (protein->sheet);
    if (protein->sse_sequence) free (protein->sse_sequence);

    
    return 0;
}
