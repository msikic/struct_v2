# include "struct.h"

/* how many vars of each type need to be read in: */
# define DOUBLES 17
# define INTS    6
# define STRINGS 4
# define CHARACTERS 2
    
int read_cmd_file (char *filename) {
    FILE * fptr, *log = stdout; 
    char line[LONGSTRING];
    char token[MAX_TOK][MEDSTRING] = {{'\0'}};
    char comment_char;
    int  max_token;
    int line_ctr, retval;
    int ctr, token_assigned;
    int errmsg ( FILE *fptr, int line_ctr, char line[LONGSTRING],
		 char * fmt, char * warnstr) ;

    double * ptrs_to_doubles[DOUBLES] = {
	&(options.sheet_cosine),
	&(options.tol), &(options.z_max_store), &(options.z_max_out),
    	&(options.F_guess_max), &(options.F_eff_max),
	&(options.z_max_corr),  &(options.z_min_compl), 
	&(options.grad_step_size), &(options.grad_stop_tol),
	&(options.far_far_away),  &(options.gap_open),
	&(options.gap_extend) , &(options.endgap),  &(options.far_away_cosine),
	&(options.H_length_mismatch_tol), &(options.H_length_mismatch_tol)};
    char * names_of_doubles[DOUBLES] = {
	"sheet_cosine",
	"tol", "z_max_store", "z_max_out",
	"F_guess_max", "F_eff_max",
	"z_max_corr", "z_min_compl",
	"grad_step_size", "grad_stop_tol",
	"far_far_away",  "gap_open",   
	"gap_extend", "endgap",  "far_away_cosine", "h_tol", "s_tol"};
    
    int * ptrs_to_ints[INTS] = {
	&(options.grid_size), &(options.number_maps_cpl),
	&(options.number_maps_out),
    	&(options.grad_max_step), &(options.exp_table_size),
	&(options.min_no_SSEs)};
    char * names_of_ints[INTS] = {
	"grid_size", "number_maps_cpl", "number_maps_out",
	"grad_max_step", "exp_table_size", "min_no_SSEs"};
    
    char * strings[STRINGS] = { options.outname, options.path, options.pdbf_tgt, options.pdbf_qry};
    char * names_of_strings[STRINGS] = {"outname", "path", "pdbf_tgt", "pdbf_qry"};

    char * ptrs_to_chars[CHARACTERS] = {&(options.chain_tgt), &(options.chain_qry)};
    char * names_of_chars[STRINGS]   = {"chain_tgt", "chain_qry"};
    
    fptr   = efopen ( filename, "r" );
    if (! fptr ) return 1;
    
 
    line_ctr = 0;
    memset ( line, 0, LONGSTRING);
    while(fgets(line,LONGSTRING,fptr)!=NULL){
	line_ctr++;
	/* tokenize */
	retval = tokenize ( token, &max_token, line, comment_char= '!' );
	switch ( retval ) {
	case  TOK_TOOMNY:
	    errmsg ( log, line_ctr, line, "\t\t %s\n", "Too many tokens.");
	    fclose (log);
	    break;
	case TOK_TOOLONG:
	    errmsg ( log, line_ctr, line, "\t\t %s\n", "Token too long.");
	    fclose (log);
	    break;
	}
	if ( max_token < 0 ) continue;
	
	
	token_assigned = 0;
	/* check first token for meaning (first token should be a keyword)*/
	for (ctr=0; ctr<DOUBLES && !token_assigned; ctr ++ ) {
	    if (  ! strcmp (token[0], names_of_doubles[ctr])  ) {
		if ( max_token < 1 ) {
		    errmsg ( log, line_ctr, line,
			     "\tKeyord %s should be followed by a real-valued number.\n",
			     token[0]);
		    return 1;
		}
		*(ptrs_to_doubles[ctr]) = atof ( token[1] );
		token_assigned = 1;
	    }
	}
	for (ctr=0; ctr<INTS && !token_assigned; ctr ++ ) {
	    if (  ! strcmp (token[0], names_of_ints[ctr])  ) {
		if ( max_token < 1 ) {
		    errmsg ( log, line_ctr, line,
			     "\tKeyord %s should be followed by an integer.\n",
			     token[0]);
		    return 1;
		}
		*(ptrs_to_ints[ctr]) = atoi ( token[1] );
		token_assigned = 1;
		break;
	    }
	}
	for (ctr=0; ctr<STRINGS && !token_assigned; ctr ++ ) {
	    if (  ! strcmp (token[0], names_of_strings[ctr])  ) {
		if ( max_token < 1 ) {
		    errmsg ( log, line_ctr, line,
			     "\tKeyord %s should be followed by a string.\n",
			     token[0]);
		    return 1;
		}
		sprintf ( strings[ctr], "%s", token[1] );
		token_assigned = 1;
		break;
	    }
	}

	for (ctr=0; ctr<CHARACTERS && !token_assigned; ctr ++ ) {
	    if (  ! strcmp (token[0], names_of_chars[ctr])  ) {
		if ( max_token < 1 ) {
		    errmsg ( log, line_ctr, line,
			     "\tKeyord %s should be followed by a character.\n",
			     token[0]);
		    return 1;
		}
		*(ptrs_to_chars[ctr]) = token[1][0];
		token_assigned = 1;
		break;
	    }
	}

	/* some hacking for switches */ 
	if ( ! token_assigned  &&  !strcmp (token[0], "verbose")  ) {
	    options.verbose = 1;
	    token_assigned  = 1;
	}

	if ( ! token_assigned  &&  !strcmp (token[0], "smith")  ) {
	    options.smith_waterman = 1;
	    token_assigned = 1;
	}

	if ( ! token_assigned  &&  !strcmp (token[0], "perp")  ) {
	    options.use_perp = 1;
	    token_assigned = 1;
	}
	if ( ! token_assigned  &&  !strcmp (token[0], "length")  ) {
	    options.use_length = 1;
	    token_assigned = 1;
	}
	if ( ! token_assigned  &&  !strcmp (token[0], "header")  ) {
	    options.print_header = 1;
	    token_assigned = 1;
	}
	if ( ! token_assigned  &&  !strcmp (token[0], "report_no_sse_overlap")  ) {
	    options.report_no_sse_overlap = 1;
	    token_assigned = 1;
	}
	if ( ! token_assigned  &&  !strcmp (token[0], "use_endgap")  ) {
	    options.use_endgap = 1;
	    token_assigned = 1;
	}
	if ( ! token_assigned  &&  !strcmp (token[0], "postp")  ) {
	    options.postprocess = 1;
	    token_assigned = 1;
	}

	
	if ( ! token_assigned){
	    errmsg ( log, line_ctr, line, "\t\t Keyword %s not recognized.\n", token[0]);
	    return 1;
	}    
	
	memset (line, 0, LONGSTRING);
    }
    fclose (fptr);
    
    /* if we are doing the postprocessing there should
       be two pdbfiles specified */
    if ( options.postprocess && (!options.pdbf_tgt[0] || !options.pdbf_qry[0]) ) {
	fprintf (log, "postp option should be accompanied by two pdb file names\n");
	fprintf (log,"(options pdbf1 and pdbf2, followed by the full path).\n");
	return 1;
    }

    
    return 0;
}

/****************************************************************************/
int errmsg ( FILE *fptr, int line_ctr, char line[LONGSTRING],
	     char * fmt, char * warnstr) {

    fprintf ( fptr, "\tError on line %3d:     %s", line_ctr, line);
    fprintf ( fptr, fmt, warnstr);
    return 0;
}
