# include "struct.h"




# define DELETE  -789


int read_pdb ( char * pdbname,  char chain, Protein * protein, int postprocess) {

    /* TODO for the moment we rely on PDB annotation
       to extract structural elements - that should be changed */
    Residue * sequence = NULL;
    Helix   * helix = NULL, *helix_new = NULL;
    Strand  * strand = NULL, *strand_new = NULL;
    FILE    * fptr = NULL;
    int     * sse_sequence = NULL;
    int     ** sse_sub_elements = NULL;
    char line[BUFFLEN];
    char oldresno[PDB_ATOM_RES_NO_LEN+2];
    /* res name: 4 digits + insertion code + \0 */
    char oldrestype [PDB_ATOM_RES_NAME_LEN+2];
    char tmp[BUFFLEN], tmp2[BUFFLEN], *auxptr;
    char atomtypes_read_in[BUFFLEN];
    char old_chain;
    int atomctr, resctr,  no_res,ctr, nonblank;
    int helixctr, no_helices, helix_length;
    int helix_begin, helix_end;
    int strandctr, no_strands, strand_length;
    int strand_begin, strand_end;
    int reading_helix, reading_strand;
    int no_helices_after_cleanup;
    int no_strands_after_cleanup;
    int sse_ctr;
    int numerical_pdb_id;
    int retval;
    int ssc2, b1, e1, b2, e2, length;
    int a_container_exists;
    int overlap = 0;
    int return_value = 0, chain_found;
    int ca_trace, elmt_ctr;
    char single_letter ( char code[]);
    
    int has_backbone (Residue * sequence, int from, int to);
    
    /* open file */
    fptr = fopen ( pdbname, "r");
    if ( !fptr ) {
	fprintf (stderr, "Cno %s.\n", pdbname);
	return ERR_NO_FILE_OR_CHAIN;
    }


    
    /********************************************/
    /********************************************/
    /* cleanup                                  */
    memset (protein, 0, sizeof(Protein) );

    /********************************************/
    /********************************************/
    /* count residues                           */
    memset (line,  0, BUFFLEN);
    memset (oldresno, 0, PDB_ATOM_RES_NO_LEN+2);
    resctr = 0;
    chain_found = 0;
    old_chain = '\0';
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	
	if ( resctr ) {
	    if ( ! strncmp(line,"TER", 3) ||  ! strncmp(line,"END", 3) ||  line[PDB_ATOM_CHAINID] != old_chain)
		break;
	}
	if ( chain  && line[PDB_ATOM_CHAINID] != chain ) continue;
	chain_found  = 1;
	
	if( ! strncmp(line,"ATOM", 4) ||  ! strncmp(line,"HETATM", 6)){
	    
	    if (  strncmp (line+PDB_ATOM_RES_NO, oldresno,  PDB_ATOM_RES_NO_LEN+1) ) {
		
		strncpy (oldresno, line+PDB_ATOM_RES_NO, PDB_ATOM_RES_NO_LEN+1);
		oldresno[PDB_ATOM_RES_NO_LEN+1] = '\0';
		/* handling the case when the chain is not given, meaning: "take the first chain" */ 
		old_chain = line[PDB_ATOM_CHAINID];
		resctr ++;
	    }
	} 
	
    }

    /* sanity: */
    if ( chain && ! chain_found) {
	fprintf (stderr, "Chain %c not found in %s\n", chain, pdbname);
	return ERR_NO_FILE_OR_CHAIN;
    }
    
    no_res = resctr;
   
    /* allocate space */
    sequence = emalloc ( no_res*sizeof (Residue));
    if ( ! sequence ) return 1;
  
    
    /*********************************************/
    /*********************************************/
    /*   read in residue numbers and atom coords */
    rewind ( fptr);
    memset (line,  0, BUFFLEN);
    old_chain = '\0';

    memset (oldresno, 0, PDB_ATOM_RES_NO_LEN+2);
    memset (oldrestype, 0, PDB_ATOM_RES_NAME_LEN+2);
    resctr= -1;
    atomctr = 0;
    ca_trace = 1;
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	
	
	if ( resctr > -1) {
	    if  (! strncmp(line,"TER", 3) ||  ! strncmp(line,"END", 3)  ||  line[PDB_ATOM_CHAINID] != old_chain)
	    break;
	}
	if ( chain  && line[PDB_ATOM_CHAINID] != chain ) continue;
	
	if( ! strncmp(line,"ATOM", 4) ||  ! strncmp(line,"HETATM", 6)){
 	   /* if it's a hydrogen - skip */
	    if ( line[PDB_ATOM_ATOM_NAME] == 'H'
		 ||  line[PDB_ATOM_ATOM_NAME+1] == 'H') continue;
	    /* adjust the counters */ 
	    if (  strncmp (line+PDB_ATOM_RES_NO, oldresno,  PDB_ATOM_RES_NO_LEN+1) ) {
		/*+1 in  PDB_ATOM_RES_NO_LEN+1 means I am including the insertion code
		  in the identifier */
		strncpy (oldresno, line+PDB_ATOM_RES_NO, PDB_ATOM_RES_NO_LEN+1);
		strncpy (oldrestype, line+PDB_ATOM_RES_NAME, PDB_ATOM_RES_NAME_LEN);
		oldresno[PDB_ATOM_RES_NO_LEN+1]   = '\0';
		oldrestype[PDB_ATOM_RES_NAME_LEN] = '\0';
		
		/* handling the case when the chain is not given, meaning: "take the first chain" */ 
		old_chain = line[PDB_ATOM_CHAINID];
		
		resctr ++;
		if ( resctr >= no_res ) {
		    fprintf (stderr, "Error reading %s: rectr:%d   no res: %d\n", pdbname, resctr, no_res);
		    return ERR_NONSENSE;
		}
		atomctr = 0;
		/* keep track of atom types we have read in */
		memset (atomtypes_read_in, 0, BUFFLEN*sizeof(char));
		atomtypes_read_in[0] = '_';
		
		sequence[resctr].no_atoms = 1;
		strncpy ( sequence[resctr].pdb_id, oldresno, PDB_ATOM_RES_NO_LEN+2);
		sequence[resctr].pdb_id[PDB_ATOM_RES_NO_LEN+1]   = '\0';
		
		strncpy ( sequence[resctr].res_type, oldrestype, PDB_ATOM_RES_NAME_LEN+1);
		sequence[resctr].res_type[PDB_ATOM_RES_NAME_LEN] = '\0';
		sequence[resctr].res_type_short  =
		    single_letter ( sequence[resctr].res_type );
		/* modified residues are ok for the purposes here */
		/* unless they are sugars or some such - deal with it below */
		/* by checking the backbone atoms */
		///if ( !sequence[resctr].res_type_short ) return 1;
	   
	    } else {
		atomctr ++;
		sequence[resctr].no_atoms = atomctr + 1;
		if ( atomctr >= MAX_NO_ATOMS ) {
		    fprintf ( stderr, "Error parsing %s: I thought every aa has < %d atoms.\n",
			      pdbname, MAX_NO_ATOMS );
		    return ERR_MAX_ATOMS;
		}
	    }
	    /* read in atom info */
	    
	    auxptr = line+ PDB_ATOM_ATOM_NAME;
	    memset ( tmp, 0, PDB_ATOM_ATOM_NAME_LEN+1);
	    /* skip initial blanks*/
	    ctr  = 0;
	    while ( !(isalpha (*(auxptr + ctr))) &&
		    (ctr <= PDB_ATOM_ATOM_NAME_LEN) ) ctr++;
	    /* copy alphanum info */
	    nonblank = 0;
	    while (  isalnum (*(auxptr +ctr))  &&  (ctr <= PDB_ATOM_ATOM_NAME_LEN) ) {
		tmp[nonblank] =  *(auxptr +ctr);
		nonblank ++;
		ctr++;
	    }

	    /* have we already seen this atom type by any chance? */
	    tmp[nonblank] = '_';
	    if ( strstr (atomtypes_read_in, tmp) ) {
		/* ahould I check for an alt location code, or just move on? */
		//printf ( " %s >> %s  //// %s\n", sequence[resctr].pdb_id, atomtypes_read_in, tmp);
		continue;
	    } else {
		sprintf (atomtypes_read_in, "%s%s", atomtypes_read_in, tmp);
	    }
	    tmp[nonblank] = '\0';
	    
	    strncpy ( sequence[resctr].atom[atomctr].type, tmp, PDB_ATOM_ATOM_NAME_LEN );

	    /* is this a backbone atom?*/
	    sequence[resctr].atom[atomctr].backbone = 0;
	    if ( nonblank == 1) {
		  sequence[resctr].atom[atomctr].backbone =
		      !(  strcmp ( tmp, "N") && strcmp ( tmp, "C") && strcmp ( tmp, "O")  );
	    } else if (  nonblank == 2) {
		if (  ! strcmp ( tmp, "CA" )) {
		    sequence[resctr].atom[atomctr].backbone = 1;
		    sequence[resctr].Ca = sequence[resctr].atom+atomctr;
		}  else  {
		    sequence[resctr].atom[atomctr].backbone = 0;
		}
	    }
	    /* check if this is Ca trace */
	    if ( strcmp ( tmp, "CA" ) ) ca_trace = 0;
	    
	    strncpy ( tmp, line+PDB_ATOM_X, PDB_ATOM_X_LEN);
	    tmp[PDB_ATOM_X_LEN] = '\0';
	    sequence[resctr].atom[atomctr].x=atof(tmp);
	    strncpy ( tmp, line+PDB_ATOM_Y, PDB_ATOM_Y_LEN);
	    tmp[PDB_ATOM_Y_LEN] = '\0';
	    sequence[resctr].atom[atomctr].y=atof(tmp);
	    strncpy ( tmp, line+PDB_ATOM_Z, PDB_ATOM_Z_LEN);
	    tmp[PDB_ATOM_Z_LEN] = '\0';
	    sequence[resctr].atom[atomctr].z=atof(tmp);
	   
	}
    }

    if ( ca_trace) return ERR_SSE_NONE|ERR_CA_TRACE;
    
    /* clean PDB id tags from spaces */
    for (resctr=0; resctr < no_res; resctr ++ ) {
	retval = string_clean (sequence[resctr].pdb_id, PDB_ATOM_RES_NO_LEN+1);
	if ( retval ) {
	    fprintf (stderr, "Error in read_pdb(): empty id string for residue %d.\n", resctr);
	    return ERR_NONSENSE;
	}
    }

    /* store the sequence and its length */
    protein->sequence = sequence;
    protein->length   = no_res;

    /******************************************************/
    /******************************************************/
    /* if this is postprocessing, we are done right here  */
    if (postprocess) {
        /* we don't need o parse SSE info if this is a postprocessing run
	   - for pieces of pdb structure w/o header -
	    we handled that elsewhere */
   
	/* close file */
	fclose(fptr);
	return 0;
    }


    /*********************************************/
    /*********************************************/
    /* read in the HELIX/STRAND annotation       */
    rewind ( fptr);
    memset (line,  0, BUFFLEN);
    old_chain = '\0';
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	
	if (resctr && ( ! strncmp(line,"TER", 3) ||  ! strncmp(line,"END", 3)) ) break;
	if( ! strncmp(line,"ATOM", 4) ||  ! strncmp(line,"HETATM", 6)){
	    break; /* we have passed the annotation section*/
	    
	} else if(  ! strncmp(line,"HELIX", 5) ){
	    if ( chain  && line[PDB_HELIX_BEGIN_CHAIN]  != chain) continue;
	    if ( ! chain ) {
		/* handling the case when the chain is not given, meaning:
		   "take the first chain" */ 
		if ( !old_chain) {
		    old_chain = line[PDB_HELIX_BEGIN_CHAIN];
		} else if ( old_chain != line[PDB_HELIX_BEGIN_CHAIN] ) {
		    continue;
		}
	    }
	    memset ( tmp, 0, PDB_ATOM_RES_NO_LEN+1);
	    memset ( tmp2, 0, PDB_ATOM_RES_NO_LEN+1);
	    strncpy (tmp,  line+PDB_HELIX_BEGIN_ID, PDB_ATOM_RES_NO_LEN);
	    strncpy (tmp2, line+PDB_HELIX_END_ID,   PDB_ATOM_RES_NO_LEN);
	    
	    helix_begin  = atoi(tmp);
	    helix_end    = atoi(tmp2);
	    helix_length = helix_end - helix_begin + 1;
	    
	    if ( helix_length >= MIN_HELIX_LENGTH ) {
                /* annotate all residues in this range as belonging to a helix */
		for (resctr=0; resctr<no_res; resctr++) {
		    /* get rid of the insertions code, for the purposes here */
		    /* if helix starts, say at 133, the inserted residue 133A belongs to it too */
		    memset ( tmp, 0, PDB_ATOM_RES_NO_LEN+1);
		    strncpy (tmp,  sequence[resctr].pdb_id, PDB_ATOM_RES_NO_LEN);
		    numerical_pdb_id = atoi(tmp);
		    if ( helix_begin <= numerical_pdb_id  &&  numerical_pdb_id <= helix_end ) {
			sequence[resctr].belongs_to_helix = 1;
		    }
		}
	    }
	    
	} else if(  ! strncmp(line,"SHEET", 5) ){
	    
	    if ( chain  && line[PDB_SHEET_BEGIN_CHAIN] != chain) continue;
	    
	    if ( ! chain ) {
		/* handling the case when the chain is not given, meaning:
		   "take the first chain" */ 
		if ( !old_chain) {
		    old_chain = line[PDB_SHEET_BEGIN_CHAIN];
		} else if ( old_chain != line[PDB_SHEET_BEGIN_CHAIN] ) {
		    continue;
		}
	    }
	    memset ( tmp, 0, PDB_ATOM_RES_NO_LEN+1);
	    memset ( tmp2, 0, PDB_ATOM_RES_NO_LEN+1);
	    strncpy (tmp,  line+PDB_SHEET_BEGIN_ID, PDB_ATOM_RES_NO_LEN);
	    strncpy (tmp2, line+PDB_SHEET_END_ID, PDB_ATOM_RES_NO_LEN);

	    strand_begin  = atoi(tmp);
	    strand_end    = atoi(tmp2);
	    strand_length = strand_end - strand_begin + 1;
	    
	    if ( strand_length >= MIN_STRAND_LENGTH ) {
		/* annotate all residues in this range as belonging to a strand */
		for (resctr=0; resctr<no_res; resctr++) {
		    /* get rid of the insertions code, for the purposes here */
		    /* if strand starts, say at 133, the inserted residue 133A belongs to it too */
		    memset ( tmp, 0, PDB_ATOM_RES_NO_LEN+1);
		    strncpy (tmp,  sequence[resctr].pdb_id, PDB_ATOM_RES_NO_LEN);
		    numerical_pdb_id = atoi(tmp);
		    if ( strand_begin <= numerical_pdb_id  &&  numerical_pdb_id <= strand_end ) {
			sequence[resctr].belongs_to_strand = 1;
		    }
		}
	    }
	    

	}
	
    }
   
    /*********************************************/
    /*********************************************/
    /* count SSEs and do the sanity checking     */
    helixctr = 0;
    strandctr = 0;
    reading_helix = 0;
    reading_strand = 0;
    for (resctr=0; resctr<no_res; resctr++) {
	
	if ( sequence[resctr].belongs_to_helix ) {
	    if ( !reading_helix ) helixctr++;
	    reading_helix = 1;
	} else {
	    reading_helix = 0;
	}
	if ( sequence[resctr].belongs_to_strand ) {
	    if ( !reading_strand ) strandctr++;
	    reading_strand = 1;
	}  else {
	    reading_strand = 0;
	}
    }
   
    no_helices  = helixctr;
    no_strands  = strandctr;

    if ( no_helices==0 && no_strands==0 ) {
	fprintf ( stdout, "%d residues, no SSEs\n", no_res);
	return ERR_SSE_NONE;
    }
    
    /*********************************************/
    /*********************************************/
    /* allocate space  and store per-SSE info    */
    helix  = emalloc ( no_helices*sizeof (Helix));
    if ( ! helix ) return 1;
    strand  = emalloc ( no_strands*sizeof (Strand));
    if ( ! strand ) return 1;
    
    sse_sequence =  emalloc ( (no_strands+no_helices)*sizeof(int) );
    if ( !sse_sequence  ) return 1;

    /* for handling difficult cases below */
    sse_sub_elements = intmatrix (no_strands+no_helices, no_strands+no_helices+1);
    if ( !sse_sequence ) return 1;

    for (helixctr=0; helixctr<no_helices; helixctr++) {
	helix[helixctr].begin = -101;
	helix[helixctr].end   = -101;
    }
	    
    for (strandctr=0; strandctr<no_strands; strandctr++) {
	strand[strandctr].begin = -101;
	strand[strandctr].end   = -101;
    }
	    
    helixctr  = -1;
    strandctr = -1;
    sse_ctr = 0;
    reading_helix = 0;
    reading_strand = 0;
    for (resctr=0; resctr<no_res; resctr++) {
	
	if ( sequence[resctr].belongs_to_helix ) {
	    
	    if ( !reading_helix ) {
		helixctr++;
		helix[helixctr].begin = resctr;
		strncpy (helix[helixctr].begin_id, sequence[resctr].pdb_id, PDB_ATOM_RES_NO_LEN+2);
		sse_sequence[sse_ctr] = 2*helixctr;
		sse_ctr++;
	    }
	    reading_helix = 1;

	} else {
	    if ( reading_helix && resctr ) {
		helix[helixctr].end = resctr-1;
		strncpy (helix[helixctr].end_id, sequence[resctr-1].pdb_id, PDB_ATOM_RES_NO_LEN+2);
	    }
	    reading_helix = 0;
	}
	
	if ( sequence[resctr].belongs_to_strand ) {
	    if ( !reading_strand ) {
		strandctr++;
		strand[strandctr].begin = resctr;
		strncpy (strand[strandctr].begin_id, sequence[resctr].pdb_id, PDB_ATOM_RES_NO_LEN+2);
		sse_sequence[sse_ctr] = 2*strandctr+1; /* modulo will tell if it's a strand */
		sse_ctr++;
	    }
	    reading_strand = 1;
	    
	} else {
	    if ( reading_strand && resctr ) {
		strand[strandctr].end = resctr-1;
		strncpy (strand[strandctr].end_id, sequence[resctr-1].pdb_id, PDB_ATOM_RES_NO_LEN+2);
	    }
	    reading_strand = 0;
	}
    } 
    /* in case the sse reaches the end of the peptide chain */
    if ( reading_strand && resctr ) {
	strand[strandctr].end = resctr-1;
 	strncpy (strand[strandctr].end_id, sequence[resctr-1].pdb_id, PDB_ATOM_RES_NO_LEN+2);
    }
    if ( reading_helix  && resctr ) {
	helix[helixctr].end  = resctr-1;
	strncpy (helix[helixctr].end_id, sequence[resctr-1].pdb_id, PDB_ATOM_RES_NO_LEN+2);
    }
	    
    /* more sanity checking (or paranoid programming?) - do all elements have beginning and end? */
    for (helixctr=0; helixctr<no_helices; helixctr++) {
	if ( helix[helixctr].begin < -100 ) {
	    fprintf (stderr, "Beginning of helix %d (?).\n", helixctr+1);
	    return 1;
	}
	if ( helix[helixctr].end < -100 ) {
	    fprintf (stderr, "End of helix %d (?).\n", helixctr+1);
	    return 1;
	}
    }
    
    for (strandctr=0; strandctr<no_strands; strandctr++) {
	if ( strand[strandctr].begin < -100 ) {
	    fprintf (stderr, "Beginning of strand %d (?).\n", strandctr+1);
	    return 1;
	}
	if ( strand[strandctr].end < -100 ) {
	    fprintf (stderr, "End of strand %d (?).\n", strandctr+1);
	    return 1;
	}
    }

    /*********************************************/
    /*********************************************/
    /* cleanup and handling of difficult cases   */

    /* 1) check for minimum length               */
    for (sse_ctr=0; sse_ctr<no_strands+no_helices; sse_ctr++) {
	if ( sse_sequence[sse_ctr]%2 ) {
	    strandctr = sse_sequence[sse_ctr]/2;
	    strand_begin = strand[strandctr].begin;
	    strand_end  =  strand[strandctr].end;
	    length = strand_end - strand_begin +1;
	    strand[strandctr].length = length;
	    if ( length <MIN_STRAND_LENGTH) strand[strandctr].begin = DELETE;
	} else {
	    helixctr = sse_sequence[sse_ctr]/2;
	    helix_begin = helix[helixctr].begin;
	    helix_end   = helix[helixctr].end;
	    length = helix_end - helix_begin +1;
	    helix[helixctr].length = length;
	    if ( length <MIN_HELIX_LENGTH) helix[helixctr].begin = DELETE;

	}
    }

    /* One way to deal with HETATM notation which sometimes may */
    /* mean "a posttranlsationally modified residue" and sometimes */
    /* it may mean "a small ligand," is to make sure sure that */
    /* the purported residue has backbone atoms. (Ligand won't have them.) */
    /* Also, the way I am fitting onto strands means they */
    /* need to have the peptide bond discernible, therefore */
    /* 2) Make sure that the backbone atoms are present.*/
    for (sse_ctr=0; sse_ctr<no_strands+no_helices; sse_ctr++) {
	
	if ( sse_sequence[sse_ctr]%2 ) {
	    strandctr = sse_sequence[sse_ctr]/2;
	    strand_begin = strand[strandctr].begin;
	    if (  strand_begin == DELETE) continue;
	
	    strand_end  =  strand[strandctr].end;
	    if ( ! has_backbone (sequence, strand_begin, strand_end) )
		strand[strandctr].begin = DELETE;
	    
	} else {
	    
	    helixctr = sse_sequence[sse_ctr]/2;
	    helix_begin = helix[helixctr].begin;
	    if ( helix_begin == DELETE ) continue;
	    
	    helix_end   = helix[helixctr].end;
	    if ( ! has_backbone (sequence, helix_begin, helix_end) )
		helix[helixctr].begin = DELETE;
	    
	    
	}
    }
    
    /* 3) resolve overlaps                      */
    a_container_exists = 0;
    for (sse_ctr=0; sse_ctr<no_strands+no_helices; sse_ctr++) {
	int strandctr1 = 0, helixctr1 = 0, ssc1_is_helix = 0;
	if ( sse_sequence[sse_ctr]%2 ) {
	    ssc1_is_helix = 0;
	    strandctr1 = sse_sequence[sse_ctr]/2;
	    b1 = strand[strandctr1].begin;
	    e1 =  strand[strandctr1].end;
	} else {
	    ssc1_is_helix = 1;
	    helixctr1 = sse_sequence[sse_ctr]/2;
	    b1 = helix[helixctr1].begin;
	    e1 = helix[helixctr1].end;
	}
	if (b1 == DELETE) continue;
	
	for (ssc2=sse_ctr+1; ssc2<no_strands+no_helices; ssc2++) {
	    int strandctr2 = 0, helixctr2 = 0, ssc2_is_helix = 0;
	    if ( sse_sequence[ssc2]%2 ) {
		ssc2_is_helix = 0;
		strandctr2 = sse_sequence[ssc2]/2;
		b2 = strand[strandctr2].begin;
		e2 = strand[strandctr2].end;
	    } else {
		ssc2_is_helix = 1;
		helixctr2 = sse_sequence[ssc2]/2;
		b2 = helix[helixctr2].begin;
		e2 = helix[helixctr2].end;
	    }
	    if (b2 == DELETE) continue;
	    
	    if (e1 <= b2 +2 ) continue; /* no problem */
	    if (b2 < b1 ) {
		fprintf (stderr, "how did this happen.\n");
		return ERR_NONSENSE;
		
	    } else if ( e2 <= e1)  {
		/* element 2 is contained in 1 */
		ctr = (++sse_sub_elements[sse_ctr][0]); /* this is counter */
		sse_sub_elements[sse_ctr][ctr] = ssc2;
		a_container_exists = 1;
		return_value |= ERR_CONTAINER;
		
	    } else if ( e1 >= b2 + MAX_ALLOWED_OVERLAP)  { /* this is arbitrary, but cut out */
		                        /* the overlapping region - it is
					   usually something phishy */
		
		return_value |= ERR_OVERLAP; /* these are all recoverable errors,
						but the caller might want to know about them */

		if ( ssc2_is_helix ) {

		    helix[helixctr2].begin = e1 + 1;
		    /* are we too short now ? */
		    helix_begin = helix[helixctr2].begin;
		    helix_end   = helix[helixctr2].end;
		    length = helix_end - helix_begin +1;
		    if ( length <MIN_HELIX_LENGTH) helix[helixctr2].begin = DELETE;
		} else {
		    strand[strandctr2].begin = e1 + 1;
		    strand_begin = strand[strandctr2].begin;
		    strand_end  =  strand[strandctr2].end;
		    length = strand_end - strand_begin +1;
		    if ( length <MIN_STRAND_LENGTH) strand[strandctr2].begin = DELETE;
		}

		
		if ( ssc1_is_helix ) {
		    helix[helixctr1].end = b2 - 1;
		    /* are we too short now ? */
		    helix_begin = helix[helixctr1].begin;
		    helix_end   = helix[helixctr1].end;
		    length = helix_end - helix_begin +1;
		    if ( length <MIN_HELIX_LENGTH) {
			helix[helixctr1].begin = DELETE;
		        /* break out of ssc2 loop */
			break;
		    }
		} else {
		    strand[strandctr1].end = b2 - 1;
		    /* are we too short now ? */
		    strand_begin = strand[strandctr1].begin;
		    strand_end  =  strand[strandctr1].end;
		    length = strand_end - strand_begin +1;
		    if ( length <MIN_STRAND_LENGTH) {
			strand[strandctr1].begin = DELETE;
			/* break out of ssc2 loop */
			break;
		    }
		}
	    }
	}
	
    }

    /* it looks there are not many cases of containers (maybe 0.3% of PDB),
       and they seem to indicate either a pathological structure, or
       typo in the header - so for now just get rid of both the container
       and containee (we will return ERR_CONTAINER, however) */
    
    if ( a_container_exists ) {
	for (sse_ctr=0; sse_ctr<no_strands+no_helices; sse_ctr++) {
	    
	    if ( ! sse_sub_elements[sse_ctr][0] ) continue;
	    
	    if ( sse_sequence[sse_ctr]%2 ) {
		strandctr = sse_sequence[sse_ctr]/2;
		b1 = strand[strandctr].begin;
		e1 = strand[strandctr].end;
	    } else {
		helixctr = sse_sequence[sse_ctr]/2;
		b1 = helix[helixctr].begin;
		e1 = helix[helixctr].end;
	    }
	    
	    if (b1 == DELETE) continue;
	    
	    for (ctr=1; ctr <= sse_sub_elements[sse_ctr][0]; ctr++ ) {
		ssc2 = sse_sub_elements[sse_ctr][ctr];
		if ( sse_sequence[ssc2]%2 ) {
		    strandctr = sse_sequence[ssc2]/2;
		    strand[strandctr].begin = DELETE;
		} else {
		    helixctr = sse_sequence[ssc2]/2;
		    helix[helixctr].begin = DELETE;
		}
	    }
	}
	
    }

    /********************************************************************/
    /********************************************************************/
    /* after the cleanup, what is the new number of helices & strands? */
    no_helices_after_cleanup = 0;
    for (helixctr=0; helixctr<no_helices; helixctr++) {
	if (helix[helixctr].begin != DELETE ) {
	    no_helices_after_cleanup ++;
	}
    }
    no_strands_after_cleanup = 0;
    for (strandctr=0; strandctr<no_strands; strandctr++) {
	if (strand[strandctr].begin != DELETE ) {
	    no_strands_after_cleanup ++;
	}
    }
     
    if ( no_helices_after_cleanup==0 && no_strands_after_cleanup==0 ) {
	return_value |=  ERR_SSE_NONE;
	return return_value;
    }
    
    /********************/
    if ( no_helices_after_cleanup != no_helices ){
	helix_new  = emalloc ( no_helices_after_cleanup*sizeof (Helix));
	if ( ! helix_new ) return 1;

	/* repack */
	int hctr_new = 0;
	for (helixctr=0; helixctr<no_helices; helixctr++) {
	    if (helix[helixctr].begin != DELETE ) {
		memcpy ( helix_new + hctr_new, helix+helixctr, sizeof (Helix));
		hctr_new++;
	    }
	}
	
	/* free the old space and reassign the pointer */
	free (helix);
	helix = helix_new;
    }
    
    /********************/
    if ( no_strands_after_cleanup != no_strands ){
	strand_new  = emalloc ( no_strands_after_cleanup*sizeof (Strand));
	if ( ! strand_new  ) return 1;
	/* repack */
	int sctr_new = 0;
	for (strandctr=0; strandctr<no_strands; strandctr++) {
	    if ( strand[strandctr].begin != DELETE ) {
		memcpy ( strand_new + sctr_new, strand+strandctr, sizeof (Strand));
		sctr_new++;
	    }
	}

	/* free the old space and reassign the pointer */
	free (strand);
	strand = strand_new;
    }
    
    /********************/
    if ( no_helices_after_cleanup != no_helices ||
	 no_strands_after_cleanup != no_strands ){
	/* update the sse_sequence */
	
	free (sse_sequence);
	sse_sequence =  emalloc ( (no_strands_after_cleanup+no_helices_after_cleanup)*sizeof(int) );
	if ( !sse_sequence  ) return 1;
	
	sse_ctr=0;
	for (resctr=0; resctr<no_res; resctr++) {
	    for (helixctr=0; helixctr<no_helices_after_cleanup; helixctr++) {
		if ( helix[helixctr].begin == resctr ) {
		    sse_sequence[sse_ctr] = 2*helixctr;
		    sse_ctr++;
		    break;
		}
	    }
	    for (strandctr=0; strandctr<no_strands_after_cleanup; strandctr++) {
		if ( strand[strandctr].begin == resctr ) {
		    sse_sequence[sse_ctr] = 2*strandctr+1;
		    sse_ctr++;
		    break;
		}
	    }
	}
    }
    
    no_strands = no_strands_after_cleanup;
    no_helices = no_helices_after_cleanup;

    /********************************************************************/
    /********************************************************************/
    overlap = 0;
    for (sse_ctr=0; sse_ctr<no_strands+no_helices; sse_ctr++) {
	if ( sse_sequence[sse_ctr]%2 ) {
	    strandctr = sse_sequence[sse_ctr]/2;
	    b1 = strand[strandctr].begin;
	    e1 =  strand[strandctr].end;
	} else {
	    helixctr = sse_sequence[sse_ctr]/2;
	    b1 = helix[helixctr].begin;
	    e1 = helix[helixctr].end;
	}
	
	for (ssc2=sse_ctr+1; ssc2<no_strands+no_helices; ssc2++) {
	    if ( sse_sequence[ssc2]%2 ) {
		strandctr = sse_sequence[ssc2]/2;
		b2 = strand[strandctr].begin;
		e2 = strand[strandctr].end;
	    } else {
		helixctr = sse_sequence[ssc2]/2;
		b2 = helix[helixctr].begin;
		e2 = helix[helixctr].end;
	    }
	    if (b2 <= e1  ) {
		length = e1 - b2 + 1;
		if ( length > MAX_ALLOWED_OVERLAP ) {
		    printf ("overlap: b1:%d  b2:%d  e1:%d  e2:%d  -- len: %d\n",
		    	    b1, b2, e1, e2, length);
		    overlap = 1;
		}
		continue; /* no problem */
	    } else if (b2 < b1 ) {
		fprintf (stderr, "how did this happen.\n");
		return ERR_NONSENSE;
		
	    } else if ( e2 <= e1)  {
		/* element 2 is contained in 1 */
		fprintf (stderr, "how can I still have a container here?.\n");
		return ERR_NONSENSE;
	    }
	}
	
    }
# if 0
    printf ("chain: *%c* no residues: %d   helices: %d   strands: %d \n",
    	    chain, no_res, no_helics, no_strands);
    for (sse_ctr=0; sse_ctr<no_strands+no_helices; sse_ctr++) {
	if ( sse_sequence[sse_ctr]%2 ) {
	    strandctr = sse_sequence[sse_ctr]/2;
	    strand_begin = strand[strandctr].begin;
	    strand_end  =  strand[strandctr].end;
	    printf ("\t sse %3d is strand %2d:  %4d -- %4d (pdb: %5s -- %5s)\n",
		    sse_ctr+1,  strandctr, strand_begin, strand_end,
		    sequence[strand_begin].pdb_id,
		    sequence[strand_end].pdb_id );
	} else {
	    helixctr = sse_sequence[sse_ctr]/2;
	    helix_begin = helix[helixctr].begin;
	    helix_end   = helix[helixctr].end;
	    printf ("\t sse %3d is helix  %2d:  %4d -- %4d (pdb: %5s -- %5s)\n",
		    sse_ctr+1,  helixctr, helix_begin, helix_end,
		    sequence[helix_begin].pdb_id,
		    sequence[helix_end].pdb_id );
	}
    }
    exit (1);
# endif
    /*fix  belongs-to field */
    for (resctr=0; resctr < no_res; resctr ++ ) {
	sequence[resctr].belongs_to_helix      = -1;
	sequence[resctr].alt_belongs_to_helix  = -1;
	sequence[resctr].belongs_to_strand     = -1;
	sequence[resctr].alt_belongs_to_strand = -1;
    }

    for (strandctr=0; strandctr<no_strands; strandctr++) {
	strand_begin = strand[strandctr].begin;
	strand_end  =  strand[strandctr].end;
	for (resctr=strand_begin; resctr<=strand_end; resctr++) {
	    /* how many elements do we already belong to ? */
	    elmt_ctr = 0;
	    if ( sequence[resctr].belongs_to_helix >= 0 ) elmt_ctr++;
	    if ( sequence[resctr].alt_belongs_to_helix >= 0 ) elmt_ctr++;
	    if ( sequence[resctr].belongs_to_strand >= 0 ) elmt_ctr++;
	    if ( sequence[resctr].alt_belongs_to_strand >= 0 ) elmt_ctr++;

	    if ( elmt_ctr > 1 ) {
		fprintf (stderr,
			 "Error in %s: residue %s belongs to more than 2 elements (!?).\n",
			 pdbname, sequence[resctr].pdb_id);
		fprintf (stderr, " %d %d %d %d \n",  sequence[resctr].belongs_to_helix,
			 sequence[resctr].alt_belongs_to_helix,   sequence[resctr].belongs_to_strand,
			 sequence[resctr].alt_belongs_to_strand );
		
		return ERR_NONSENSE;
	    }
	    if ( sequence[resctr].belongs_to_strand < 0 ) {
		sequence[resctr].belongs_to_strand = strandctr;
	    } else if  (sequence[resctr].alt_belongs_to_strand < 0 ) {
		sequence[resctr].alt_belongs_to_strand = strandctr;
		
	    } 
	}
    }

    for (helixctr=0; helixctr<no_helices; helixctr++) {
	helix_begin = helix[helixctr].begin;
	helix_end  =  helix[helixctr].end;
	for (resctr=helix_begin; resctr<=helix_end; resctr++) {
	    /* how many elements do we already belong to ? */
	    elmt_ctr = 0;
	    if ( sequence[resctr].belongs_to_helix >= 0 ) elmt_ctr++;
	    if ( sequence[resctr].alt_belongs_to_helix >= 0 ) elmt_ctr++;
	    if ( sequence[resctr].belongs_to_strand >= 0 ) elmt_ctr++;
	    if ( sequence[resctr].alt_belongs_to_strand >= 0 ) elmt_ctr++;

	    if ( elmt_ctr > 1 ) {
		fprintf (stderr,
			 "Error in %s: residue %s belongs to more than 2 elements (!?).\n",
			 pdbname, sequence[resctr].pdb_id);
		fprintf (stderr, " %d %d %d %d \n",  sequence[resctr].belongs_to_helix,
			 sequence[resctr].alt_belongs_to_helix,   sequence[resctr].belongs_to_strand,
			 sequence[resctr].alt_belongs_to_strand );
		
		
		return ERR_NONSENSE;
	    }
	    
	    if ( sequence[resctr].belongs_to_helix < 0 ) {
		sequence[resctr].belongs_to_helix = helixctr;
	    } else if  (sequence[resctr].alt_belongs_to_helix < 0 ) {
		sequence[resctr].alt_belongs_to_helix = helixctr;
	    }
	}
    }

    free_imatrix (sse_sub_elements);



  
    /*****************************************************************/
    /* the pointers filled in here */
    protein->helix    = helix;
    protein->strand   = strand;
    protein->no_helices   = no_helices;
    protein->no_strands   = no_strands;
    protein->sse_sequence = sse_sequence;
   
    /* close file */
    fclose (fptr);
   
    return return_value;
}


/*****************************************************************/
/*****************************************************************/
/*****************************************************************/
int has_backbone (Residue * sequence, int from, int to ) {

    int resctr, atom_ctr, has_bb;
    int has_C, has_CA, has_N;
    Residue * residue;
    
    /* require that all residues have N, CA, and C */

    has_bb = 1;
    for (resctr=from; resctr<=to; resctr++) {
	
	residue = sequence+resctr;
	if ( residue->no_atoms < 3 ) {
	    has_bb = 0;
	    break;
	}
	has_C  = 0;
	has_CA = 0;
	has_N  = 0;
	for (atom_ctr = 0; atom_ctr < residue->no_atoms; atom_ctr++) {
	    if ( ! strcmp (residue->atom[atom_ctr].type, "CA")) {
		has_CA = 1;
	    } else if ( ! strcmp (residue->atom[atom_ctr].type, "C") ) {
		has_C = 1 ;
	    }  else if ( ! strcmp (residue->atom[atom_ctr].type, "N") ) {
		has_N = 1 ;
	    }
	  
	}
	if ( ! (has_CA && has_C && has_N) ) {
	    has_bb = 0;
	    break;
	}
	
    }

    return has_bb;

}
