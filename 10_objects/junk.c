    /* for all mapped blocks calculate similarity as exp (-d/5) */
    element_ctr_1 = -1;
    for (resctr1=0; resctr1<no_res_1; resctr1++) {
	
	memcpy (tmp, protein1->sequence[resctr1].pdb_id, PDB_ATOM_RES_NO_LEN*sizeof(char) );
	num_pdb_id = atoi (tmp);

	/* break if we are outside of the last element */
	if ( num_pdb_id > element_1_end [descr1->no_of_elements-1]) break; 
	
	if ( element_ctr_1 <0 ||  num_pdb_id>  element_1_end  [element_ctr_1] ) element_ctr_1++;
	
	if ( num_pdb_id >= element_1_begin[element_ctr_1] ) {
	    /* residue with number  resctr1 belongs to the element element_ctr_1 */
	    /* now, according to our mapping, what does element_ctr_1 map
	       onto in the other structure? (is there an element that it maps to?)*/
	    element_ctr_2 = map->x2y[element_ctr_1];
	    if (element_ctr_2 < 0) {
		/* actually, in this case 
		continue;
	    }

	    /* find the representative atom (CA) for  residue resctr1*/

	    /* translate to the frame with CM at the origin */

	    /* rotate the translation vector for element number element_ctr_2 */
	    
	    /* which residues belong to element_ctr_2 in the other structure? */
	    for (resctr2=element_2_begin[element_ctr_2]; resctr2<element_2_[element_ctr_2]; resctr2++) {
		/* find the representative atom (CA) for  residue resctr2*/

		/* translate to the frame with CM at the origin */

		/* rotate */

		/* translate out */
		
		/* find distance: */
		d = 0;
		for (i=0; i<3; i++) {
		    aux = ca1[i] - ca2[i];
		    d += aux*aux;
		}
		d = sqrt (d);
		similarity[resctr1][resctr2] = exp (-d/5.0);
	    }
	}
    }
