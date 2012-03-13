
#include <stdio.h>
#include <string.h>

#include "pdb.h"
#include "protein.h"

#include "stride.h"

#include "calcstride.h"

static int findResidue(char * resnum, char * restype,
		       Residue* residues, int bgn, int n_residues,
		       const char * pdbfile) {
  int r;
  for (r = bgn; r < n_residues; r++) {
    if (strcmp(resnum, residues[r].pdb_id) == 0) {
      if (strcmp(restype, residues[r].res_type) != 0)
	fprintf(stderr,
		"WARNING: calcStride(): %s: res-type mismatch: '%s' '%s'\n",
		pdbfile, restype, residues[r].res_type);
      break;
    }
  }
  if (r == n_residues) {
    fprintf(stderr,
	    "WARNING: calcStride(): %s: residue not found: '%s' (%s)\n",
	    pdbfile, resnum, restype);
    /*
    for (r = bgn; r < n_residues; r++)
      fprintf(stderr, "'%s' (%s)\n", residues[r].pdb_id, residues[r].res_type);
    */
    return -1;
  }
  return r;
}    


const char * calcStride(Residue* residues, int n_residues,
			FILE* fptr, const char * pdbfile,
			char chain) {
  StrideData sd;
  const char * errmsg;
  char chains[2];
  StrideSse* sse;
  int r, s;
  int bgn, end;
  int n_sses_assigned;
  HelixType helix_type;


  if (chain == '\0')
    return "calcStride(): chain == '\0'";
  chains[0] = (chain != ' ') ? chain : '\0';
  chains[1] = '\0';

  strideInit(&sd);
  if (sd.errmsg)
    return sd.errmsg;

  strideRead(&sd, chains, fptr, pdbfile);
  if (sd.errmsg) {
    errmsg = sd.errmsg;
    strideDone(&sd);
    return errmsg;
  }

  /*
  for (r = 0; r < n_residues; r++)
    fprintf(stderr, "'%s'\n", residues[r].pdb_id); 
  */

  n_sses_assigned = 0;

  for (s = 0; s < sd.n_helices; s++) {
    sse = &(sd.helices[s]);

    if (sse->chain_id != chain)
      continue;

    bgn = findResidue(sse->bgn_resnum, sse->bgn_restype,
		      residues, 0, n_residues, pdbfile);
    if (bgn < 0)
      continue;
    end = findResidue(sse->end_resnum, sse->end_restype,
		      residues, bgn, n_residues, pdbfile);
    if (end < 0)
      continue;

    switch (sse->subtype) {
    case 'H':
      helix_type = HELIX_TYPE_ALPHA;
      break;
    case 'G':
      helix_type = HELIX_TYPE_310;
      break;
    case 'I':
      helix_type = HELIX_TYPE_PI;
      break;
    default:
      helix_type = HELIX_TYPE_UNDEF;
    }

    if ((end - bgn + 1) >= MIN_HELIX_LENGTH) {
      for (r = bgn; r <= end; r++) {
	residues[r].belongs_to_helix = 1;
	residues[r].helix_type = helix_type;
      }
      n_sses_assigned++;
    }
  }

  for (s = 0; s < sd.n_strands; s++) {
    sse = &(sd.strands[s]);

    if (sse->chain_id != chain)
      continue;

    bgn = findResidue(sse->bgn_resnum, sse->bgn_restype,
		      residues, 0, n_residues, pdbfile);
    if (bgn < 0)
      continue;
    end = findResidue(sse->end_resnum, sse->end_restype,
		      residues, bgn, n_residues, pdbfile);
    if (end < 0)
      continue;

    if ((end - bgn + 1) >= MIN_STRAND_LENGTH) {
      for (r = bgn; r <= end; r++) {
	residues[r].belongs_to_strand = 1;
	residues[r].helix_type = HELIX_TYPE_UNDEF;
      }
      n_sses_assigned++;
    }
  }

  strideDone(&sd);
  if (sd.errmsg)
    return sd.errmsg;

  if (n_sses_assigned == 0)
    return "calcStride(): no SSEs assigned";

  return 0;
}
