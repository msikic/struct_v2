
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>

#include "stride.h"

int main(int argc, char* argv[]) {
  char * pdbfile;
  char * chains;
  FILE *file;
  StrideData sd;
  int i;
  StrideSse* sse;

  if ((argc != 2) && (argc != 3)) {
    fprintf(stderr, "Usage: %s pdbfile [chains]\n", argv[0]);
    exit(EXIT_FAILURE);
  }
  pdbfile = argv[1];
  chains = "";
  if (argc == 3)
    chains = argv[2];
  
  strideInit(&sd);
  if (sd.errmsg) {
    fprintf(stderr, "strideInit() failed: %s\n", sd.errmsg);
    exit(EXIT_FAILURE);
  }

  file = fopen(pdbfile, "r");
  if (file == 0) {
    fprintf(stderr, "Unable to open '%s': %s\n", pdbfile, strerror(errno));
    exit(EXIT_FAILURE);
  }

  strideRead(&sd, chains, file, pdbfile);
  if (sd.errmsg) {
    fprintf(stderr, "strideRead() failed: %s\n", sd.errmsg);
    exit(EXIT_FAILURE);
  }

  if (fclose(file) != 0) {
    fprintf(stderr, "Unable to close '%s'\n", pdbfile);
    exit(EXIT_FAILURE);
  }

  for (i = 0; i < sd.n_helices; i++) {
    sse = &(sd.helices[i]);
    printf("HELIX(%c) '%c' %4s(%3s) %4s(%3s)\n",
	   sse->subtype,
	   sse->chain_id,
	   sse->bgn_resnum,
	   sse->bgn_restype,
	   sse->end_resnum,
	   sse->end_restype);
  }

  for (i = 0; i < sd.n_strands; i++) {
    sse = &(sd.strands[i]);
    printf("STRAND(%c) '%c' %4s(%3s) %4s(%3s)\n",
	   sse->subtype,
	   sse->chain_id,
	   sse->bgn_resnum,
	   sse->bgn_restype,
	   sse->end_resnum,
	   sse->end_restype);
  }

  strideDone(&sd);
  if (sd.errmsg) {
    fprintf(stderr, "strideDone() failed: %s\n", sd.errmsg);
    exit(EXIT_FAILURE);
  }

  return EXIT_SUCCESS;
}
