

#include <exception>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <unistd.h>
using namespace std;

#include "pdb.h"
#include "protein.h"
#include "readpdb.h"

#include "mlexception.h"
#include "geomdescr.h"


static void printUsage(const char * progname) {
  fprintf(stderr, "Usage: %s [-j] <name>  <in file>  <out file>  [<chain>]\n",
	  progname);
}

int main(int argc, char *argv[]) {
  const char * progname = argv[0];

  bool jmolpoints = false;
  for (;;) {
    int c = getopt (argc, argv, "j");
    if (c == -1)
      break;
    switch (c) {
    case 'j':
      jmolpoints = true;
      break;
    case '?':
      printUsage(progname);
      exit(1);
    default:
      fprintf(stderr, "getopt() error\n");
      exit(1);
    }
  }
  int n_options = argc - optind;
  if (!((n_options == 3) || (n_options == 4))) {
    printUsage(progname);
    exit(1);
  }
  const char *name = argv[optind+0];
  const char *infile = argv[optind+1];
  const char *outfile = argv[optind+2];
  char chain = (n_options > 3) ? argv[optind+3][0] : '\0';

  Protein protein;
  int readpdb_status = readPDB(infile, chain, &protein);
  if (readpdb_status & (ERR_SSE_NONE|ERR_NO_FILE_OR_CHAIN|ERR_NONSENSE)) {
    fprintf(stderr, "Parsing failed for %s: %d\n", name, readpdb_status);
    exit(readpdb_status);
  }

  try {

    geomdescrCalc(&protein, jmolpoints);
  
    geomdescrPrint(&protein, name, outfile);

  }
  catch (MLException& e) {
    fprintf(stderr, "Line fitting failed: ");
    fprintf(stderr, "%s", e.what());
    fprintf(stderr, "\n");
    exit(e.retval);
  }

  return readpdb_status;
}
