
#include <exception>
#include <iostream>
#include <cstdio>
#include <cstring>
#include <cmath>
using namespace std;

#include "pdb.h"
#include "protein.h"

#include "mlexception.h"
#include "vec.h"
#include "line.h"
#include "linefit.h"
#include "geomdescr.h"

#define HELIX_CONST  1
#define STRAND_CONST 2

//enum SSEType {HELIX, STRAND};
//enum HelixType {UNKNOWN, THREETEN, ALPHA, PI};

static void caExtract(Residue* residues, int bgn, int end,
		      bool ca_valid[], Vec ca_pos[]) {
  for (int r = bgn; r <= end; r++) {
    ca_valid[r - bgn] = false;
    for (int a = 0; a < residues[r].no_atoms; a++) {
      Atom* atom = &(residues[r].atom[a]);
      if (strcmp(atom->type, "CA") == 0) {
	ca_valid[r - bgn] = true;
	ca_pos[r - bgn] = Vec(atom->x, atom->y, atom->z);
	break;
      }
    }
  }
}

/*
static HelixType cenGuessHelixType(int n_ca, bool ca_valid[], Vec ca_pos[]) {
  return ALPHA;
}
*/

static void cenCalcHelix( bool* valid, Vec* pos,
			  HelixType helix_type, Vec ca_pos[3]) {
  // Radii
  // collagen helix 1.6
  // 3-10 helix     1.9
  // alpha helix    2.3
  // pi helix       2.8

  double radius = 2.3;
  //double radius = 2.7;
  if (helix_type == HELIX_TYPE_310)
    radius = 1.9;
  if (helix_type == HELIX_TYPE_PI)
    radius = 2.8;
    
  // Should implement some checking here
  *valid = true;

  Vec spoke = radius * normalize(((ca_pos[2] + ca_pos[0]) / 2.0) - ca_pos[1]);
  *pos = ca_pos[1] + spoke;

  //double cos80 = cos((80.0 / 180.0) * PI);
  //*pos = ((((ca_pos[2] + ca_pos[0]) / 2.0) - ca_pos[1]) / (1 + cos80))
  //  + ca_pos[1];  
}


static void cenCalcStrand(bool* valid, Vec* pos, Vec ca_pos[3]) {

  // Should implement some checking here
  *valid = true;
  *pos = (ca_pos[0] + (2.0 * ca_pos[1]) + ca_pos[2]) / 4.0;
}


static void cenCalc(int sse_type,  int n,
		    bool ca_valid[], HelixType helix_type[], Vec ca_pos[],
		    bool cen_valid[], Vec cen_pos[]) {
  double ca_dist_min = 3.6;
  double ca_dist_max = 4.0;
  //double ca_dist_min = 3.0;
  //double ca_dist_max = 4.5;


  // The default
  for (int i = 0; i < n; i++)
    cen_valid[i] = false;

  // If sanity checks OK then dispatch to SSE specific function
  for (int i = 1; i < (n - 1); i++) {
    if ((!ca_valid[i-1]) ||
	(!ca_valid[i])   ||
	(!ca_valid[i+1]))
      continue;
    // C-alpha distance for trans peptide bond is around 3.8A
    double ca_dist;
    ca_dist = (ca_pos[i] - ca_pos[i-1]).len();
    if ((ca_dist > ca_dist_max) || (ca_dist < ca_dist_min))
      continue;
    ca_dist = (ca_pos[i+1] - ca_pos[i]).len();
    if ((ca_dist > ca_dist_max) || (ca_dist < ca_dist_min))
      continue;

    if (sse_type == HELIX) {
      cenCalcHelix(&(cen_valid[i]), &(cen_pos[i]),
		   helix_type[i], &(ca_pos[i-1]));
      continue;
    }

    if (sse_type == STRAND) {
      cenCalcStrand(&(cen_valid[i]), &(cen_pos[i]), &(ca_pos[i-1]));
      continue;
    }
  }
}


static void drawJmolPoints(const char * prefix, int id,
			   int n_points, Vec points[]) {

  for (int b = 0; b < (n_points - 1); b++) {
    Vec* bgn = &(points[b]);
    Vec* end = &(points[b+1]);
    cout << "draw " 
	 << prefix << id << "X" << b 
	 << " cylinder";
    cout << " {" << bgn->x << " " << bgn->y << " " << bgn->z << "}";
    cout << " {" << end->x << " " << end->y << " " << end->z << "}";
    cout << " width 0.2;";
    cout << "color $"
	 << prefix << id << "X" << b 
	 << " white" << endl;
  }

}


void geomdescrCalc(Protein *protein, bool jmolpoints) throw(MLException) {

  Residue* residues = protein->sequence;
  int n_residues = protein->length;

  //if (jmolpoints) {
  //  cout << "load 9pai.pdb" << endl
  //	 << "spacefill OFF" << endl
  //	 << "wireframe OFF" << endl
  //	 << "color GROUP" << endl
  //	 << "backbone" << endl;
  //}

  for (int h=0; h < protein->no_helices; h++) {
    Helix *helix = &(protein->helix[h]);
    helix->n_lines = 0;
    helix->lines = new Line[0];

    int bgn_res = helix->begin - 1;
    if (bgn_res < 0)
      bgn_res = 0;
    int end_res = helix->end + 1;
    if (end_res >= n_residues)
      end_res = n_residues - 1;

    int n_ca = end_res - bgn_res + 1;
    bool ca_valid[n_ca];
    Vec ca_pos[n_ca];
    HelixType ca_helix_type[n_ca];

    caExtract(residues, bgn_res, end_res, ca_valid, ca_pos);
    for (int c = 0; c < n_ca; c++)
      ca_helix_type[c] = residues[bgn_res + c].helix_type;
    
    bool cen_valid[n_ca];
    Vec cen_pos[n_ca];

    cenCalc(HELIX, n_ca, ca_valid, ca_helix_type, ca_pos,
	    cen_valid, cen_pos);

    int n_points = 0;
    for (int c = 0; c < n_ca; c++) {
      if (cen_valid[c]) {
	n_points++;
      }
    }
    int ca_index[n_points];
    Vec points[n_points];
    
    int p = 0;
    for (int c = 0; c < n_ca; c++) {
      if (cen_valid[c]) {
	ca_index[p] = c;
	points[p] = cen_pos[c];
	p++;
      }
    }

    int min_line_len = 3;
    if (n_points < min_line_len)
      continue;

    if (jmolpoints)
      drawJmolPoints("helix", h, n_points, points);

    LineFit linefit(min_line_len, n_points, points);

    int n_lines;
    Line* lines = linefit.find(&n_lines, n_points, points);
    if (lines == 0) {
      cerr << "Unable to fit points to line: " << n_points << endl;
      continue;
    }

    helix->n_lines = n_lines;
    delete[] ((Line*) helix->lines);
    helix->lines = new Line[n_lines];
    copyToArray((Line*) helix->lines, n_lines, lines);

    for (int i = 0; i < n_lines; i++) {
      Line* line = &(((Line*) (helix->lines))[i]);
      line->bgn = ca_index[line->bgn] + bgn_res;
      line->end = ca_index[line->end] + bgn_res;
    }
  }

  for (int s=0; s < protein->no_strands; s++) {
    Strand *strand = &(protein->strand[s]);
    strand->n_lines = 0;
    strand->lines = new Line[0];

    int bgn_res = strand->begin - 1;
    if (bgn_res < 0)
      bgn_res = 0;
    int end_res = strand->end + 1;
    if (end_res >= n_residues)
      end_res = n_residues - 1;

    int n_ca = end_res - bgn_res + 1;
    bool ca_valid[n_ca];
    Vec ca_pos[n_ca];

    caExtract(residues, bgn_res, end_res, ca_valid, ca_pos);

    if (n_ca == 3) {
      if (!(ca_valid[0] && ca_valid[1] && ca_valid[2]))
	continue;
      strand->n_lines = 1;
      delete[] ((Line*) strand->lines);
      strand->lines = new Line[1];
      Line* line = (Line*) strand->lines;
      line->setTSD(1, 0.0);
      line->bgn = strand->begin;
      line->end = strand->end;
      line->pos = (ca_pos[0] + ca_pos[1] + ca_pos[2]) / 3.0;
      line->dir = normalize(ca_pos[2] - ca_pos[0]);
      continue;
    }

    if (n_ca == 4) {
      if (!(ca_valid[0] && ca_valid[1] && ca_valid[2] && ca_valid[3]))
	continue;
      strand->n_lines = 1;
      delete[] ((Line*) strand->lines);
      strand->lines = new Line[1];
      Line* line = (Line*) strand->lines;
      line->setTSD(1, 0.0);
      line->bgn = strand->begin;
      line->end = strand->end;
      line->pos = (ca_pos[0] + ca_pos[1] + ca_pos[2] +ca_pos[3]) / 4.0;
      line->dir = normalize((ca_pos[2] - ca_pos[0]) + (ca_pos[3] - ca_pos[1]));
      continue;
    }

    bool cen_valid[n_ca];
    Vec cen_pos[n_ca];

    cenCalc(STRAND, n_ca, ca_valid, 0, ca_pos, cen_valid, cen_pos);

    int n_points = 0;
    for (int c = 0; c < n_ca; c++) {
      if (cen_valid[c]) {
	n_points++;
      }
    }
    int ca_index[n_points];
    Vec points[n_points];
    
    int p = 0;
    for (int c = 0; c < n_ca; c++) {
      if (cen_valid[c]) {
	ca_index[p] = c;
	points[p] = cen_pos[c];
	p++;
      }
    }

    int min_line_len = 3;
    if (n_points < min_line_len)
      continue;

    if (jmolpoints)
      drawJmolPoints("strand", s, n_points, points);

    LineFit linefit(min_line_len, n_points, points);

    int n_lines;
    Line* lines = linefit.find(&n_lines, n_points, points);
    if (lines == 0) {
      cerr << "Unable to fit points to line: " << n_points << endl;
      continue;
    }

    strand->n_lines = n_lines;
    delete[] ((Line*) strand->lines);
    strand->lines = new Line[n_lines];
    copyToArray((Line*) strand->lines, n_lines, lines);

    for (int i = 0; i < n_lines; i++) {
      Line* line = &(((Line*) (strand->lines))[i]);
      line->bgn = ca_index[line->bgn] + bgn_res;
      line->end = ca_index[line->end] + bgn_res;
    }
  }

}


/////////////
//  Print  //
/////////////

static void printLines(FILE* fptr, Residue* residues,
		       const char * sse_name, int sse_const,
		       int n_lines, Line* lines) {
  for (int i = 0; i < n_lines; i++) {
    Line* line = &(lines[i]);
    fprintf (fptr, "%-6s  %3d  %3d %5d  %5s %5s",
	     sse_name, sse_const,
	     0, // strand->sheet_number
	     line->end - line->bgn + 1, 
	     residues[line->bgn].pdb_id,
	     residues[line->end].pdb_id);
    fprintf(fptr, " %8.2lf %8.2lf %8.2lf",
	    line->dir.x, line->dir.y, line->dir.z);
    fprintf(fptr, " %8.2lf %8.2lf %8.2lf",
	    line->pos.x, line->pos.y, line->pos.z);
    fprintf(fptr, "\n");
  }
}

void geomdescrPrint(Protein *protein, const char *name, const char *outfile) 
  throw(MLException) {

  FILE* fptr = fopen(outfile, "w");
  if (fptr == 0)
    throw(MLException("geomdescrPrint(): unable to open file"));

  fprintf (fptr, "name: %4s\n", name);
  fprintf (fptr, "number of residues: %d\n", protein->length);

  int n_helices = 0;
  for (int h = 0; h < protein->no_helices; h++)
    n_helices += protein->helix[h].n_lines;
  fprintf (fptr, "number of helices:  %d\n", n_helices);

  int n_strands = 0;
  for (int h = 0; h < protein->no_strands; h++)
    n_strands += protein->strand[h].n_lines;
  fprintf (fptr, "number of strands:  %d\n", n_strands);

  int n_tot_sse = protein->no_strands + protein->no_helices;
  for (int s = 0; s < n_tot_sse; s++) {
    // Strand
    if ((protein->sse_sequence[s] % 2) == 1) {
      int index = protein->sse_sequence[s] / 2;
      int n_lines = protein->strand[index].n_lines;
      Line *lines = (Line*) (protein->strand[index].lines);
      printLines(fptr, protein->sequence,
		 "STRAND", STRAND_CONST,
		 n_lines, lines);
    }
    // Helix
    else {
      int index = protein->sse_sequence[s] / 2;
      int n_lines = protein->helix[index].n_lines;
      Line *lines = (Line*) (protein->helix[index].lines);
      printLines(fptr, protein->sequence,
		 "HELIX", HELIX_CONST,
		 n_lines, lines);
    }
  }

  fprintf(fptr, "#\n");
  fclose (fptr);
}


/*
static HelixType cenGuessHelixType(int n_ca, bool ca_valid[], Vec ca_pos[]) {

  int counts[3] = {0, 0, 0};
  double dists[3] = {0.0, 0.0, 0.0};
  for (int sep = 3; sep <= 5; sep++) {
    for (int ca = 0; ca < (n_ca - sep); ca++) {
      if ((!ca_valid[ca]) || (!ca_valid[ca + sep]))
	continue;
      double dist = (ca_pos[ca + sep] - ca_pos[ca]).len();
      if (dist > 20.0)
	continue;
      dists[sep - 3] += dist;
      counts[sep - 3]++;
    }
  }

  for (int i = 0; i < 3; i++) {
    if (counts[i] == 0) {
      cout << "Defaulted to ALPHA" << endl;
      return ALPHA;
    }
    dists[i] /= (double) counts[i];
    cerr << dists[i] << " " << counts[i] << endl;
  }

  double smallest_dist = dists[0];
  int smallest_sep = 3;
  for (int i = 1; i < 3; i++) {
    if (dists[i] < smallest_dist) {
      smallest_dist = dists[i];
      smallest_sep = i + 3;
    }
  }

  switch (smallest_sep) {
  case 3:
    cout << "Guessed THREETEN" << endl;
    return THREETEN;
  case 4:
    cout << "Guessed ALPHA" << endl;
    return ALPHA;
  case 5:
    cout << "Guessed PI" << endl;
    return PI;
  }

  cerr << "cenGuessHelixType(): failed" << endl;
  return UNKNOWN;
}
*/
