
// W.A. Sherman, 2010
// An overview of how this all works:
//
// 1st:
// LineFit(int min_line_len, int n_points, Vec points[]);
// Among other things, this allocates the cache of
// previously calculated results.
// - Cache the fit of each line segment as seg_cache[bgn][end]
// - Cache the best list of segments fitting
//   the remaining points as list_cache[n_lines][bgn]
//
// 2nd:
// LineFit.find(int* n_lines, int n_points, Vec points[]);
// The core loop is:
//     n_lines = 1;
//     for (;;) {
//         Line* lines = findList(0, n_lines);
//         if (n_lines >= n_lines_max)
//             return lines;
//         if (lines->getTSD() <= tsd_max)
//             return lines;
//         n_lines++;
//     }
// Basically, first try to fit all the points with
// just one line - if the total squared deviation ("TSD")
// is OK then we're done.
// Otherwise, try to fit all the points with
// two lines - if TSD is OK we're done.
// ...and so on until we reach the maximum possible
// number of lines - then we just return the best fit
// regardless of what the TSD is.
//
// 3rd:
// LineFit.findList(int bgn1, int n_lines)
// The core of the routine is:
//     if (n_lines == 1) {
//         return LineFit.findSeg(bgn1, n_points - 1);
//     }
//     else {
//         while (bgn2 <= bgn2_max) {
//             seg  = findSeg(bgn1, bgn2 - 1);
//             list = findList(bgn2, n_lines-1);
//             lines = append(seg, list);
//             if (lines.tsd < best_lines.tsd)
//                 best_lines = lines;
//             bgn2++;
//         }
//         return best_lines;
//     }
// The basic idea is that we split the points
// and fit the first half to a single line
// and then find the best fit to the second half
// but using only (n_lines - 1).
// We iterate over all possible splits and return
// the optimal split
// Note that LineFit.findList() is called recursively
// and that when we get to one line left we need to fit
// all the way to the end (use up the remaining points).
// Also, because the results are cached, in most
// cases it's not necessary to recurse all the way
// down to (n_lines==1).
//
// 4th:
// LineFit.findSeg(int bgn, int end);
// This one is simple: if the segment is already
// calculated then just return that - otherwise
// call Line.calcSeg() to do the heavy lifting
// (and then also store the result in the cache).
//
// 5th:
// LineFit.calcSeg(int bgn, int end);
// This is where we actually fit lines to points.
// For fitting an individual line segment, we want
// the distances to be perpendicular to the line
// so this is actually a non-trivial problem.
// The solution is called Total Least-Squares (or ODR)
// The first step is to translate to the centroid
// Then we use and algorithm based on eigen values
// to get the direction (can also use SVD algorithm).
//
// 6th:
// LineFit.getXXX();
// LineFit.setXXX();
// Finally, in order to do our dynamic programming
// we need to be able to show how the path is all
// linked up. We use some linked lists that are
// indexed by number of lines (n_lines).
// The get/set routines make keeping track of this
// a bit easier (e.g we can use the index of
// the total number of lines in the list).

#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>
using namespace std;

#include "mlexception.h"
#include "vec.h"
#include "line.h"
#include "linefit.h"

LineFit::LineFit(int min_line_len, int n_points, Vec points[])
  throw(MLException) {

  if (n_points < min_line_len)
    throw(MLException("n_points < min_line_len"));

  this->min_line_len = min_line_len;
  this->n_points = n_points;
  this->points = points;

  n_lines_max = n_points / min_line_len;

  list_cache = new Line*[n_lines_max * n_points];
  for (int nl = 0; nl < n_lines_max; nl++) {
    for (int p = 0; p < n_points; p++) {
      list_cache[nl*n_points + p] = 0;
    }
  }
  
  seg_cache  = new Line*[n_points * n_points];
  for (int bgn = 0; bgn < n_points; bgn++) {
    for (int end = 0; end < n_points; end++) {
      seg_cache[bgn*n_points + end] = 0;
    }
  }
}

LineFit::~LineFit() {
  for (int bgn = 0; bgn < n_points; bgn++) {
    for (int end = 0; end < n_points; end++) {
      Line* line = seg_cache[bgn*n_points + end];
      if (line != 0)
	delete line;
    }
  }
  delete[] list_cache;
  delete[] seg_cache;
}

void LineFit::checkCachedList(int n_lines, int bgn) throw(MLException) {
  if (n_lines < 1)
    throw(MLException("LineFit::checkCachedList(): n_lines < 1"));
  if (n_lines > n_lines_max)
    throw(MLException("LineFit::checkCachedList(): n_lines > n_lines_max"));
  if (bgn < 0)
    throw(MLException("LineFit::checkCachedList(): bgn < 0"));
  if (bgn >= n_points)
    throw(MLException("LineFit::checkCachedList(): end >= n_points"));
}

Line* LineFit::getCachedList(int n_lines, int bgn) throw(MLException) {
  checkCachedList(n_lines, bgn);
  return list_cache[(n_lines - 1) * n_points + bgn];
}

Line* LineFit::setCachedList(int n_lines, int bgn, Line* lines)
  throw(MLException) {
  checkCachedList(n_lines, bgn);
  list_cache[(n_lines - 1) * n_points + bgn] = lines;
  return lines;
}

void LineFit::checkCachedSeg(int bgn, int end) throw(MLException) {
  if (bgn < 0)
    throw(MLException("LineFitCache::checkSeg(): bgn < 0"));
  if (end >= n_points)
    throw(MLException("LineFitCache::checkSeg(): end >= n_points"));
  if (bgn >= end)
    throw(MLException("LineFitCache::checkSeg(): bgn >= end"));
}

Line* LineFit::getCachedSeg(int bgn, int end) throw(MLException) {
  checkCachedSeg(bgn, end);
  return seg_cache[bgn*n_points + end];
}

Line* LineFit::setCachedSeg(int bgn, int end, Line* line) throw(MLException) {
  checkCachedSeg(bgn, end);
  seg_cache[bgn*n_points + end] = line;
  return line;
}

extern "C" void dsyev_(char * jobz, char * uplo, int* N,
		       double * A, int * leading_dim,
		       double * eigenvalues,
		       double *workspace, int *workspace_size,
		       int * retval);

Line* LineFit::calcSeg(int bgn, int end) throw(MLException) {

  // Don't bother to sanity check bgn/end
  // because expect to be called from get()

  // Allocate line (the return val)
  Line* line = new Line(n_lines_max, bgn, end);

  // Special case: only 2 points
  if ((end - bgn) == 1) {
    line->setTSD(1, 0.0);
    line->pos = (points[bgn] + points[end]) / 2.0;
    line->dir = (points[end] - points[bgn]).normalize();
    return line;
  }

  int lenseg = end - bgn + 1;
  
  // Centroid
  Vec pos = Vec(0.0, 0.0, 0.0);
  for (int p = bgn; p <= end; p++)
    pos += points[p];
  pos /= (double) lenseg; 

  // find the "moments of inertia"
  double I[3][3] = {{0.0, 0.0, 0.0},
		    {0.0, 0.0, 0.0},
		    {0.0, 0.0, 0.0}};
  for (int p = bgn; p <= end; p++) {
    double my_point[3] = {points[p].x - pos.x,
			  points[p].y - pos.y,
			  points[p].z - pos.z};
    for (int x = 0; x < 3; x++) {  /* modulo = circular permutation */
      I[x][x] += my_point[(x+1)%3]*my_point[(x+1)%3] +
	my_point[(x+2)%3]*my_point[(x+2)%3];
      for (int y = x+1; y < 3; y++) { /* offdiag elements */
	I[x][y] -= my_point[x]*my_point[y];
      }
    }
  }
  for (int x = 0; x < 3; x++) { 
    for (int y = x+1; y < 3; y++) {
      I[y][x] =  I[x][y];
    }
  }

  /*****************************************/
  /* diagonalize I[][], pick the direction
     with the smallest moment of inertia,
     and rotate back to the initial frame  */
  /*****************************************/
  char jobz = 'V'; /* find evalues and evectors */
  char uplo = 'L'; /* amtrix is stored as lower (fortran convention) */
  int  N = 3; /* the order of matrix */
  int leading_dim = N;
  int retval;
  double A[N*N];
  double eigenvalues[N];
  double workspace[3*N];
  int workspace_size = 3*N;

  for (int x = 0; x < 3; x++) {
    for (int y=0; y < 3; y++) {
      A[x*3+y] = I[x][y];
    }
  }
   
  dsyev_ ( &jobz, &uplo, &N, A,  &leading_dim, eigenvalues,
	   workspace, &workspace_size, &retval);
  if ( retval ) {
    cerr << "Dsyev  error: " << retval << "\n";
    throw(MLException("eigen decompfailed"));
  }

  double tsd = eigenvalues[0];
  Vec dir = Vec(A[0], A[1], A[2]);

  // Force dir N->C
  if (dot(dir, points[end] - points[bgn]) < 0) {
    dir *= -1.0;
  }

  // Brute force TSD
  double tsd_bf = 0.0;
  for (int p = bgn; p <= end; p++)
    tsd_bf += cross(points[p] - pos, dir).lensq();
  double dif = (tsd - tsd_bf) / (tsd + tsd_bf);
  if (abs(dif) > 0.1) {
    cerr << "Problem with TSD calcs" << endl;
    cerr << tsd << " " << tsd_bf << endl;
    //throw(MLException("problem with TSD calcs"));
  }

  line->setTSD(1, tsd);
  line->pos = pos;
  line->dir = dir;

  return line;
}

Line* LineFit::findSeg(int bgn, int end) throw(MLException) {
  Line* line = getCachedSeg(bgn, end);
  if (line == 0)
    line = setCachedSeg(bgn, end, calcSeg(bgn, end));
  return line;
}

Line* LineFit::findList(int bgn1, int n_lines) throw (MLException) {

  if (n_lines < 1)
    throw(MLException("findList(): n_lines < 1"));
  if (n_lines > n_lines_max)
    throw(MLException("findList(): n_lines > max"));
  if ((n_points - bgn1) < (n_lines * min_line_len)) {
    ostringstream oss;
    oss << "findList(): "
	<< "bgn: " << bgn1 << " n_pts: " << n_points << " n_lns: " << n_lines;
    throw(MLException(oss));
  }

  Line* best_lines = getCachedList(n_lines, bgn1);
  if (best_lines != 0)
    return best_lines;

  // Do (n_lines == 1) separately to
  // insure that (end == (n_points - 1))
  if (n_lines == 1) {
    best_lines = findSeg(bgn1, n_points - 1);
  }
  else {

    double best_tsd = 0.0;

    int bgn2 = bgn1 + min_line_len;
    int bgn2_max = n_points - (n_lines - 1) * min_line_len;
    while (bgn2 <= bgn2_max) {

      Line* lines =
	findSeg(bgn1, bgn2 - 1)->append(n_lines-1, findList(bgn2, n_lines-1));

      double tsd = lines->getTSD(n_lines);
      if ((best_lines == 0) || (tsd < best_tsd)) {
	best_tsd = tsd;
	best_lines = lines;
      }

      bgn2++;
    }

  }

  setCachedList(n_lines, bgn1, best_lines);
  return best_lines;
}

Line* LineFit::find(int* n_lines, int n_points, Vec points[])
  throw(MLException) {
  //double rmsd_max = 1.0;
  double rmsd_max = 10.0;
  double msd_max = rmsd_max * rmsd_max;
  double tsd_max = ((double) n_points) * msd_max;

  *n_lines = 1;
  for (;;) {
    Line* lines = findList(0, *n_lines);
    if ((*n_lines) >= n_lines_max) // Not going to do any better
      return lines;
    if ((lines->getTSD(*n_lines)) <= tsd_max)
      return lines;
    (*n_lines)++;
  }

  throw(MLException("LineFit::find(): shouldn't be able to get here"));
  return 0;
}
