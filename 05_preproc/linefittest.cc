
#include <iostream>
#include <cmath>
using namespace std;

#include "mlexception.h"
#include "vec.h"
#include "linefit.h"

int main(int argc, char *argv[]) {

  int n_points = 5;
  Vec points[n_points];
  for (int p = 0; p < n_points; p++) {
    points[p] = Vec(0.1 * ((double) (p - 2)),
		    0.2 * ((double) (p - 2)),
		    0.3 * ((double) (p - 2)));
  }

  points[2] += Vec(0.1, 0.0, 0.0);

  try {

    LineFitCache linefitcache(n_points, points);
    LineFit linefit = linefitcache.get(0, 2);
    linefit.print();

  }
  catch (MLException& e) {
    cerr << e.what() << endl;
    return e.retval;
  }

  cout << "success" << endl;
  return 0;
}
