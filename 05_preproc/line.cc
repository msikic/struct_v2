
#include <iostream>
#include <sstream>
using namespace std;

#include "mlexception.h"
#include "vec.h"
#include "line.h"


Line::Line() {
  n_lines_max = 0;
  bgn = 0;
  end = 0;
  pos = Vec(0.0, 0.0, 0.0);
  dir = Vec(0.0, 0.0, 0.0);
  next = new Line*[0];
  tsd = new double[0];
}

Line::Line(int n_lines_max, int bgn, int end) {
  this->n_lines_max = n_lines_max;
  this->bgn = bgn;
  this->end = end;

  pos = Vec(0.0, 0.0, 0.0);
  dir = Vec(0.0, 0.0, 0.0);

  next = new Line*[n_lines_max];
  tsd = new double[n_lines_max];
  for (int nl = 0; nl < n_lines_max; nl++) {
    next[nl] = (Line*) 0;
    tsd[nl] = -1.0;
  }
}

Line::~Line() {
  delete[] next;
  delete[] tsd;
}

double Line::getTSD(int n_lines) {
  return tsd[n_lines - 1];
}

Line* Line::getNext(int n_lines) {
  return next[n_lines - 1];
}

void Line::setTSD(int n_lines, double val) {
  tsd[n_lines - 1] = val;
}

void Line::setNext(int n_lines, Line* lines) {
  next[n_lines - 1] = lines;
}

Line* Line::append(int n_lines, Line* lines) throw(MLException) {
  if (lines == 0) {
    ostringstream oss;
    oss << "Line::append(): lines == null: n_lines == " << n_lines;
    throw(MLException(oss));
  }
  if (n_lines < 1)
    throw(MLException("Line.append(): n_lines < 1")); 
  if ((n_lines + 1) > n_lines_max)
    throw(MLException("Line.append(): (n_lines + 1) > n_lines_max")); 
  setTSD(n_lines + 1, getTSD(1) + lines->getTSD(n_lines));
  setNext(n_lines + 1, lines);
  return this;
}

void Line::print() {
  cout << "tsd: " << tsd << endl;
  cout << "pos: " << pos.x << " " << pos.y << " " << pos.z << endl;
  cout << "dir: " << dir.x << " " << dir.y << " " << dir.z << endl;
}

void copyToArray(Line* array, int n_lines, Line* lines) {
  int a = 0;
  int nl = n_lines;
  while (nl > 0) {
    if (lines == 0) {
      ostringstream oss;
      oss << "copyToArray(): premature null: "
	  << "n_lns: " << n_lines << ", a: " << a ;
      throw(MLException(oss));
    }
    array[a] = *lines;
    a++;
    lines = lines->getNext(nl);
    nl--;

  }
}

//void free(int n_lines, Line* lines) {
//  while (n_lines > 0) {
//    Line* tmp = lines;
//    lines = lines->getNext(n_lines);
//    n_lines--;
//    delete tmp;
//  }
//}
//
//
//Line::Line(const Line& line) {
//  n_lines_max = line.n_lines_max;
//
//  rmsd_seg = line.rmsd_seg;
//  pos = line.pos;
//  dir = line.dir;
  
//  next = new Line[n_lines_max];
//  rmsd_tot = new double[n_lines_max];
//  for (int nl = 0; nl < n_lines_max; nl++) {
//    next[nl] = line.next[nl];
//    rmsd_tot[nl] = line.rmsd_tot[nl];
//  }
//}
//
//Line& Line::operator=(const Line& other) {
//  if (this != &other) {
//    n_lines_max = other.n_lines_max;
//
//    rmsd_seg = other.rmsd_seg;
//    pos = other.pos;
//    dir = other.dir;
//
//    delete[] next;
//    delete[] rmsd_tot;
//    next = new Line[n_lines_max];
//    rmsd_tot = new double[n_lines_max];
//    for (int nl = 0; nl < n_lines_max; nl++) {
//      next[nl] = other.next[nl];
//      rmsd_tot[nl] = other.rmsd_tot[nl];
//    }
//  }
//}
