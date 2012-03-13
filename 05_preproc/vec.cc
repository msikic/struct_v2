
#include <cmath>

#include "vec.h"

Vec::Vec() {
  x = 0.0;
  y = 0.0;
  z = 0.0;
}

Vec::Vec(double x, double y, double z) {
  this->x = x;
  this->y = y;
  this->z = z;
}

Vec& Vec::operator+=(Vec a) {
  x += a.x;
  y += a.y;
  z += a.z;
  return *this;
}

Vec& Vec::operator-=(Vec a) {
  x -= a.x;
  y -= a.y;
  z -= a.z;
  return *this;
}

Vec& Vec::operator*=(double s) {
  x *= s;
  y *= s;
  z *= s;
  return *this;
}

Vec& Vec::operator/=(double s) {
  x /= s;
  y /= s;
  z /= s;
  return *this;
}

Vec& Vec::normalize() {
  double s = sqrt(x*x + y*y + z*z);
  x /= s;
  y /= s;
  z /= s;
  return *this;
}

double Vec::len() {
  return sqrt(x*x + y*y + z*z);
}

double Vec::lensq() {
  return (x*x + y*y + z*z);
}

Vec operator+(Vec a, Vec b) {
  Vec r = a;
  return r += b;
}

Vec operator-(Vec a, Vec b) {
  Vec r = a;
  return r -= b;
}

Vec operator*(double s, Vec v) {
  Vec r = v;
  return r *= s;
}

Vec operator*(Vec v, double s) {
  Vec r = v;
  return r *= s;
}

Vec operator/(Vec v, double s) {
  Vec r = v;
  return r /= s;
}

double dot(Vec a, Vec b) {
  return (a.x * b.x + a.y * b.y + a.z * b.z);
}

Vec cross(Vec a, Vec b) {
  return Vec(a.y * b.z - a.z * b.y,
	     a.z * b.x - a.x * b.z,
	     a.x * b.y - a.y * b.x);
}

Vec normalize(Vec a) {
  Vec r = a;
  return r.normalize();
}

double lensq(Vec& a) {
  return a.lensq();
}
