
class Vec {
 public:
  double x;
  double y;
  double z;
  Vec();
  Vec(double x, double y, double z);
  Vec& operator+=(Vec a);
  Vec& operator-=(Vec a);
  Vec& operator*=(double s);
  Vec& operator/=(double s);
  Vec& normalize();
  double len();
  double lensq();
};

Vec operator+(Vec a, Vec b);
Vec operator-(Vec a, Vec b);
Vec operator*(double s, Vec v);
Vec operator*(Vec v, double s);
Vec operator/(Vec v, double s);
double dot(Vec a, Vec b);
Vec cross(Vec a, Vec b);
Vec normalize(Vec a);
double lensq(Vec& a);
