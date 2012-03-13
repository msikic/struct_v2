

float Torsion(float *Coord1, float *Coord2, float *Coord3, float *Coord4);
float Dist(float *Coord1, float *Coord2);
float Ang(float *Coord1, float *Coord2, float *Coord3);
void PHI(CHAIN *Chain, int Res);
void PSI(CHAIN *Chain, int Res);
void OMEGA(CHAIN *Chain, int Res);
void Place123_X(float *Coord1, float *Coord2, float *Coord3,
		float Dist3X, float Ang23X, float *CoordX);
float VectorProduct(float *Vector1, float *Vector2, float *Product);
void Project4_123(float *Coord1, float *Coord2, float *Coord3, float *Coord4, 
		  float *Coord_Proj4_123);
double GetAtomRadius(char *AtomType);
