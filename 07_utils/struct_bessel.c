#include <struct.h>

 /*----------------------------------------------------------------------------*/
 /*-General Information--------------------------------------------------------*/
 /*                                                                            */
 /*   This function computes the modified Bessel Function of the first         */
 /*   kind and zero order. Taken from the Numerical Recipes book.              */
 /*                                                                            */
 /*----------------------------------------------------------------------------*/
 /*-Interface Information------------------------------------------------------*/
 double mBesselF (double x  )   {   /*  I   Argument to compute Io(x).         */
 /*----------------------------------------------------------------------------*/
   double ax, bessel;
   double y;
 
   if ((ax = fabs (x)) < 3.75)
     {
       y = x / 3.75;
       y *= y;
       bessel = 1.0 + y * (3.5156229 + y * (3.0899424 + y * (1.2067492
                                                             + y * (0.2659732 +
                                                                    y *
                                                                    (0.360768e-1
                                                                     +
                                                                     y *
                                                                     0.45813e-2)))));
     }
   else
     {
       y = 3.75 / ax;
       bessel = (exp (ax) / sqrt (ax)) * (0.39894228 + y * (0.1328592e-1
                                                            +
                                                            y * (0.225319e-2 +
                                                                 y *
                                                                 (-0.157565e-2
                                                                  +
                                                                  y *
                                                                  (0.916281e-2
                                                                   +
                                                                   y *
                                                                   (-0.2057706e-1
                                                                    +
                                                                    y *
                                                                    (0.2635537e-1
                                                                     +
                                                                     y *
                                                                     (-0.1647633e-1
                                                                      +
                                                                      y *
                                                                      0.392377e-2))))))));
     }
 
   return (bessel);
 }
