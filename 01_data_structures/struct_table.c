# include "struct.h"


double int_table [TABLE_SIZE][TABLE_SIZE];

double exp_table [TABLE_SIZE];


double lookup ( double cos_alpha, double cos_beta ) {

    int a, b;

    a = (int) ( (1.0 - cos_alpha)*NR_POINTS/2 );
    b = (int) ( (1.0 - cos_beta)*NR_POINTS/2 );

    if (a > NR_POINTS) a = NR_POINTS; 
    if (a < 0) a = 0; 
    if (b > NR_POINTS) b = NR_POINTS; 
    if (b < 0) b = 0;
    
    return int_table[a][b];
}

int set_up_exp_table () {
    int i;
    for (i=0; i< TABLE_SIZE; i++)
	exp_table[i] = exp (-(double)i*MAX_EXP_VALUE/(TABLE_SIZE-1));

    return 0;
}
