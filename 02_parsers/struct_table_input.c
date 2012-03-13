# include "struct.h"

int read_integral_table (char * filename ) {

    int a, b;
    int line_ctr;
    double value, delta = 0.0;;
    char line[BUFFLEN];
    FILE * fptr;
    /* int_table defined in struct_table.c and declared in struct.h */
    memset (int_table[0], 0, TABLE_SIZE*TABLE_SIZE*sizeof(double) );

    fptr = efopen (filename, "r" ) ;
    if ( !fptr ) return 1;

    line_ctr = 0;
    memset (line,  0, BUFFLEN);
    while(fgets(line, BUFFLEN, fptr)!=NULL){
	line_ctr ++;
	if ( line[0] == '%' ) { /* comment */
	    if (strstr (line, "delta" ) ) {
		sscanf (line, " %*s  %*s   %lf", &delta );
	    }
	    continue;
	}
	if ( sscanf (line, "%d %d  %*s   %*s   %lf ",
		     &a, &b, &value ) == 3 ) {
	    if ( a >= TABLE_SIZE ||  b >= TABLE_SIZE  ) {
		fprintf (stderr, "Error reading in the lookup table:\n");
		fprintf (stderr, "line %d: %s\n", line_ctr, line);
		return 1;
	    }
	    int_table [a][b] = value;
	    int_table [b][a] = value;
	}
	memset (line,  0, BUFFLEN);
    }
    fclose (fptr);

    if ( !delta ) {
	fprintf (stderr, "Error reading in the lookup table:\n");
	fprintf (stderr, "no delta found\n");
	return 1;
    }

    options.alpha = delta;
    
    return 0;
    
}
