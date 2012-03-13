# ifndef _UTILS_H
# define _UTILS_H
# include <stdio.h>


/******************************/
/*   tokenizer                */
/******************************/

# define TOK_TOOMNY  1 /* tokenazier error codes */
# define TOK_TOOLONG 2
# define MAX_TOK 30  /* max number of tokens per line in the commandfile */
# define LONGSTRING  250
# define MEDSTRING   100
# define SHORTSTRING  25

int  cluster_counter (int  no_of_things,  int **neighbors, 
		      int * no_of_clusters, int ** clusters);
int tokenize ( char  token[MAX_TOK][MEDSTRING], int * max_token,
	       char * line, char comment_char);




int      array_qsort (int * sorted_pos, double * sa, int sequence_length );
char   **chmatrix(int rows, int columns);
double **dmatrix(int rows, int columns);
double ***d3matrix(int rows, int columns, int floors);
void *   emalloc(int	size);
FILE *   efopen(char * name, char * mode);
void     free_matrix(void **m);
void     free_cmatrix(char **m);
void     free_imatrix(int **m);
void     free_dmatrix(double **m);
void     free_d3matrix(double  ***m);
int    **intmatrix(int rows, int columns);
int intmatrix_init(int **matrix, int rows, int columns, int val);
int      string_clean ( char* string, int length);

# endif
