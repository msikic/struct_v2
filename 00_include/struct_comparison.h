# include "struct.h"

# define  ERROR            0
# define  NO_MATCH         1
# define  MATCH_TO_PARENT  2
# define  FULL_MATCH       3
# define  PIECEWISE_MATCH  4
# define  PARTIAL_MATCH    5

# define  MATCH_TABLE_SIZE 6


extern int  match_algebra[][];
extern int  match_algebra_set;

typedef struct Comparison_node {

    struct Comparison_node *left[2], *right[2], *parent;
    /* two versions of left and right (to handle swaps)
       - a bifurcating binary tree - what a beaut */ 
    int type;
    int id;
    int match_flag;
    Node *node1, * node2; /* pointers to the trees being compared */
    
    double ** matrix1, ** matrix2;
    int       size1, size2;
    double  rot[3][3]; /* TODO do I really nead an orthogonal tfm here? */
    double  transl[3];
    double rmsd[2];
    
} Comparison_node;


typedef struct Comparison_tree{

    Comparison_node *root;
    Comparison_node *current_node; /* to be used while building the tree */
    int size; /* leaves + inner nodes */
    int alloc_size; /* mem space allocated at root */
    int no_of_leaves;
    int bifurcates;
    int next_available_node;
    double **matrix1, ** matrix2;

} Comparison_tree;



