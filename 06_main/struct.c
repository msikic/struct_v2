#include "struct.h"
#include <stdlib.h>


Options options;

int set_default_options() {
    /* set the default options */
    memset(&options, 0, sizeof (Options));


    options.min_no_SSEs = 4;

    options.sheet_cosine = 1.0; /*min cos angle between two strands
				 in a sheet to be represented by the
				 same vector */
    options.alpha = 0.0; /* gaussian width for the scoring fn F */
    /* note that it is set from the input table */
    /* rather than from the cmd file        */
    options.tol /* tol in F for the claim that the query and  */
            /* the target are one and the same structure  */
            = 0.1;
    options.z_max_store /* z-score cutoff to store a map */
            = -2.0;
    options.z_max_out /* z-score cutoff for map output  */
            = -2.0;

    options.F_guess_max /* max  value of F for an intial guess    */
            = -1.5;
    options.F_eff_max /* max  value of F after gradient descent */
            = -3.0;
    options.z_max_corr /* max value of z for which two maps are */
            /* considered correlated */
            = -2;
    options.z_min_compl /* min value of z for which two maps are */
            /* considered complementary */
            = 2;
    options.grad_step_size/* step size in gradient descent */
            = 0.0001;
    options.grad_stop_tol /* stopping precision in gradient descent */
            = 1.e-4;

    options.far_far_away /* nonsense distance for the pairwise alignment*/
            = -10;
    options.gap_open /* gap penalties for the pairwise alignment */
            = 0.0;
    options.gap_extend
            = 0.0;
    options.endgap
            = 0.0;
    options.far_away_cosine /* minimum cosine for F_effective estimate*/
            = 0.8;
    options.grid_size /* minimal number of points for the sphere grid */
            = 400;
    options.number_maps_cpl /* number of top scoring maps among which to   */
            /* look for a complement */
            = 10;
    options.number_maps_out /* number of top scoring maps to output     */
            = 5;
    options.grad_max_step /* max number of steps in  gradient descent */
            = 100;
    options.exp_table_size /* size of the lookup table for the exp funcion */
            = TABLE_SIZE; /* note: changing exp table size not implemented */
    options.smith_waterman = 1;
    options.use_perp = 0;

    options.verbose = 0;
    options.print_header = 0;
    options.report_no_sse_overlap = 1;

    /* penalizing the length mismatch */
    options.use_length = 0;
    options.H_length_mismatch_tol = 10.0;
    options.S_length_mismatch_tol = 5.0;

    /* path to the integral table */
    memset(options.path, 0, BUFFLEN);
    sprintf(options.path, "%s", INTEGRAL_TABLE);
    
    options.sw_score = 1; // database search calculation type
    options.score_out = 2; // 0 - SW, 1 - Hungarian, 2 - both - choose better one
    options.threshold_distance = 25; 
    

    return 0;
}

/***************************************************/
/***************************************************/
/***************************************************/

/***************************************************/
int main(int argc, char * argv[]) {

    Descr qry_descr = {
        {0}
    };
    Descr tgt_descr = {
        {0}
    };
    clock_t CPU_time_begin, CPU_time_end;
    int retval, qry_done, tgt_done;
    int db_ctr, db_effective_ctr;
    int user_defined_name;
    FILE * qry_fptr = NULL, * tgt_fptr = NULL, * digest = NULL;

//    Score score;


    //int compare(Descr *descr1, Descr *descr2, Score * score);
    int compare(Descr *descr1, Descr *descr2, Score * score, Score * score_hung);
    int read_cmd_file(char *filename);

    if (argc < 3) {
        fprintf(stderr,
                "Usage: %s <db file> <qry file> [<parameter file>].\n",
                argv[0]);
        exit(1);
    }
    if (!(qry_fptr = efopen(argv[2], "r"))) return 1;
    if (!(tgt_fptr = efopen(argv[1], "r"))) return 1;

    /* set defaults: */
    set_default_options();

    /* change them with the cmd file, if the cmd file given */
    if (argc == 4) {
        if (read_cmd_file(argv[3])) return 1;
    }

    /* read in the table of integral values */
    /* the array int_table in struct_table.c */
    if (read_integral_table(options.path)) {
        fprintf(stderr, "In data file  %s.\n\n", options.path);
        exit(1);
    }
    set_up_exp_table();

    user_defined_name = options.outname[0];


    /*********************************/
    /* loop over the query database :*/
    qry_done = 0;
    retval = -1;
    db_effective_ctr = 0;
    CPU_time_begin = clock();

    while (!qry_done) {
        retval = get_next_descr(qry_fptr, &qry_descr);
        if (retval == 1) {
            continue;
        } else if (retval == -1) {
            qry_done = 1;
            continue;
        }

        /* digest file for larger scale comparisons */
        if (!digest) {
            if (!user_defined_name) {
                sprintf(options.outname, "%s.struct_out",
                        qry_descr.name);
            }

            // ************ added by Mile
            // output name in postprocessing consists of query and target name
            retval = get_next_descr(tgt_fptr, &tgt_descr);
            if (options.postprocess) {
                sprintf(options.outname, "%s_%s.struct_out", qry_descr.name, tgt_descr.name);
            }


            // ************* end by Mile

            digest = efopen(options.outname, "w");
            if (!digest) exit(1);
            if (options.print_header) {
                fprintf(digest, "%% columns: \n");
                fprintf(digest, "%% query, target: structure names\n");
                fprintf(digest, "%% geom_z:  z score for the orientational match \n");
                fprintf(digest, "%% <dL>:    average length mismatch for matched SSEs \n");
                fprintf(digest, "%% T:       total score assigned to matched SSEs \n");
                fprintf(digest, "%% frac:    T divided by the number of matched SSEs \n");
                fprintf(digest, "%% GC_rmsd: RMSD btw geometric centers of matched SSEs (before postprocessing) \n");
                fprintf(digest, "%% A:       (after postprocessing) the alignment score \n");
                fprintf(digest, "%% aln_L:   (after postprocessing) the alignment length \n\n");
                fprintf(digest, "%% %6s%6s %6s %6s  %6s %6s %6s %6s %6s %6s \n",
                        "query ", "target ", "geom_z", "<dL>", "  T  ", "frac",
                        "GC_rmsd", "rmsd  ", "A  ", "aln_L  ");
            }

        } else {
            /* otherwise write to the same old digest file */
        }

        /* loop over the database :*/


        // Added by Mile - using FOR instead of WHILE - parallelization

        int tgt_counter = 0;
        int i;
        int *retval_array;
        Descr *tgt_descr_array;
         
        rewind(tgt_fptr);
        tgt_done = 0;

        /*
         * Counting number of successful targets
         */
        while (!tgt_done) {
            retval = get_next_descr(tgt_fptr, &tgt_descr);
            if (retval == 0 || retval == 1) {
                tgt_counter++;
            } else if (retval == -1) {
                tgt_done = 1;
            }

        }
        /*
         * Initialization of a Descr array (array of targets) - easy parallelization
         */
        
        rewind(tgt_fptr);
        tgt_descr_array = (Descr *) calloc(tgt_counter, sizeof(Descr));
        if (tgt_descr_array == NULL) {
            printf("malloc return NULL!\n");
        }
        retval_array = (int *) calloc(tgt_counter, sizeof(int));    
        if (retval_array == NULL) {
            printf("malloc return NULL!\n");
        }
        
        
        /*
         * Storing targets a returning values
         */
        
        for(i = 0; i < tgt_counter; ++i) {
            retval = get_next_descr(tgt_fptr, &tgt_descr_array[i]);
            retval_array[i] = retval;
        }
        
        // Added by Mile - end


        rewind(tgt_fptr);
        //	tgt_done = 0;
        db_ctr = 0;
        db_effective_ctr = 0;
        if (!user_defined_name) CPU_time_begin = clock();
        retval = -1;

        /*
                while ( ! tgt_done) {
         */
        
        // Start of parallelization
        if (options.postprocess) omp_set_num_threads(1);
        else omp_set_num_threads(6);
        
        #pragma omp parallel // num_threads(1)
        {

            #pragma omp for
            for (i = 0; i < tgt_counter; ++i) { // Added by Mile

                
                int retval = retval_array[i];
                /*
                 * Two scores one for Smith Waterman, another for Hungarian in database search phase
                 */
                Score score;
                Score score_hung;
                Descr tgt_descr = tgt_descr_array[i];
/*
                printf("%s %d\n", tgt_descr.name, retval);
*/
                
 //               Descr qry_descr = qry_descr;
                
                
                #pragma omp atomic
                db_ctr++; // atomic
                
/*
                retval = get_next_descr(tgt_fptr, &tgt_descr);
*/
                if (retval == 1) {
                    continue;
                } else if (retval == -1) {
                    //  tgt_done = 1;
                    printf("Error!!!!\n");
                    exit(1);
                    // added by Mile
                } else {

                    /* min number of elements */
                    int helix_overlap =
                            (qry_descr.no_of_helices < tgt_descr.no_of_helices) ?
                            qry_descr.no_of_helices : tgt_descr.no_of_helices;
                    int strand_overlap =
                            (qry_descr.no_of_strands < tgt_descr.no_of_strands) ?
                            qry_descr.no_of_strands : tgt_descr.no_of_strands;
                    double fraction_assigned;
                    int query_size = qry_descr.no_of_strands + qry_descr.no_of_helices;
                    int target_size = tgt_descr.no_of_strands + tgt_descr.no_of_helices;
                    if (helix_overlap + strand_overlap >= options.min_no_SSEs) {

                        #pragma omp atomic
                        db_effective_ctr++; // atomic

                        /* here is the core of the operation: */
                        retval = compare(&tgt_descr, &qry_descr, &score, &score_hung);
                        if (retval) {
                            printf(" error comparing  db:%s   query:%s   \n",
                                    tgt_descr.name, qry_descr.name);
                            exit(retval);
                        }

 
                        /*
                         * Output score. Can be based:
                         * - only on SW alignment during the database search
                         * - only on Hungarian algorithm during the database search
                         * - on combination depending on the postprocessing score
                         */  
                                
                        
                        switch (options.score_out) {
                            case 0: // SW
                                if (query_size > target_size) {
                                   fraction_assigned = score.total_assigned_score / target_size;
                                } else {
                                   fraction_assigned = score.total_assigned_score / query_size;
                                }
                                retval =  print_score(digest, &qry_descr, &tgt_descr, &score, fraction_assigned, 1);  
                                break;
                            case 1: // Hungarian
                                if (query_size > target_size) {
                                    fraction_assigned = score_hung.total_assigned_score / target_size;
                                } else {
                                    fraction_assigned = score_hung.total_assigned_score / query_size;
                                }
                                retval =  print_score(digest, &qry_descr, &tgt_descr, &score_hung, fraction_assigned, 1);
                                break;
                            case 2: // either SW or Hungarian depends on score
                                if (score.total_assigned_score > score_hung.total_assigned_score) {
                                    if (query_size > target_size) {
                                        fraction_assigned = score.total_assigned_score / target_size;
                                    } else {
                                        fraction_assigned = score.total_assigned_score / query_size;
                                    }
                                    retval =  print_score(digest, &qry_descr, &tgt_descr, &score, fraction_assigned, 1);
                                } else {
                                    if (query_size > target_size) {
                                        fraction_assigned = score_hung.total_assigned_score / target_size;
                                    } else {
                                        fraction_assigned = score_hung.total_assigned_score / query_size;
                                    }
                                    retval =  print_score(digest, &qry_descr, &tgt_descr, &score_hung, fraction_assigned, 1);
                                    
                                }
                                break;
                        }
                            
                            
                        
                        if (retval) {
                            printf("error in printing to output file\n");
                            exit(retval);
                        }

                    } else if (options.report_no_sse_overlap) {
                        retval =  print_score(digest, &qry_descr, &tgt_descr, &score, fraction_assigned, 0);
                        if (retval) {
                            printf("error in printing to output file\n");
                            exit(retval);
                        }
                    }
                }
                /*
                if (options.postprocess) tgt_done = 1; // for now, we postprocess only
                                                    one pair of structures (not structure against database) */
               // if (options.postprocess) break; // added by Mile tricky but I think it should work even without it
            }

        }

        // Added by Mile
        // Memory cleaning
        for(i = 0; i < tgt_counter; ++i) {
            descr_shutdown ( &tgt_descr_array[i] );
        }
        free(tgt_descr_array);
        free(retval_array);    
        
        
        // End added by Mile
        
        if (!user_defined_name && db_effective_ctr) {
            CPU_time_end = clock();
            fprintf(digest, "done   CPU:  %10.3lf s\n", (double) (CPU_time_end - CPU_time_begin) / CLOCKS_PER_SEC);
            fflush(digest);
        }

        if (!user_defined_name) {
            fclose(digest);
            digest = NULL;
        } /* otherwise we keep writing into the saem digest file */

        if (options.postprocess) qry_done = 1; /* for now, we postprocess only
						one pair of structures (not structure against database) */

    }

    if (digest) {
        CPU_time_end = clock();
        fprintf(digest, "done   CPU:  %10.3lf s\n", (double) (CPU_time_end - CPU_time_begin) / CLOCKS_PER_SEC);
        fflush(digest);
    }
    if (options.verbose) {
        printf("\n\nlooked at %d db entries.\n",
                db_effective_ctr);
        printf("the output written to %s.\n\n", options.outname);
    }
    /**************************************************/
    /* housekeeping, good for tracking memory leaks   */ if (digest) fclose(digest);
//    map_consistence(0, 0, NULL, NULL, NULL, NULL, NULL);
//    compare(NULL, NULL, NULL);
    descr_shutdown(&qry_descr);
    descr_shutdown(&tgt_descr);

    fclose(qry_fptr);
    fclose(tgt_fptr);

    return 0;

}

/**************************************************************************************/
/**************************************************************************************/

/**************************************************************************************/
int compare(Descr * descr1, Descr *descr2, Score *score, Score *score_hung) {
    /* descr 1 is db or target, descr2 is query - important in postprocessing*/
    int retval, retval1, retval2;
    int map_ctr, map_ctr_hung, best_ctr;
    int NX, NY; //NX_eff, NY_eff;
    Representation X_rep = {0}, Y_rep = {0};

    //static Map * map = NULL;
    Map * map = NULL;
    //static int map_max = MAP_MAX * 9;
    int map_max = MAP_MAX * 9;
    //static int best_max = MAP_MAX;
    int best_max = MAP_MAX;
    // static int *map_best = NULL;
    int *map_best = NULL;
    // static int NX_allocated = 0, NY_allocated = 0;
    
    /* TODO: go back to map lineage description */

    int construct_representation(double ** x, int *x_type, int *x_length,
            int * NX_effective, Descr * descr,
            int ** represents, int * is_rep_by);
    int match_clustering(Representation* X_rep, Representation* Y_rep,
            Map * map, int map_max,
            int * map_ctr, int * map_best, int best_max, int parent_map);
    int recursive_map_out(Map * map, int map_ctr,
            Descr * descr1, Descr * descr2,
            Protein *protein1, Protein *protein2,
            int depth);
    int rec_map_out_for_postproc(Map * map, int map_ctr,
            Representation *X_rep, Representation *Y_rep, int depth);


    /****************************************/
    /*  arrays for immediate consumption    */
    /****************************************/
    /* needs to come first to figure out NX */
    rep_initialize(&X_rep, descr1);
    rep_initialize(&Y_rep, descr2);


    /****************************************/
    /*  initialization                      */
    /****************************************/
    /*  shorthands   */
    NX = X_rep.N_full;
  //  NX_eff = X_rep.N_compact;
    NY = Y_rep.N_full;
//    NY_eff = Y_rep.N_compact;


    if ((map = init_map(NX, NY, map_max)) == NULL) return 1;
    
    map_best = emalloc(map_max * sizeof (int));
    if (!map_best) return 1;
    


    /****************************************/
    /*  look for complementary maps         */
    /****************************************/
    map_ctr = 0;
    map_ctr_hung = 0;
    retval = 1;
    
    Map * map_hung;
    if ((map_hung = init_map(NX, NY, map_max)) == NULL) return 1;
    int *map_best_hung = NULL;
    map_best_hung =  emalloc(map_max * sizeof (int));
    
    /*
     * Structure matching based on secondary structure orientation
     * Two algorithms used:
     * - Smith Waterman
     * - Hungarian
     */
    
    retval1 = complement_match(&X_rep, &Y_rep, map, map_max,
                &map_ctr, map_best, best_max, -1);
    options.sw_score = 0;
    retval2 = complement_match(&X_rep, &Y_rep, map_hung, map_max,
                &map_ctr_hung, map_best_hung, best_max, -1);
    
    
    if (retval1 || retval2) {
        rep_shutdown(&X_rep);
        rep_shutdown(&Y_rep);
    // added by Mile
        clean_map(map, map_max);
        free(map_best);
        
        clean_map(map_hung, map_max);
        free(map_best_hung);
        
    //added by Mile
        
        
        return retval1;
    }

    /****************************************/
    /*   output the best case stats         */
    /****************************************/
    memset(score, 0, sizeof (Score));
    memset(score_hung, 0, sizeof (Score));
    
    best_ctr = 0;
    while (map_best[best_ctr] > -1) {
        best_ctr++;
    } // TODO - check cases when hungarian returns score and SW does not

    if (best_ctr) {
        assign_score(map, map_best, &X_rep, &Y_rep, score);
        assign_score(map_hung, map_best_hung, &X_rep, &Y_rep, score_hung);

        /****************************************/
        /*     postprocess - find the actual tf */
        /* 	   & mapping on the bb level    */
        /****************************************/
        if (options.postprocess) {
            Protein qry_structure = {0};
            Protein tgt_structure = {0};
            /* input CA coordinates for the vectors
               -- for now the input should be a little different
               if we are expecting to postprocess -- we will use the
               original pdb */
            /* the last arg is 1 in postprocessing */
            retval = read_pdb(options.pdbf_qry, options.chain_qry, &qry_structure, 1);
            if (retval) {
                fprintf(stderr, "Error reading %s, chain %c; retval %d\n",
                        options.pdbf_qry, options.chain_qry, retval);
                return 1;
            }
            retval = read_pdb(options.pdbf_tgt, options.chain_tgt, &tgt_structure, 1);
            if (retval) {
                fprintf(stderr, "Error reading %s, chain %c; retval %d\n",
                        options.pdbf_tgt, options.chain_tgt, retval);
                return 1;
            }
            /* for now,  we will just postprocess the  best map */

            best_ctr = 0;
            map_ctr = map_best[best_ctr];
            int map_ctr_hung = map_best[best_ctr];

            postprocess(descr1, &tgt_structure, &X_rep, descr2, &qry_structure, 
                    &Y_rep, map + map_ctr, score);
            postprocess(descr1, &tgt_structure, &X_rep, descr2, &qry_structure, 
                    &Y_rep, map_hung + map_ctr_hung, score_hung);
            
            
            rec_map_out_for_postproc(map, map_ctr, &X_rep, &Y_rep, 0); 
            rec_map_out_for_postproc(map_hung, map_ctr_hung, &X_rep, &Y_rep, 0);
            
            if(score_hung->total_assigned_score > score->total_assigned_score) {
                options.score_out = 1;
            }

            if (options.verbose) {
                if (options.score_out == 0){
                    recursive_map_out(map, map_ctr, descr1, descr2,
                                &qry_structure, & tgt_structure, 0);
                }
                else {
                    recursive_map_out(map_hung, map_ctr_hung, descr1, descr2,
                                &qry_structure, & tgt_structure, 0);
                } 
            }
            protein_shutdown(&qry_structure);
            protein_shutdown(&tgt_structure);

        } else {

            /****************************************/
            /*     output the maps                  */
            /****************************************/
            if (options.verbose) {
                best_ctr = 0;
                while (best_ctr < options.number_maps_out && best_ctr < map_max
                        && (map_ctr = map_best[best_ctr]) > -1
                        && map[map_best[best_ctr]].z_score < options.z_max_out) {
                    recursive_map_out(map, map_ctr, descr1, descr2, NULL, NULL, 0);
                    best_ctr++;
                }
            }
            // TODO - do something for Hungarian too. Currently without maps
        }
    } /* end if best_ctr */

    /****************************************/
    /*  shutting down                       */
    /****************************************/
    rep_shutdown(&X_rep);
    rep_shutdown(&Y_rep);

    // added by Mile
    free(map_best);
    clean_map(map, map_max); 
    
    clean_map(map_hung, map_max);
    free(map_best_hung);
    //added by Mile
        
    return 0;

}

/**************************************************************************/
double quat_rmsd(double parent_map_q[4], double current_map_q[4]) {

    double rmsd = 0.0;
    double aux;
    int i;

    for (i = 0; i < 4; i++) {
        aux = parent_map_q[i] - current_map_q[i];
        rmsd += aux*aux;
    }
    rmsd /= 4;

    return sqrt(rmsd);
}

/**
 * 
 * @param map
 * @param map_best
 * @param X_rep
 * @param Y_rep
 * @param score
 * @return 
 */
int assign_score(Map *map, int *map_best, Representation *X_rep, Representation *Y_rep, Score *score) {
    int sub_map_ctr, i;
    int NX = X_rep->N_full;
    Map * current_map;
    score -> z_score = map[map_best[0]].z_score;
    score -> rmsd = map[map_best[0]].rmsd;
    

    score -> total_assigned_score = map[map_best[0]].assigned_score;

    if (score -> total_assigned_score > 0) {
        double avg_length_mismatch = 0;
        int j, map_size = 0;
        for (i = 0; i < NX; i++) {
            j = map[map_best[0]].x2y[i];
            if (j < 0) continue;
            avg_length_mismatch += fabs(X_rep->length[i] - Y_rep->length[j]);
            map_size++;
        }
        if (map_size) avg_length_mismatch /= map_size;
        score -> avg_length_mismatch = avg_length_mismatch;

    } else {
        score -> avg_length_mismatch = 0;
    }
    score -> number_of_maps = 1;

    current_map = map + map_best[0];
    score -> w_submap_z_score = map[map_best[0]].z_score;
    sub_map_ctr = 0;
    while (current_map->submatch_best && current_map->submatch_best[0] > -1) {

        current_map = map + current_map->submatch_best[0];
        if (current_map->z_score >= 0) continue;
        score -> total_assigned_score += current_map->assigned_score;
        score -> w_submap_z_score *= current_map->z_score;
        score -> number_of_maps++;
        if (sub_map_ctr < 3) {
            score->q_rmsd[sub_map_ctr] = quat_rmsd((map + map_best[0])->q, current_map->q);
        }

        sub_map_ctr++;
    }

    score -> w_submap_z_score = fabs(score -> w_submap_z_score);
    return 0;
}

Map * init_map(int NX, int NY, int map_max) {
    int map_ctr;
    Map *map;
    map = emalloc(map_max * sizeof (Map)); /* TODO is this enough? */
    if (!map) return NULL;
    for (map_ctr = 0; map_ctr < map_max; map_ctr++) {
        if (initialize_map(map + map_ctr, NX, NY)) return NULL;
    }

    for (map_ctr = 0; map_ctr < map_max; map_ctr++) {
        (map + map_ctr)->x2y_size = NX;
        (map + map_ctr)->y2x_size = NY;
    }
    return map;
}


/**
 * 
 * @param map
 * @param map_max
 * @return 
 */

int clean_map(Map *map, int map_max) {
    int map_ctr;
    for (map_ctr = 0; map_ctr < map_max; map_ctr++) {
        if (free_map(map + map_ctr)) return 1;
    }
    free(map);
    return 0;
}
