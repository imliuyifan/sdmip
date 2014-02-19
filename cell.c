/************************************************************************\
**
 ** cell.c
 **
 **
 **   solve_cell()
 **   new_cell()
 **   free_cell()
 **   write_cell()
 **
 ** History:
 **   18 Feb 1992 - <Jason Mai>  - created.
 **   20 Mar 1992 - <Jason Mai>  - created.
 **
 \************************************************************************/

#include <time.h>
#include "prob.h"
#include "cell.h"
#include "soln.h"
#include "quad.h"
#include "utility.h"
#include "theta.h"
#include "subprob.h"
#include "sigma.h"
#include "testout.h"
#include "solver.h"
#include "improve.h"
#include "optimal.h"
#include "memory.h"
#include "omega.h"
#include "lambda.h"
#include "log.h"
#include "master.h"
#include "cuts.h"
#include "batch.h"
#include <limits.h> 

FILE *fix;
int compare_function(const void  *a,const void *b) {
    u_type *x = (u_type *) a;
    u_type *y = (u_type *) b;
    if (x->u < y->u) return -1;
    else if (x->u > y->u) return 1; return 0;
}

/************************************************************************\
** This function represents the SD algorithm, as solved for a
 ** single cell.  It creates temporary data structures required for
 ** the solution, then proceeds with the iterative algorithm.
 ** It terminates when an optimal solution has been discovered,
 ** or when the limit on the number of iterations has been reached.
 ** It assumes that the cell passed to it has already been allocated
 ** and initialized, and likewise that the problem has been initialized.
 \************************************************************************/
void solve_cell(sdglobal_type* sd_global, cell_type *cell, prob_type *prob,
		vector x_k, char *fname)
{
	soln_type *soln;
	int omeg_idx;
	int i, j, cnt;
	double obj;
	BOOL new_omega;
	double phi[4]; /*  SS messing around  */
	double conf_int[2]; /* JH 3/13/98 */
	int num_subproblems = 0; /* JH 3/12/98 */
	FILE *f2_out; /* JH 3/12/98 */
    FILE *time_sample; /*Yifan 2013.06.30*/
    FILE *time_benders, *soln_benders; /* modified by Yifan 2013.10.06 */
	FILE *master_dual = NULL; /*Yiifan 2012-09-11*/
	FILE *batch_dual = NULL; /*Yifan 2012-09-11*/
	FILE *master_obj = NULL; /*Yifan 2012-09-23*/
	FILE *batch_obj = NULL; /*Yifan 2012-09-23*/
	FILE *batch_d = NULL; /* modified by Yifan 2012.09.23 */
	/*  FILE          *junk; */
#ifdef DEBUG
	struct tms total_init_time, total_term_time; /* zl 08/18/04 */
#endif
	clock_t total_start_time, total_end_time; /* zl 06/29/04 */
	clock_t iter_start_time, iter_end_time; /* zl 06/29/04 */
	clock_t argmax_start, argmax_end; /* zl 06/30/04 */
	FILE *time_file, *obj_file; /* zl 06/29/04 */
	char time_fname[NAME_SIZE * 2], obj_fname[NAME_SIZE * 2]; /* modified by Yifan to store longer names */
	/* zl 06/30/04 */
	BOOL incumb_change = FALSE;
    time_t start, stop;
    double diff, total_acc=0;
//  int batch_id;
	/*  double        sigma = 1.0;  _sigma_ here is a parameter to control the
	 effect of the quadratic term in regularized 
	 QP method. zl */
    

#ifdef DEBUG
	int idx;
#endif

#ifdef TRACE
	printf("Inside solve_cel\n");
#endif

	/* Open the file to record run time spent on different procedures.
	 zl, 06/29/04. */
	strcpy(time_fname, fname);
	strcat(time_fname, ".time.out");
	time_file = fopen(time_fname, "w+");

	/* open the file to record the UB, LB, and error per iteration. 
	 zl, 06/30/04. */
	strcpy(obj_fname, fname);
	strcat(obj_fname, ".obj.out");
	obj_file = fopen(obj_fname, "w+");
	fprintf(obj_file, "Iter \tincumb_est \tcandid_est \terror\n");
#ifdef DEBUG
	times(&total_init_time); /* zl, 08/18/04. */
#endif
	total_start_time = clock(); /* zl, 06/29/04. */
    start = time(NULL);

	phi[2] = 0;
	soln = new_soln(sd_global, prob, x_k);
	if (!(cell->master = new_master(prob->master, cell->cuts,
			prob->num->max_cuts, x_k)))
		err_msg("Copy", "solve_cell", "cell->master");

	/* Initialize data structure for Batch-Mean problem at the begining of the first iteration Yifan 2012-09-10 */
	if (prob->current_batch_id == 0)
	{
		if (!(sd_global->batch_problem = new_batch_problem(prob->master,
				prob->num->max_cuts)))
			err_msg("Copy", "solve_cell", "cell->subprob");
		new_batch_incumb(sd_global, prob, x_k);

		/* modified by Yifan 2013.02.15 */
		sd_global->bcuts = new_bcuts(prob, BATCH_SIZE, sd_global->bcuts);
		/* modified by Yifan 2013.05.05 */
		sd_global->bfcuts = new_bcuts(prob, BATCH_SIZE, sd_global->bfcuts);
		sd_global->bfcuts_pool = new_bcuts(prob, BATCH_SIZE,
				sd_global->bfcuts_pool);

		sd_global->quad_v = arr_alloc(BATCH_SIZE, double);
		/* modified by Yifan 2012.10.05 */
		sd_global->Obj_lb = arr_alloc(BATCH_SIZE,double);
	}

	if (!(cell->subprob = new_subprob(prob->subprob)))
		err_msg("Copy", "solve_cell", "cell->subprob");

#ifdef DEBUG
	printf("cell->master->lp = %d\n", cell->master->lp);
	printf("cell->subprob->lp = %d\n", cell->subprob->lp);
	printf("Writing problems\n");
	print_contents(cell->master, "master.out");
	print_contents(cell->subprob, "subprob.out");
	printf("Wrote our version\n");
	print_problem(cell->master, "master.mps");
	print_problem(cell->subprob, "subprob.mps");
	printf("Wrote CPLEX's version\n");
#endif

	/* In regularized QP method, we construct the Q matrix for the quadratic
	 term of the regularized master problem. Also we need to change the rhs 
	 and lower bounds from 'x' to 'd' here. zl */
	if (sd_global->config.MASTER_TYPE == SDQP)
	{
		construct_QP(prob, cell, cell->quad_scalar); //modified by Yifan 09/14/2011
		change_rhs(prob, cell, soln);
		change_bounds(prob, cell, soln);
	}

	fprintf(time_file, "Iter\t %s\t %s\t %s\t %s\t %s\t", "Iter_time", "Master",
			"Subprob", "Full_Test", "argmax");

	fprintf(time_file, "%s\t %s\t %s\t %s\t %s\n", "Iter_acc", "Master_acc",
			"Subprob_acc", "Full_Test_acc", "argmax_acc");

	/*Code below are added for evaluation!!!!!*/

#if 0
	fix=fopen("50samples_output.txt", "r");

	soln->mip_x[0]=0.0;
	for (i=1; i<=prob->num->mast_cols; i++)
	{
		fscanf(fix, "%lf",&soln->mip_x[i]);
		soln->incumb_x[0]+=soln->mip_x[i];
	}
	evaluate_inc(sd_global, cell,prob, soln, soln->mip_x, fname,conf_int,0);
	fclose(fix);
#endif

	/*Code above are added for evaluation!!!!!*/

	/*
	 while (!optimal(prob, cell, soln, phi) && cell->k < prob->num->iter)
	 */
	/* modified by zl, 06/30/04. 
	 while (!(soln->optimality_flag) && cell->k < prob->num->iter)
	 */
	/* modified by zl, 08/10/04. */
	write_statistics(sd_global, prob, cell, soln);

	while (TRUE)
	{
		//for (batch_id=0; batch_id<1; batch_id++) {

		iter_start_time = clock(); /* zl, 06/29/04. */
		/* The argmax procedures involve in two functions, namely,
		 form_new_cut and form_incumb_cut; solve_subprob() function is 
		 called in both solve_cell() in cell.c and in form_incumb_cut() 
		 in cuts.c. Therefore, we need to initiallize soln_subprob_iter 
		 and  _argmax_iter_ to zero at the beginning of every iteration, 
		 then accumulate the time spent on them throughout the iteration. 
		 zl, 06/30/04. */
		soln->run_time->soln_subprob_iter = 0.0;
		soln->run_time->argmax_iter = 0.0;

		/* Reset _smpl_test_flag_ to FALSE at the start of each iteration.
		 zl, 08/17/04. */
		soln->smpl_test_flag = FALSE;

		/* Yifan 03/23/2012 Reset the incumbent chage flag at the begining of each iteration*/
		/* So we assume the incumbent stays the same at the begining at the begining of each iteration*/
		incumb_change = FALSE;
		soln->incumbent_change = FALSE;

		/*print_vect(soln->Master_pi,prob->num->mast_rows+prob->num->max_cuts+1,"master_pi_cell.c");*//*debugging Yifan*/

		/* 1. Optimality Test */
		soln->optimality_flag = optimal(sd_global, prob, cell, soln, phi);
		/* Be careful here! Without this check the iteration will be off by
		 one from the original code. zl, 08/09/08.
		 And the number of optimality tests may be off. zl, 08/10/04. */
		//printf("This is opt-flag: %d, s->flag: %d \n",soln->optimality_flag, *soln->dual_statble_flag);
		/*Update code in this "if statement" for Batch-mean Yifan 2012-09-05*/
		if ((soln->optimality_flag && *soln->dual_statble_flag)
				|| cell->k >= prob->num->iter)
		{
			sd_global->ck[prob->current_batch_id] = cell->k;
			/* Record SD Lower Bound estimate modified by Yifan 2012.10.05 */
			sd_global->Obj_lb[prob->current_batch_id] = soln->incumb_est;

			record_master_dual_and_obj(sd_global, prob, soln, master_dual,
					master_obj);

			/* change_rhs_back(prob, cell, soln);*//*Keep the decisions as d instead of chaning them to x*/
			/* write_prob(cell->master, "final.lp"); */
			/* write_prob(sd_global->batch_problem, "final-1cut-a.lp");*/
            
            /* modified by Yifan 2013.9.18 */
			// update_batch_rhs(sd_global, prob, cell, soln, prob->current_batch_id);
			// update_batch_bounds(sd_global, prob, cell, soln, prob->current_batch_id);
            int ctype_idx;
            for (ctype_idx = 0; ctype_idx < prob->num->mast_cols; ctype_idx++) {
                sd_global->batch_problem->typex[ctype_idx + prob->current_batch_id * prob->num->mast_cols] = prob->master->typex[ctype_idx];
            }
            sd_global->batch_problem->typex[BATCH_SIZE*prob->num->mast_cols + prob->current_batch_id] = prob->master->typex[prob->num->mast_cols];

			save_batch_incumb(sd_global, prob, cell, soln,
					prob->current_batch_id);

			/*Group cuts structure under bcuts structure  Yifan 2013/01/17*/

			/* modified by Yifan 2013.02.15 */
			sd_global->bcuts->batch[prob->current_batch_id] = cell->cuts;
			sd_global->bfcuts->batch[prob->current_batch_id] =
					cell->feasible_cuts_added;
			sd_global->bfcuts_pool->batch[prob->current_batch_id] =
					cell->feasible_cuts_pool;

			for (i = 0; i < prob->num->max_cuts; i++)
			{
				add_cut_to_batch(sd_global, cell->cuts->val[i], prob, cell,
						soln, prob->current_batch_id, i);
			}

			if (prob->current_batch_id == BATCH_SIZE - 1)
			{

				/* Setup the whole problem as a quadratic problem */
                /* modified by Yifan 2013.9.18 */
				//construct_batch_QP(sd_global, prob, cell, 1.0 / BATCH_SIZE);
                copy_type(sd_global->batch_problem, sd_global->batch_problem->typex);
				/* modified by Yifan 2013.05.05 */
				if (0)
				{
					for (cnt = 0; cnt < BATCH_SIZE; cnt++)
					{
						for (i = 0; i < sd_global->bfcuts->batch[cnt]->cnt; i++)
						{
							for (j = 0; j < BATCH_SIZE; j++)
							{
								add_fcut_to_batch(sd_global,
										sd_global->bfcuts->batch[cnt]->val[i],
										prob, cell, soln, j, sd_global->ck[j]);
							}
						}
					}
				}
				/*
				 else
				 {
				 j=0;
				 for (cnt = 0; cnt < BATCH_SIZE; cnt++) {
				 for (i = 0; i < sd_global->bfcuts->batch[cnt]->cnt; i++) {
				 add_fcut_to_batch(sd_global->bfcuts->batch[cnt]->val[i], prob, cell, soln, cnt, ck[j]);
				 }
				 }
				 }*/

				add_batch_equality(sd_global, prob, cell, soln);

				/*Get batch dual and obj and corresponding statistics*/
				process_batch_prob(sd_global, prob, soln, &obj, batch_dual,
						batch_d, batch_obj);
			}

			/*record the quadratic sclar used in this replication*/
			sd_global->quad_v[prob->current_batch_id] = cell->quad_scalar;

			break;
		}

		/* 2. Generate Omega */
		soln->omega->k = cell->k;
      
#ifdef OMEGA_FILE
      omeg_idx = get_observ(sd_global, soln->omega, prob->num, &new_omega);
#else
      omeg_idx = generate_observ(sd_global, soln->omega, prob->num,
                                 &new_omega, &(sd_global->config.RUN_SEED));
#endif
      
#ifdef REC_OMEGA
      int i;
      /* Yifan 05/10/2012 Scenario Reader */
      FILE *omg;
      omg = fopen("sampleOmegas_run", "a");
      for(i=1; i <= sd_global->omegas.num_omega; i++){
        /* Note that RT start from 1 but mean start from 0*/
        fprintf(omg, "%f\t", soln->omega->RT[i]+sd_global->omegas.mean[i-1]);
      }
      fprintf(omg, "\n");
      fclose(omg);
#endif

#ifdef DEBUG
		for(idx=0;idx<soln->omega->most;idx++)
		if (valid_omega_idx(soln->omega, idx))
		print_omega(soln->omega, prob->num, idx);
#endif

		cell->k++;
		/* JH 
		 junk = fopen("Junk","a"); 
		 fprintf(junk,"....... k = %d \n", cell->k); 
		 fclose(junk); 
		 */

		/* 3. Solve Subproblem with candidate */
		if (!solve_subprob(sd_global, prob, cell, soln, soln->candid_x,
				omeg_idx))
		{
			cplex_err_msg(sd_global, "Subproblem", prob, cell, soln);
			break;
		}

		/* 4. Enter Feasibility Mode */
		if (prob->subprob->feaflag == FALSE)
		{
			resolve_infeasibility(sd_global, cell, prob, soln, omeg_idx,
					new_omega, phi);
		}

		/* 5. Form and update optimality cut */
		/* Recording the time spent on "argmax" procedures. zl, 06/30/04. */
		argmax_start = clock();

		form_new_cut(sd_global, prob, cell, soln, omeg_idx, new_omega);

		argmax_end = clock();
		soln->run_time->argmax_iter += ((double) (argmax_end - argmax_start))
				/ CLOCKS_PER_SEC;

		++num_subproblems;

		if (cell->cuts->cnt > 0 && soln->incumbent_change == TRUE)
		{
			soln->incumb_cut = cell->cuts->cnt - 1;
			soln->incumbent_change = FALSE;
		}

		/* 6. Update feasibility cut */
		if (cell->feasible_cuts_pool->cnt > 0 && prob->subprob->feaflag == TRUE)
		{
			/* Only update feasibility cuts when sub problem is feasible and we have a new omega */
			/* Since it is already updated in resolve_infeasibility if the subproblem is feasible*/
			if (new_omega == TRUE)
			{
				form_fea_cut(sd_global, prob, cell, soln, omeg_idx,
						soln->candid_x, new_omega);
			}
		}

		/* 7. Form incumbent cut */
		/* 8. Check improvement */

		if (cell->k > 1 && cell->cuts->cnt > 0)
		{
			if (cell->k - soln->last_update == prob->tau)
			{
				++num_subproblems; /* JH 3/12/98 */
				/* Because a solve_subprob function is called in form_incumb_cut, 
				 we have to go into the function to record the time spent on 
				 argmax procedures without counting solution time for solving 
				 subproblem repeatedly. zl, 06/30/04. */
				incumb_change = form_incumb_cut(sd_global, prob, cell, soln,
						omeg_idx);

				if (incumb_change == TRUE)
				{
					printf("incumbent will be replaced!\n");
					replace_incumbent(sd_global, prob, cell, soln, phi);
					/*
					 new_incumbent(prob, cell, soln, soln->candid_est);
					 construct_QP(prob, cell, cell->quad_scalar);
					 change_rhs(prob, cell, soln);
					 change_bounds(prob, cell, soln);*/
					/*???? quad scalar updates????*/

				}
			}
#ifdef TEST_SOLN
			if (1)
			{
				printf("c->k:%d\t Before Checking improvement!\n",cell->k);
				printf("Candidate_est is: %f\n",soln->candid_est);
				print_vect(soln->candid_x, prob->num->mast_cols, "Candidate X");
				printf("Incumbent_est is: %f\n",soln->incumb_est);
				print_vect(soln->incumb_x, prob->num->mast_cols, "Incumbent X");
			}
#endif
			/* check_improvement(prob, cell, soln); */
			if (incumb_change == FALSE)
			{
#ifdef TEST_SOLN
				printf("\nChecking Improvement...\n");
#endif
				check_improvement(sd_global, prob, cell, soln, phi); /*  SS messing around  */
			}
		}
		else
		{
			cell->quad_scalar = sd_global->config.MIN_QUAD_SCALAR;
			construct_QP(prob, cell, cell->quad_scalar);
		}

		/* Yifan 03/12/2012 Test for incumbent changes*/
#ifdef TEST_SOLN
		if (1)
		{
			printf("\nc->k:%d\n",cell->k);
			printf("Candidate_est is: %f\n",soln->candid_est);
			print_vect(soln->candid_x, prob->num->mast_cols, "Candidate X");
			printf("Incumbent_est is: %f\n",soln->incumb_est);
			print_vect(soln->incumb_x, prob->num->mast_cols, "Incumbent X");
			printf("__________________________________________\n");
		}
#endif

		/* 9. Solve QP master */
		/* Modified for regularized QP method. zl */
		if (sd_global->config.MASTER_TYPE == SDLP)
		{
			if (!solve_master(sd_global, prob, cell, soln))
			{
				cplex_err_msg(sd_global, "Master", prob, cell, soln);
				break;
			}
		}
		else
		{
			if (!solve_QP_master(sd_global, prob, cell, soln))
			{
				cplex_err_msg(sd_global, "QP_Master", prob, cell, soln);
				break;
			}
		}

		write_statistics(sd_global, prob, cell, soln);

		/*thin_data(prob, cell, soln);*/

		/* Yifan 03/12/2012 Test for incumbent changes*/
		if (0)
		{
			printf("Candidate_est is: %f\n", soln->candid_est);
			print_vect(soln->candid_x, prob->num->mast_cols, "Candidate X");
			printf("Incumbent_est is: %f\n", soln->incumb_est);
			print_vect(soln->incumb_x, prob->num->mast_cols, "Incumbent X");
		}

		/* Record UB, LB, and error, iteratively. zl, 06/30/04. */
		if (cell->k >= sd_global->config.MIN_ITER)
			fprintf(obj_file, "%d \t%f \t%f \t%f\n", cell->k, soln->incumb_est,
					soln->candid_est,
					(soln->incumb_est - soln->candid_est) / soln->incumb_est);

		/* At the end of an iteration, we accumulate the soln_suprob_accum 
		 and argmax_accum. zl, 06/30/04. */
		soln->run_time->argmax_accum += soln->run_time->argmax_iter;
		soln->run_time->soln_subprob_accum += soln->run_time->soln_subprob_iter;

		iter_end_time = clock(); /* zl, 06/29/04. */
		soln->run_time->iteration_time = ((double) (iter_end_time
				- iter_start_time)) / CLOCKS_PER_SEC;
		soln->run_time->iteration_accum += soln->run_time->iteration_time;

		fprintf(time_file, "%d\t %lf\t %lf\t %lf\t %lf\t %lf\t", cell->k,
				soln->run_time->iteration_time,
				soln->run_time->soln_master_iter,
				soln->run_time->soln_subprob_iter,
				soln->run_time->full_test_iter, soln->run_time->argmax_iter);
		fprintf(time_file, "%lf\t %lf\t %lf\t %lf\t %lf\n",
				soln->run_time->iteration_accum,
				soln->run_time->soln_master_accum,
				soln->run_time->soln_subprob_accum,
				soln->run_time->full_test_accum, soln->run_time->argmax_accum);

		//}
	}

	total_end_time = clock(); /* zl, 06/29/04. */
    stop = time(NULL);
    diff = difftime(stop, start);
    total_acc += diff;
    time_benders = fopen("time_benders", "a");
    fprintf(time_benders, "sd elapsed time:\t%f\t%f\n", diff, total_acc);
    fclose(time_benders);
    
#ifdef DEBUG
	times(&total_term_time); /* zl, 08/18/04. */
#endif
	soln->run_time->total_time = ((double) (total_end_time - total_start_time)) / CLOCKS_PER_SEC;
	/* zl, 08/18/04. */
    print_vect(soln->incumb_x, prob->num->mast_cols, "x_k from the incumbent problem");
    //evaluate_inc(sd_global, cell, prob, soln, soln->incumb_x, fname, conf_int, 0);

    
#if 1
    /* modified by Yifan 2013.10.02 */
    double obj_val,max_height, max_height2;
    int ssn_idx;
    max_height = max_cut_height(sd_global, cell->cuts, soln->incumb_x, cell, prob->num);
    /*soln->mip_x[1]= 4;
    soln->mip_x[2]= 4;
    soln->mip_x[3]= 7;
    soln->mip_x[4]= 4;
    max_height2 = max_cut_height(sd_global, cell->cuts, soln->mip_x, cell, prob->num);
    evaluate_inc(sd_global, cell, prob, soln, soln->mip_x, fname, conf_int, 0);*/

    destruct_QP(prob, cell);
    change_rhs_back(prob, cell, soln);
    change_bounds_back(prob, cell, soln);
    copy_arr(soln->mip_x, soln->candid_x, prob->num->mast_cols);
    add_box_bounds(sd_global, prob, cell, soln);
    copy_type(cell->master, prob->master->typex);
    print_problem(cell->master, "FirstMIP.lp");
    set_intparam(NULL, PARAM_SCRIND, CPX_OFF);
    CPXmipopt(env, cell->master->lp);
    CPXgetobjval (env, cell->master->lp, &obj_val);
    get_x(cell->master, soln->mip_x+1, 0, prob->num->mast_cols);
    soln->mip_x[0] = one_norm(soln->mip_x+1, prob->num->mast_cols);
    print_vect(soln->mip_x, prob->num->mast_cols, "MIP-X");
    
    for (ssn_idx=0; ssn_idx<prob->num->mast_cols; ssn_idx++) {
        soln->mip_x[ssn_idx]=0;
    }
    u_group_type u_group;
    for (ssn_idx=0; ssn_idx<89; ssn_idx++) {
        u_group.uvalue[ssn_idx].u = soln->incumb_x[1+ssn_idx];
        u_group.uvalue[ssn_idx].u_idx = 90+ssn_idx;
    }
    qsort(u_group.uvalue, sizeof(u_group.uvalue)/sizeof(*u_group.uvalue), sizeof(*u_group.uvalue), compare_function);
    double ssn_sum = 0, temp;
    for (ssn_idx=88; ssn_idx>=0; ssn_idx--) {
        if(ssn_sum > prob->master->rhsx[1]) break;
        soln->mip_x[u_group.uvalue[ssn_idx].u_idx] = 1;
        soln->mip_x[u_group.uvalue[ssn_idx].u_idx-89] = soln->incumb_x[u_group.uvalue[ssn_idx].u_idx-89];
        CPXgetcoef(env, cell->master->lp, 1, u_group.uvalue[ssn_idx].u_idx-1, &temp);
        ssn_sum += temp;
    }
    soln->mip_x[0] = one_norm(soln->mip_x+1, prob->num->mast_cols);
    print_vect(soln->mip_x, prob->num->mast_cols, "Adjusted-MIP-X");
    
    /* modified by Yifan 2013.10.04 */
#if 0
	fix=fopen("500samples_output.txt", "r");
    
	soln->mip_x[0]=0.0;
	for (i=1; i<=prob->num->mast_cols; i++)
	{
		fscanf(fix, "%lf",&soln->mip_x[i]);
		soln->mip_x[0]+=soln->mip_x[i];
	}
    sd_global->config.EVAL_SEED1 = prob->eval_seed;
	evaluate_inc(sd_global, cell,prob, soln, soln->mip_x, fname,conf_int,0);
	fclose(fix);
#endif
    
    //evaluate_inc(sd_global, cell, prob, soln, soln->mip_x, fname, conf_int, 0);
    int mip_idx = 1;
    char **mip_names;
    mip_names = arr_alloc(2, string);
    mip_names[0] = arr_alloc(NAME_SIZE, char);
    
    /* modified by Yifan 2013.9.26 */
    remove_all_cuts(sd_global, prob, cell, soln);
    print_problem(cell->master, "after_rv_all_cuts.lp");
    
    /* Allocate memory for the incumbent IP solution 2013.10.01 Yifan */
    double * incumbent_x_k;
    double pred_improve = 0.0, alpha = 0.1;
    double actual_improve = 0.0;
    int num_flip = 1;
    if (!(incumbent_x_k = arr_alloc(prob->num->mast_cols+1, double)))
        err_msg("Allocation", "solve_cell", "incumbent_x_k");
    copy_arr(incumbent_x_k, soln->mip_x, prob->num->mast_cols);
    sd_global->config.EVAL_SEED1 = prob->eval_seed;
    evaluate_inc(sd_global, cell, prob, soln, incumbent_x_k, fname, conf_int, 0);
    
    soln_benders = fopen("soln_benders", "a");
    fprintf(soln_benders, "Iteration: %d\n", 0);
    fprint_vect(soln_benders, incumbent_x_k, prob->num->mast_cols, "MIP_X:");
    fclose(soln_benders);
    
    start = time(NULL);
    while (mip_idx<=1000) {
        update_cont_bounds(sd_global, prob, cell, soln, incumbent_x_k, alpha, mip_idx);
        // omeg_idx = generate_observ(sd_global, soln->omega, prob->num, &new_omega, &(sd_global->config.RUN_SEED));
       /* if (mip_idx%50==0) {
            for (omeg_cnt = 0; omeg_cnt < soln->omega->most; omeg_cnt++) {
                if (!solve_subprob(sd_global, prob, cell, soln, incumbent_x_k,
                                   omeg_cnt))
                {
                    cplex_err_msg(sd_global, "Subproblem", prob, cell, soln);
                    break;
                }
                stochastic_updates(sd_global, cell, cell->lambda, cell->sigma, soln->delta, soln->omega, prob->num, prob->Rbar, prob->Tbar, cell->subprob, soln->Pi, omeg_cnt, FALSE);
            }
        }*/

        form_ip_cut(sd_global, prob, cell, soln, omeg_idx, FALSE);

        
        //add_box_bounds(sd_global, prob, cell, soln);
        add_flip_bounds(sd_global, prob, cell, soln,incumbent_x_k, &num_flip);
        if (mip_idx>1) {
            max_height = max_cut_height(sd_global, cell->cuts, soln->mip_x, cell, prob->num)+CxX(prob->c, soln->mip_x, prob->num->mast_cols);
            actual_improve = max_height2 - max_height;
            //if (actual_improve > 0.2 * pred_improve) {
            if (max_height < max_height2) {
                copy_arr(incumbent_x_k, soln->mip_x, prob->num->mast_cols);
                num_flip = min(num_flip*2, 89);
                alpha = min(alpha * 2.0, 0.99);
                //sd_global->config.EVAL_SEED1 = prob->eval_seed;
                //evaluate_inc(sd_global, cell, prob, soln, incumbent_x_k, fname, conf_int, 0);
            }
            else
            {
                num_flip = max(num_flip/2, 1);
                alpha = max(alpha*0.5, 0.1);
                //sd_global->config.EVAL_SEED1 = prob->eval_seed;
                //evaluate_inc(sd_global, cell, prob, soln, soln->mip_x, fname, conf_int, 0);
            }
        }
        strcpy(mip_names[0], "MIP    ");
        mip_names[0][3] = '0' + mip_idx / 10000 % 10;
        mip_names[0][4] = '0' + mip_idx / 1000 % 10;
        mip_names[0][5] = '0' + mip_idx / 100 % 10;
        mip_names[0][6] = '0' + mip_idx / 10 % 10;
        mip_names[0][7] = '0' + mip_idx / 1 % 10;
        mip_names[0][8] = '.';
        mip_names[0][9] = 'l';
        mip_names[0][10] = 'p';
        mip_names[0][11] = '\0';
        //print_problem(cell->master, mip_names[0]);
        max_height2 = max_cut_height(sd_global, cell->cuts, incumbent_x_k, cell, prob->num)+CxX(prob->c, incumbent_x_k, prob->num->mast_cols);
        set_intparam(NULL, PARAM_SCRIND, CPX_OFF);
        CPXmipopt(env, cell->master->lp);
        CPXgetobjval (env, cell->master->lp, &obj_val);
        pred_improve = max_height2 - obj_val;
        get_x(cell->master, soln->mip_x+1, 0, prob->num->mast_cols);
        soln->mip_x[0] = one_norm(soln->mip_x+1, prob->num->mast_cols);
        printf("Iteration%d\t:",mip_idx);
        //print_vect(soln->mip_x, prob->num->mast_cols, "MIP-X");
        printf("Upper Bound: %f\n", max_height2);
        printf("Lower Bound: %f\n", obj_val);
        printf("Gap: %f\%%\n",DBL_ABS((max_height2-obj_val)/max_height2)*100);
        if (mip_idx%100==0) {
            stop = time(NULL);
            diff = difftime(stop, start);
            total_acc += diff;
            time_benders = fopen("time_benders", "a");
            fprintf(time_benders, "%d Benders' solve:\t%f\t%f\n", mip_idx, diff,total_acc);
            fclose(time_benders);
            
            soln_benders = fopen("soln_benders", "a");
            fprintf(soln_benders, "Iteration: %d\n", mip_idx);
            fprint_vect(soln_benders, incumbent_x_k, prob->num->mast_cols, "MIP_X:");
            fclose(soln_benders);
            sd_global->config.EVAL_SEED1 = prob->eval_seed;
            evaluate_inc(sd_global, cell, prob, soln, incumbent_x_k, fname, conf_int, 0);
            start = time(NULL);
        }
        mip_idx++;
    }
    mem_free(mip_names[0]);
    mem_free(mip_names);
    //sd_global->config.EVAL_SEED1 = prob->eval_seed;
    //evaluate_inc(sd_global, cell, prob, soln, incumbent_x_k, fname, conf_int, 0);
    mem_free(incumbent_x_k);



#endif
#ifdef DEBUG
	printf("init: utime = %lld, stime = %lld, cutime = %lld, cstime = %lld\n",
			total_init_time.tms_utime, total_init_time.tms_stime,
			total_init_time.tms_cutime, total_init_time.tms_cstime);

	printf("term: utime = %lld, stime = %lld, cutime = %lld, cstime = %lld\n",
			total_term_time.tms_utime, total_term_time.tms_stime,
			total_term_time.tms_cutime, total_term_time.tms_cstime);

	printf("total: utime = %f, stime = %f, cutime = %f, sutime = %f\n",
			((double)(total_term_time.tms_utime-total_init_time.tms_utime))/CLK_TCK,
			((double)(total_term_time.tms_stime-total_init_time.tms_stime))/CLK_TCK,
			((double)(total_term_time.tms_cutime-total_init_time.tms_cutime))/CLK_TCK,
			((double)(total_term_time.tms_cstime-total_init_time.tms_cstime))/CLK_TCK);
#endif

	/*  JH summarizing output for Arlie and Brett ... 3/12/98 */

	/* removed by zl for online publication. */
#ifdef DEBUG
	printf(" Printing Summary \n");
#endif
	f2_out = fopen("Summary.out", "a");
	fprintf(f2_out, "\n\n r = %lf tau = %d epsilon = %lf  \n",
			sd_global->config.R, prob->tau, sd_global->config.EPSILON);
	fprintf(f2_out, "Iterations: %d   Subproblems %d \n", cell->k,
			num_subproblems);
	fprintf(f2_out, "SD Est. Inc. Value: %lf \n", soln->incumb_est);
	fprintf(f2_out, " total time = %lf \n", soln->run_time->total_time); /*modified by Yifan 10/04/2011*/
  
    /* modified by Yifan 2013.06.30 */
    time_sample = fopen("time_sample.out", "a"); 
    fprintf(time_sample, "%lf\t%d\n",soln->run_time->total_time,cell->k);
    fclose(time_sample);

    
	for (i = 1; i <= prob->num->mast_cols; ++i)
		x_k[i] = soln->incumb_x[i];

	/* Try an evaluation at the end of the 1st replication */
	/* Should be removed during Batch-Mean Run!!!*/

#if 0
	fix=fopen("incumb-pgp2.txt", "r");

	for (j=1; j<=1; j++)
	{
		soln->incumb_x[0]=0.0;
		for (i=1; i<=prob->num->mast_cols; i++)
		{
			fscanf(fix, "%lf",&soln->incumb_x[i]);
			soln->incumb_x[0]+=soln->incumb_x[i];
		}
		if (sd_global->config.EVAL_FLAG==1)
		{	evaluate_inc(cell,prob, soln, soln->incumb_x, fname,conf_int);
		}
	}

	fclose(fix);
#endif

	if (sd_global->config.EVAL_FLAG == 1)
	{
		sd_global->config.EVAL_SEED1 = prob->eval_seed;
		evaluate_inc(sd_global, cell, prob, soln, soln->incumb_x, fname,
				conf_int, 0);
		fprintf(f2_out, "Inc. Value 0.95 CI: [%lf , %lf] \n", conf_int[0],
				conf_int[1]);
	}
  
  if (prob->current_batch_id == 2) {
    double spn[3],spn_max;
    vector *rd;
    calc_var(sd_global, sd_global->Obj_lb, &(sd_global->obj_mean), &(sd_global->obj_stdev), 3);
    if (DBL_ABS(sd_global->obj_stdev/sd_global->obj_mean) < sd_global->config.TOLERANCE) {
      sd_global->obj_flag = 1;
    }
    if(!(rd = (vector *) mem_calloc (3, sizeof(vector))))
      err_msg("Allocation", "solve_cell", "rd");
    for (i=0; i<3; i++) {
      rd[i] = arr_alloc(prob->num->mast_cols + 1, double);
    }
    
    /* Calculate the mean solution from first 3 replications */
    for (i=1; i<=prob->num->mast_cols; i++) {
      soln->incumb_avg[i]=0.0;
      for (j=0; j<3; j++) {
        soln->incumb_avg[i]+=sd_global->batch_incumb->incumb_x[j][i];
      }
      soln->incumb_avg[i]/=3.0;
    }
    soln->incumb_avg[0] = one_norm(soln->incumb_avg+1, prob->num->mast_cols);
    
    for (i = 0 ; i < 3 ; i++) {
      rdiff(sd_global, sd_global->batch_incumb->incumb_x[i], soln->incumb_avg, prob->num->mast_cols, rd[i]);
      spn[i] = sup_norm(rd[i], prob->num->mast_cols);
    }
    
    spn_max = spn[0];
    for (i = 1; i < 3; i++) {
      if (spn[i] > spn_max)
        spn_max = spn[i];
    }
    
    for (i = 0; i < 3; i++) {
      mem_free(rd[i]);
    }
    mem_free(rd);
    
    if (spn_max < sd_global->config.TOLERANCE) {
      sd_global->average_flag = 1;
    }
    
  }
  
  /* Calculate mean and stdev of early terminated mean solution and print the output */
  if (prob->current_batch_id == 2 && (sd_global->average_flag == 1 || sd_global->obj_flag == 1 ) && sd_global->config.OVERRIDE == 1) {
    /* Calculate mean and stdard deviation of duals and reduced costs for mean solution*/
    calc_mean_stdev(sd_global->batch_incumb->R_Master_pi, soln->R_Master_pi_mean, soln->R_Master_pi_stdev, prob->num->mast_rows, 3);
    calc_mean_stdev(sd_global->batch_incumb->R_Master_dj, soln->R_Master_dj_mean, soln->R_Master_dj_stdev, prob->num->mast_cols, 3);
    calc_var(sd_global, sd_global->Obj_lb, &soln->Obj_lb_mean, &soln->Obj_lb_stdev, 3);
    
    /* Evaluate the mean solution */
    sd_global->config.EVAL_SEED1 = prob->eval_seed;
    evaluate_inc(sd_global, cell, prob, soln, soln->incumb_avg, fname, conf_int, 3);
    print_detailed_soln(sd_global, soln, prob, fname, 3);
  }
  

	/* modified by Yifan 2012.09.23 */
	/* Evaluate batch x and average x of replicatiing x 2012.09.23 */

	if (prob->current_batch_id == BATCH_SIZE - 1)
	{

		/* modified by Yifan 2013.02.15 */
		/* calc_lowerbound(prob, soln); */

		batch_d = fopen("Batch_x.out", "a");
		soln->incumb_d[0] = one_norm(&soln->incumb_d[1], prob->num->mast_cols);
		sd_global->config.EVAL_SEED1 = prob->eval_seed;
		evaluate_inc(sd_global, cell, prob, soln, soln->incumb_d, fname,
				conf_int, 1);
		fprintf(batch_d, "Batch_x - Inc. Value 0.95 CI: [%lf , %lf] \n",
				conf_int[0], conf_int[1]);
		sd_global->config.EVAL_SEED1 = prob->eval_seed;
		evaluate_inc(sd_global, cell, prob, soln, soln->incumb_avg, fname,
				conf_int, 2);
		fprintf(batch_d, "Incumb_avg - Inc. Value 0.95 CI: [%lf , %lf] \n",
				conf_int[0], conf_int[1]);
		fclose(batch_d);
		print_detailed_soln(sd_global, soln, prob, fname, 1);
		print_detailed_soln(sd_global, soln, prob, fname, 2);
	}

	if (sd_global->config.DETAILED_SOLN == 1)
	{
		print_detailed_soln(sd_global, soln, prob, fname, 0);
	}

#ifdef DEBUG
	printf(" Ending Summary \n");
#endif


	fclose(f2_out);

	/* modified by Yifan 2013.05.05 */
	/* print_soln(prob, cell, soln, fname);*/

	fclose(obj_file);
	fclose(time_file);

	free_subprob(cell->subprob);
	free_master(cell->master);
	free_soln(prob, cell, soln);
}

/************************************************************************\
** This function creates a brand spankin new cell structure.
 ** By that I mean that the structures within the cell, like
 ** lambda, sigma, and cuts, do not have any data already initialized.
 ** For a single cell problem, or for first generation ancestors,
 ** this is standard; but for the parallel implementation, the new cell 
 ** must be read in from a file, and its individual structures will 
 ** not be empty.
 \************************************************************************/
cell_type *new_cell(sdglobal_type* sd_global, prob_type *p, int id_num)
{
	cell_type *c;
	int length;

#ifdef TRACE
	printf("Inside new_cell\n");
#endif

	if (!(c = (cell_type *) mem_malloc (sizeof(cell_type))))
		err_msg("Allocation", "new_cell", "cell");
	/* Yifan 06/18/2012 batch mean */

	c->members = arr_alloc(1, int); /* Make room for one member */
	c->members[0] = id_num; /* That member is itself */
	c->num_members = 1; /* Establish the number of members as 1 */
	c->quad_scalar = sd_global->config.MIN_QUAD_SCALAR; /* The quadratic scalar, 'sigma'. zl */
	c->id_num = id_num; /* Give id # according to parameter */
	c->P = 1.0; /* This cell has all the probability */
	c->N = 0; /* Cell starts with no ancestor cuts */
	c->k = 0; /* Cell starts with no cuts of its own */
	/* Equivalently, starts at 0th iteration */
	c->LP_cnt = 0; /* Cell starts with no LP solved. zl */
	c->LP_test = 0; /* Cell starts with no LP solved for
	 testing. zl */
	c->opt_mode = TRUE;
	c->incumb_infea = FALSE;
	c->fea_count = 0;

	/* modified by Yifan 2013.02.15 */
	length = p->num->iter + p->num->iter / p->tau + 1;

	/* Yifan 06/18/2012 batch mean */
	c->cuts = new_cuts(p->num->iter, p->num->mast_cols, 0);
	c->lambda = new_lambda(length, 0, p->num->rv_rows, p->coord);
	c->sigma = new_sigma(length, p->num->nz_cols, 0, p->coord);
	c->theta = new_theta(0);

	/* Yifan 03/04/2012 Updated for Feasibility Cuts*/
	c->feasible_cuts_pool = new_cuts(p->num->iter, p->num->mast_cols, 0);
	c->feasible_cuts_added = new_cuts(p->num->iter, p->num->mast_cols, 0);
	c->feasible_lambda = new_lambda(length, 0, p->num->rv_rows, p->coord);
	c->feasible_sigma = new_sigma(length, p->num->nz_cols, 0, p->coord);
	c->feasible_theta = new_theta(0);
	/* Yifan 03/04/2012 Updated for Feasibility Cuts*/

	return c;
}

/************************************************************************\
** This function releases all the memory reserved for the cell
 ** data structure, including the cel itself.
 \************************************************************************/
void free_cell(cell_type *c, num_type *num)
{

#ifdef TRACE
	printf("Inside free_cell\n");
#endif

#if 0
	printf("c->lambda->cnt is : %d\n",c->lambda->cnt);
	printf("c->sigma->cnt is : %d\n",c->sigma->cnt);
	printf("c->theta->cnt is : %d\n",c->theta->cnt);
	printf("c->feasible_lambda->cnt is : %d\n",c->feasible_lambda->cnt);
	printf("c->feasible_sigma->cnt is : %d\n",c->feasible_sigma->cnt);
	printf("c->feasible_theta->cnt is : %d\n",c->feasible_theta->cnt);
#endif
	/* Avoid freeing the cuts structure for compromise problem Yifan 2013/01/17 */
	//free_cuts(c->cuts);
	free_lambda(c->lambda);
	free_sigma(c->sigma);
	free_theta(c->theta);

	/* modified by Yifan 2013.05.05 */
	/* Yifan 03/04/2012 Updated for Feasibility Cuts*/
	/* free_cuts(c->feasible_cuts_pool); */
	/* modified by Yifan 2013.04.15 */
	//free_cuts(c->feasible_cuts_added);
	/* Yifan 03/11/12 Since feasible cuts added are cuts_type pointers*/
	/* mem_free(c->feasible_cuts_added->val);
	 mem_free(c->feasible_cuts_added);*/

	free_lambda(c->feasible_lambda);
	free_sigma(c->feasible_sigma);
	free_theta(c->feasible_theta);
	/* Yifan 03/04/2012 Updated for Feasibility Cuts*/

	mem_free(c->members);
	mem_free(c);
}

/************************************************************************\
** This function prints out the contents of a cell to a file.  This 
 ** essentially includes lambda, sigma, and the cuts, as well as the
 ** counters N, k, the cell's number id_num, and the probability P.
 \************************************************************************/
void write_cell(prob_type *p, cell_type *c, char *fname)
{
	FILE *fptr;
	char filename[NAME_SIZE * 2];

	strcpy(filename, fname);
	strcat(filename, ".dat");
	fptr = fopen(filename, "w");

	write_sigma(fptr, c->sigma, p->num);
	write_lambda(fptr, c->lambda, p->num);
	/*
	 write_cuts(fptr, c->cuts);
	 */

	fclose(fptr);
}
