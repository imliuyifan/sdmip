/***********************************************************************\
** sd.c
 ** 
 ** This program solves stochastic linear programs using stochastic
 ** decomposition.  The possible usages are:
 **
 **   sd	-- Prompts the user for one problem and solves it.
 **
 **   sd n -- Solves _n_ problems named "prob0000", "prob0001", etc.
 **
 **   sd n m -- Solves _n_ problems (as above) starting with "prob000m".
 **
 **   sd fname obj s1 s2 -- Solves the problem specified by _fname_ having
 **                         an objective sense of _obj_ (1 for minimization),
 **                         using random number seeds _s1_ and _s2_.
 **
 **
 ** This file contains the following functions:
 **   main()
 **   err_msg()
 **   parse_cmd_line()
 **
 \***********************************************************************/

#include "prob.h"
#include "cell.h"
#include "soln.h"
#include "quad.h"
#include <time.h>
#include "solver.h"
#include "utility.h"
#include "input.h"
#include "log.h"
#include "cuts.h"
#include "sdglobal.h"
#ifdef SD_win
#include <windows.h>
#endif
void sd_create_output_folder(char *buffer1, char *buffer2, char *fname);
#ifdef SD_win
void sd_mv_files(char *fname);
#else
void sd_mv_files(char *buffer1, char *buffer2, char *fname);
#endif

/* Begin SD */
int main(int argc, char *argv[])
{
	int num_rv; /* Total number of rv coefficients */
	int num_cipher; /* Number of ints needed to encode an omega */
	sd_small row, col; /* Location of master/subproblem breakpoint */
	double *x_k; /* Initial candidate, from original problem */
	double *original_x_k; /* Keep a copy of mean value solution so that multiple relication can run without solving mean value problem again */
	one_problem *probptr; /* Original LP, stored in CPLEX format */
	identity *ident; /* Identities of row and column names */
# if 0
	one_problem *orig; /* trial of original. zl */
	identity *id; /* trial of ident. zl */
#endif
	char fname[NAME_SIZE];/* Name of file containing the problem */
	int num_probs; /* Number of problems to be solved */
	int n, start; /* # of the current problem & 1st problem */
	int objsen; /* 1 for minimization, -1 for maximization */
# if 0
	FILE *sol_out; /*For writing solution to a file. zl */
# endif
	BOOL read_seeds; /* True if random # seeds should be read */
	BOOL read_iters; /* True if MIN_ITER and MAX_ITER should be read. zl 06/18/02 */
	double objective;
	int cnt = 0;
    sd_long seed1[BATCH_SIZE], seed2[BATCH_SIZE];
//	double eps[3];
//	int scan_len[3];
	int idx = 0;
	char buffer1[128], buffer2[128];
    
    sdglobal_type *sd_global;
    sd_global = (sdglobal_type *) mem_malloc(sizeof(sdglobal_type)); /* All global variables */
	sd_global->MEM_USED = 0;
	sd_global->Abar = 0;
	sd_global->MALLOC = 0;
    sd_global->average_flag = 0;
    sd_global->obj_flag = 0;

	/* Initialize the CPLEX environment. zl */
	open_Solver();

#ifdef CAL_CHECK 
	/* Open the file to trace the modification to regularized QP. zl */
	g_FilePointer = fopen("quad.out", "w+");
#endif

	/* Prepare to start the algorithm */
	log_start(sd_global);
	/* modified by zl. */
# if 0
	sol_out = fopen("sd.sol","w");
# endif

	printf("\nBeginning SD...\n\n");
	parse_cmd_line(sd_global, argc, argv, fname, &objsen, &num_probs, &start,
			&read_seeds, &read_iters);
	/* Load the solution settings */
	if (!load_config(sd_global, read_seeds, read_iters))
		return 1;
#ifdef SD_unix
    sd_create_output_folder(buffer1,buffer2,fname);
#endif

  
  //		eps[0] = 0.01;
  //		eps[1] = 0.001;
  //		eps[2] = 0.0001;
  
  //		scan_len[0] = 64;
  //		scan_len[1] = 256;
  //		scan_len[2] = 512;

	if (sd_global->config.MULTIPLE_REP == 1)
	{

		if (sd_global->config.AUTO_SEED == 0)
		{
			seed1[0] = sd_global->config.RUN_SEED1;
			seed1[1] = sd_global->config.RUN_SEED2;
			seed1[2] = sd_global->config.RUN_SEED3;
			seed1[3] = sd_global->config.RUN_SEED4;
			seed1[4] = sd_global->config.RUN_SEED5;
			seed1[5] = sd_global->config.RUN_SEED6;
			seed1[6] = sd_global->config.RUN_SEED7;
			seed1[7] = sd_global->config.RUN_SEED8;
			seed1[8] = sd_global->config.RUN_SEED9;
			seed1[9] = sd_global->config.RUN_SEED10;
			seed1[10] = sd_global->config.RUN_SEED11;
			seed1[11] = sd_global->config.RUN_SEED12;
			seed1[12] = sd_global->config.RUN_SEED13;
			seed1[13] = sd_global->config.RUN_SEED14;
			seed1[14] = sd_global->config.RUN_SEED15;
			seed1[15] = sd_global->config.RUN_SEED16;
			seed1[16] = sd_global->config.RUN_SEED17;
			seed1[17] = sd_global->config.RUN_SEED18;
			seed1[18] = sd_global->config.RUN_SEED19;
			seed1[19] = sd_global->config.RUN_SEED20;
			seed1[20] = sd_global->config.RUN_SEED21;
			seed1[21] = sd_global->config.RUN_SEED22;
			seed1[22] = sd_global->config.RUN_SEED23;
			seed1[23] = sd_global->config.RUN_SEED24;
			seed1[24] = sd_global->config.RUN_SEED25;
			seed1[25] = sd_global->config.RUN_SEED26;
			seed1[26] = sd_global->config.RUN_SEED27;
			seed1[27] = sd_global->config.RUN_SEED28;
			seed1[28] = sd_global->config.RUN_SEED29;
			seed1[29] = sd_global->config.RUN_SEED30;
		}
		else
		{
			generate_seed(seed1, seed2);
		}

#if 0
		/* A */
		seed1[0] = 9495518635394380;
		seed1[1] = 4650175399072632;
		seed1[2] = 6070772756632709;
		seed1[3] = 5451675876709589;
		seed1[4] = 5285327724846206;
		seed1[5] = 5588857889468088;
		seed1[6] = 1098833779416153;
		seed1[7] = 6192593982049265;
		seed1[8] = 4756774140130874;
		seed1[9] = 6784592265109609;
		seed1[10] = 9728429908537680;
		seed1[11] = 1163479388309571;
		seed1[12] = 3279282318700126;
		seed1[13] = 8773753208032360;
		seed1[14] = 9337302665697748;
		seed1[15] = 4415169667296773;
		seed1[16] = 4220432037464045;
		seed1[17] = 3554548844580680;
		seed1[18] = 1814300451929103;
		seed1[19] = 5339672949292608;
		seed1[20] = 5638710736762732;
		seed1[21] = 3154245808720589;
		seed1[22] = 2414929536171258;
		seed1[23] = 7998609999427572;
		seed1[24] = 7080145164625719;
		seed1[25] = 3648862740490586;
		seed1[26] = 7772725003305823;
		seed1[27] = 5982768791029230;
		seed1[28] = 1395182510837913;
		seed1[29] = 3735836402047426;
#endif

#if 0
		/* B */
		seed1[0] = 8879657642464524;
		seed1[1] = 4317459304174907;
		seed1[2] = 1499740298834250;
		seed1[3] = 8272809468603661;
		seed1[4] = 9321928632105101;
		seed1[5] = 8879657642464524;
		seed1[6] = 1646307759053034;
		seed1[7] = 1397125657640682;
		seed1[8] = 3146928660304649;
		seed1[9] = 6086062973158789;
		seed1[10] = 4261811376433110;
		seed1[11] = 5160431490422796;
		seed1[12] = 7210299483505433;
		seed1[13] = 2742341912700425;
		seed1[14] = 1085010081252686;
		seed1[15] = 8513449869606798;
		seed1[16] = 7093281297971938;
		seed1[17] = 7988825411001281;
		seed1[18] = 4183664541491746;
		seed1[19] = 3145719174690472;
		seed1[20] = 7565122826024890;
		seed1[21] = 5245385869406164;
		seed1[22] = 2209547377191484;
		seed1[23] = 9707622650545090;
		seed1[24] = 3276474213926122;
		seed1[25] = 3808908035978675;
		seed1[26] = 7200786232249811;
		seed1[27] = 3531095045851544;
		seed1[28] = 8536356961121783;
		seed1[29] = 4742397086462006;
#endif

#if 0
		/* C */
		seed1[0] = 8654069019158467;
		seed1[1] = 4785819713548623;
		seed1[2] = 7975937892764316;
		seed1[3] = 8255008427953745;
		seed1[4] = 8355970923248727;
		seed1[5] = 3835388830495228;
		seed1[6] = 3075928230436703;
		seed1[7] = 4284039266515684;
		seed1[8] = 7473499074508686;
		seed1[9] = 3670333759489016;
		seed1[10] = 2594208670938457;
		seed1[11] = 5727262445660711;
		seed1[12] = 4759600567514105;
		seed1[13] = 8478969444457220;
		seed1[14] = 5343940687542743;
		seed1[15] = 2814273299089456;
		seed1[16] = 7989132520505531;
		seed1[17] = 9986441112094174;
		seed1[18] = 8894979824116654;
		seed1[19] = 8691287178579158;
		seed1[20] = 8275759857380763;
		seed1[21] = 9976220018928872;
		seed1[22] = 4309613200835884;
		seed1[23] = 7874709394061938;
		seed1[24] = 7893592451466247;
		seed1[25] = 7108083357568830;
		seed1[26] = 8598564187996089;
		seed1[27] = 1348854001145810;
		seed1[28] = 3690791960107163;
		seed1[29] = 1627698145573959;
#endif

#if 0
		/* D */
		seed1[0] = 3893690411122860;
		seed1[1] = 5316036536755625;
		seed1[2] = 4572587611768476;
		seed1[3] = 9749041628135898;
		seed1[4] = 9272577758000225;
		seed1[5] = 6488432224836808;
		seed1[6] = 8196962071784000;
		seed1[7] = 5068262667087135;
		seed1[8] = 7921215418744271;
		seed1[9] = 7354688790019294;
		seed1[10] = 2869256240883701;
		seed1[11] = 4035953281440061;
		seed1[12] = 8880023876362669;
		seed1[13] = 7963210821217343;
		seed1[14] = 9932123249461067;
		seed1[15] = 9527048146300785;
		seed1[16] = 9667120111877512;
		seed1[17] = 7789520796672265;
		seed1[18] = 7562396690748148;
		seed1[19] = 6314775191587918;
		seed1[20] = 6511551822768524;
		seed1[21] = 3697588857961819;
		seed1[22] = 9518273207359016;
		seed1[23] = 8406157326186076;
		seed1[24] = 8186069662915543;
		seed1[25] = 1892784470925108;
		seed1[26] = 9047130629420280;
		seed1[27] = 2324213564395905;
		seed1[28] = 5555551107274368;
		seed1[29] = 9048361604101956;
#endif

#if 0
		/* E */
		seed1[0] = 6356668455071302;
		seed1[1] = 2855485413436641;
		seed1[2] = 3291598047866632;
		seed1[3] = 4741521067472913;
		seed1[4] = 7017703696683477;
		seed1[5] = 1292702333137961;
		seed1[6] = 5029588156470839;
		seed1[7] = 6705513938623459;
		seed1[8] = 8161200147085033;
		seed1[9] = 4194219324364584;
		seed1[10] = 6772355259593134;
		seed1[11] = 7109000029919588;
		seed1[12] = 2127430095908921;
		seed1[13] = 2110825389619457;
		seed1[14] = 6784662275630456;
		seed1[15] = 2329009110634238;
		seed1[16] = 7388820752166162;
		seed1[17] = 2186884426129148;
		seed1[18] = 8206146410102893;
		seed1[19] = 1659601989241207;
		seed1[20] = 6925496597308666;
		seed1[21] = 1769875597674400;
		seed1[22] = 1761840495746583;
		seed1[23] = 6935240926453844;
		seed1[24] = 8841187624726444;
		seed1[25] = 4805222323630005;
		seed1[26] = 8588410407071933;
		seed1[27] = 1191114992601797;
		seed1[28] = 6229834395460784;
		seed1[29] = 9974531442625448;
#endif

	}
	else
	{
		if (sd_global->config.AUTO_SEED == 0)
		{
			seed1[0] = sd_global->config.RUN_SEED1;
		}
		else
		{
			generate_seed(seed1, seed2);
		}
	}

	seed2[0] = sd_global->config.EVAL_SEED1;

	/*
	 printf("zl_sd In sd.c, main(), after load_config(), MIN = %d, MAX = %d\n", sd_global->config.MIN_ITER, sd_global->config.MAX_ITER);
	 */

	/* Solve all of the problems */
	for (n = 0; n < num_probs; n++)
	{
		/* Find the name of the current problem */
		if (argc >= 2 && argc < 5)
			filename_number(fname, 4, 1000, n + start);
		/* Copy the number _n+start_ into _fname_ starting at position 4 */

		/* If one of these fails, we'll have a lot of unfreed memory! */
#if 0
		if(load_core_cpx(&orig, &id, fname, objsen))
		{ /* zl trial. */
			if(load_core(&probptr, &ident, fname, objsen))
			{
#endif
		if (load_core_cpx(&probptr, &ident, fname, objsen))
		{ /* zl trial. */
			if (load_stoch(sd_global, probptr, ident, fname))
			{
				if (load_time(&row, &col, ident, fname))
				{
					/*** JH tying min_iter to first stage dimension 3/12/98 ***/
					if (sd_global->config.MIN_ITER
							< sd_global->config.ITER_FACT * col)
						sd_global->config.MIN_ITER = sd_global->config.ITER_FACT
								* col;

					/*
					 printf("zl_sd In sd.c, main(), after load_time, MIN = %d, MAX = %d\n", sd_global->config.MIN_ITER, sd_global->config.MAX_ITER);
					 */
#ifdef OMEGA_FILE
                  strcpy(sd_global->omegas.file_name, probptr->name);
                  strcat(sd_global->omegas.file_name, "Omegas");
                  sd_global->fptrOMEGA = fopen(sd_global->omegas.file_name, "r");
#endif
                  
					/* Solve the original combined problem */
					change_solver_primal(probptr); //added by Yifan to change solver to primal simplex
					setup_problem(probptr);
                    write_prob(probptr, "orig.lp");
					solve_problem(sd_global, probptr);

#if 0
					/* Write out the problem for checking purpose while doing
					 regularized QP method. Will be removed later. zl */
					write_prob(probptr, "orig.lp");
#endif

					/* Get the solution from the mean problem */
					if (!(x_k = arr_alloc(probptr->mac+1, double)))
						err_msg("Allocation", "main", "x_k");
					/* Allocate memory for the copy of mean value solution 04/25/2013 Yifan */
					if (!(original_x_k = arr_alloc(probptr->mac+1, double)))
						err_msg("Allocation", "main", "x_k");

					get_primal(x_k, probptr, probptr->mac);
					/* Make a copy of the mean value solution Yifan 04/25/2013 */
					copy_arr(original_x_k, x_k, probptr->mac);

					objective = get_objective(probptr);

					/* First stage solutions to the mean value problem Yifan*/
					print_vect(x_k, probptr->mac, "x_k from the mean problem");

					printf("Objective from the mean problem : %f\n", objective);
					get_lower_bound(sd_global, probptr, row, col);
					remove_problem(probptr);

#if 0    
					printf("~1\n");

					/* Save the mean solution to a data file */
					fprintf(sol_out, "Initial Solution: %d\n", n);
					for (j = 1; j <= col; ++j)
					if (x_k[j] != 0) fprintf(sol_out, " %s = %lf \n",
							probptr->cname[j-1], x_k[j]);

					printf("~2\n");
#endif 

					/* Now split the problem and solve it with SD */
					num_rv = sd_global->omegas.num_omega;
					num_cipher = sd_global->omegas.num_cipher;
					printf("\nnum_rv = %d;\n", num_rv);
					printf("num_cipher = %d;\n\n", num_cipher);

					if (sd_global->config.MULTIPLE_REP == 1)
					{
						for (idx = 1; idx < 2; idx++)
						{/* this layer control the tolerance of SD runs */
							for (cnt = 0; cnt < BATCH_SIZE; cnt++)
							{ /* this layer control the seed used in each SD replication */
								/*sd_global->config.EPSILON = eps[idx];
								 sd_global->config.SCAN_LEN = scan_len[idx];*/
								sd_global->config.EVAL_SEED1 = seed2[0];
								sd_global->config.RUN_SEED = seed1[cnt];

								/* Take the mean value solution as the initial candidate solution 04/25/2013 Yifan */
								copy_arr(x_k, original_x_k, probptr->mac);
								solve_SD(sd_global, probptr, x_k, num_rv,
										num_cipher, row, col, fname, cnt);
                              if (cnt == 2 && (sd_global->average_flag == 1 || sd_global->obj_flag == 1 ) && sd_global->config.OVERRIDE == 1) {
                                break;
                              }

							}
						}
					}
					else
					{
                        sd_global->config.RUN_SEED = seed1[0];
						solve_SD(sd_global, probptr, x_k, num_rv, num_cipher,
								row, col, fname, cnt);
					}

#if 0
					printf("~3\n");

					/* Save SD's solution to the same file, for comparison */
					fprintf(sol_out, "Recommended Solution: %d \n", n);
					for (j = 1; j <= col; ++j)
					if (x_k[j] != 0) fprintf(sol_out, " %s = %lf \n",
							ident->col_name[j-1], x_k[j]);

					printf("~4\n");
#endif
					mem_free(x_k);
					/* Release the copy of mean value solution Yifan 04/25/2013 */
					mem_free(original_x_k);
				} /* End of load_time. */
				else
				{
					printf("File error in load_time, problem #%d", n);
					return 1;
				}
			} /* End of load_stoch. */
			else
			{
				printf("File error in load_stoch, problem #%d", n);
				return 1;
			}
		} /* End of load_core_cpx.  zl trial. */
		else
		{
			printf("Failed in load_core_cpx, problem #%d", n);
			return 1;
		}

#if 0
	} /* End of load_core. */
	else
	{	printf("File error in load_core, problem #%d", n); return 1;}
} /* End of load_core_cpx.  zl trial. */
else
{	printf("Failed in load_core_cpx, problem #%d", n); return 1;}
#endif

		/* Clean up */
		free_ident(ident);
        free_omegas(sd_global);
      
#ifdef OMEGA_FILE
      fclose(sd_global->fptrOMEGA);
#endif


#if 0
		printf("~5\n");
		free_ident(id); /* zl_free */
		printf("~6\n");
		free_one_prob(orig); /* zl_free */
		printf("~7\n");
#endif
	}

	printf("\nEnding SD...\n\n");
	/* modified by zl*/
# if 0
	fclose(sol_out);
# endif

	log_stop(sd_global);

#ifdef CAL_CHECK 
	/* Close the file tracing the modification to regularized QP. zl */
	fclose(g_FilePointer);
#endif

	/*Yifan 2012-09-14*/
	mem_free(sd_global->batch_incumb->incumb_x);
	mem_free(sd_global->batch_incumb->R_Master_pi);
	/* modified by Yifan 2012.10.05 */
	mem_free(sd_global->batch_incumb->R_Master_dj);
	/* modified by Yifan 2012.10.05 */
	mem_free(sd_global->batch_incumb);
	mem_free(sd_global->Obj_lb);
	mem_free(sd_global->quad_v);

	/* modified by Yifan 2013.02.15 */
	free_bcuts(sd_global->bcuts);
	/* modified by Yifan 2013.05.05 */
	free_bcuts(sd_global->bfcuts_pool);
	/* modified by Yifan 2013.05.05 Free sd_global->bfcuts in the following manner
	 since sd_global->bfcuts only store pointers that point to sd_global->bfcuts_pool */
  for (idx=0; idx<BATCH_SIZE; idx++) {
    if (sd_global->bfcuts->batch[idx]!= NULL) {
      if (sd_global->bfcuts->batch[idx]->val!=NULL) {
        mem_free(sd_global->bfcuts->batch[idx]->val);
      }
    }
    if (sd_global->bfcuts->batch[idx]!=NULL) {
      mem_free(sd_global->bfcuts->batch[idx]);
    }
  }
  
    mem_free(sd_global->bfcuts->batch);
	mem_free(sd_global->bfcuts);

	free_one_prob(sd_global->batch_problem);
	free_one_prob(probptr);

	/* Release the CPLEX environment. zl */
	close_Solver();
#ifdef SD_win
    sd_mv_files(fname);
#else
    sd_mv_files(buffer1, buffer2, fname);
#endif
	return 1;
}

/****************************************************************************\
**   This is a stub, for use in place of the real optimality check.
 **
 BOOL optimal()
 {
 char	ch;

 #ifdef RUN
 return FALSE;
 #else
 
 return FALSE;

 printf("\n\n    Is this solution optimal? (y/n)");
 getchar();
 scanf("%c", &ch);
 if (ch == 'y' || ch == 'Y')
 return TRUE;
 else return FALSE;
 #endif

 }
 \****************************************************************************/
void sd_create_output_folder(char *buffer1, char *buffer2, char *fname)
{
    int status;
    strcpy(buffer1, "mkdir ./sdoutput");
	status = system(buffer1);
    if(status == -1){
        printf("system() call fails.\n");
        exit(1);
    }
	strcat(buffer1, "/");
	strcat(buffer1, fname);
	status = system(buffer1);
    if(status == -1){
        printf("system() call fails.\n");
        exit(1);
    }
	strcpy(buffer1, "./sdoutput/");
	strcat(buffer1, fname);
}

#ifdef SD_win
void sd_mv_files(char *fname)
{
    TCHAR path[BUFFER_SIZE];
    TCHAR inst_dir[BUFFER_SIZE];
    TCHAR buff[BUFFER_SIZE];
    char file_name[BUFFER_SIZE];
    TCHAR L_file_name[BUFFER_SIZE];
    
    if(!GetCurrentDirectory(BUFFER_SIZE, path))
        printf("GetCurrentDirectory() failed!\n");
    
    printf("Your current directory is: %S\n", path);
    
    wcscat(path,L"\\sdoutput");
    printf("Your SD Output directory is: %S\n", path);
    CreateDirectory(path,NULL);
    
    mbstowcs(inst_dir,fname,NAME_SIZE);
    printf("Your instance's name is: %S\n", inst_dir);
    wcscat(path,L"\\");
    wcscat(path,inst_dir);
    printf("Your instance SD Result directory is: %S\n", path);
    CreateDirectory(path,NULL);
    
  	wcscpy(buff, path);
	wcscat(buff, L"\\after_coef_change.lp");
	MoveFile(L"after_coef_change.lp", buff);
    
	wcscpy(buff, path);
	wcscat(buff, L"\\Batch_dual.out");
	MoveFile(L"Batch_dual.out", buff);
    
	wcscpy(buff, path);
	wcscat(buff, L"\\Batch_Obj.out");
	MoveFile(L"Batch_Obj.out", buff);
    
	wcscpy(buff, path);
	wcscat(buff, L"\\Batch_x.out");
	MoveFile(L"Batch_x.out", buff);
    
	wcscpy(buff, path);
	wcscat(buff, L"\\before_coef_change.lp");
	MoveFile(L"before_coef_change.lp", buff);
    
	wcscpy(buff, path);
	wcscat(buff, L"\\eval.dat");
	MoveFile(L"eval.dat", buff);
    
	wcscpy(buff, path);
	wcscat(buff, L"\\final_master.lp");
	MoveFile(L"final_master.lp", buff);
    
	wcscpy(buff, path);
	wcscat(buff, L"\\final-batch-prob.lp");
	MoveFile(L"final-batch-prob.lp", buff);
    
	wcscpy(buff, path);
	wcscat(buff, L"\\incumb.out");
	MoveFile(L"incumb.out", buff);
    
	wcscpy(buff, path);
	wcscat(buff, L"\\Master_Dual.out");
	MoveFile(L"Master_Dual.out", buff);
    
	wcscpy(buff, path);
	wcscat(buff, L"\\Master_Obj.out");
	MoveFile(L"Master_Obj.out", buff);
    
	wcscpy(buff, path);
	wcscat(buff, L"\\orig.lp");
	MoveFile(L"orig.lp", buff);
    
	wcscpy(buff, path);
	wcscat(buff, L"\\Summary.out");
	MoveFile(L"Summary.out", buff);
    
	wcscpy(buff, path);
	wcscat(buff, L"\\time_sample.out");
	MoveFile(L"time_sample.out", buff);
    
	strcpy(file_name, "\\");
	strcat(file_name, fname);
	strcat(file_name, ".detailed_rep_soln.out");
	mbstowcs(L_file_name,file_name,BUFFER_SIZE);
	wcscpy(buff, path);
	wcscat(buff, L_file_name);
	strcpy(file_name, fname);
	strcat(file_name, ".detailed_rep_soln.out");
	mbstowcs(L_file_name,file_name,BUFFER_SIZE);
	MoveFile(L_file_name, buff);
    
	strcpy(file_name, "\\");
	strcat(file_name, fname);
	strcat(file_name, ".detailed_soln.out");
	mbstowcs(L_file_name,file_name,BUFFER_SIZE);
	wcscpy(buff, path);
	wcscat(buff, L_file_name);
	strcpy(file_name, fname);
	strcat(file_name, ".detailed_soln.out");
	mbstowcs(L_file_name,file_name,BUFFER_SIZE);
	MoveFile(L_file_name, buff);
    
	strcpy(file_name, "\\");
	strcat(file_name, fname);
	strcat(file_name, ".obj.out");
	mbstowcs(L_file_name,file_name,BUFFER_SIZE);
	wcscpy(buff, path);
	wcscat(buff, L_file_name);
	strcpy(file_name, fname);
	strcat(file_name, ".obj.out");
	mbstowcs(L_file_name,file_name,BUFFER_SIZE);
	MoveFile(L_file_name, buff);
    
	strcpy(file_name, "\\");
	strcat(file_name, fname);
	strcat(file_name, ".time.out");
	mbstowcs(L_file_name,file_name,BUFFER_SIZE);
	wcscpy(buff, path);
	wcscat(buff, L_file_name);
	strcpy(file_name, fname);
	strcat(file_name, ".time.out");
	mbstowcs(L_file_name,file_name,BUFFER_SIZE);
	MoveFile(L_file_name, buff);
}
#else
void sd_mv_files(char *buffer1, char *buffer2, char *fname)
{
    int status;
    strcpy(buffer2, "mv *.out *.lp *.dat ");
	strcat(buffer2, buffer1);
	status = system(buffer2);
    if(status == -1){
        printf("system() call fails.\n");
        exit(1);
    }
}
#endif



