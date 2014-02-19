#include "prob.h"
#include "cell.h"
#include "soln.h"
#include "theta.h"
#include "testout.h"
#include "solver.h"
#include "sdglobal.h"

void print_contents(one_problem *test, char *filename)
{
	/* int i, n=1, block=50;*/
	FILE *out;

#ifdef TRACE
	printf("Inside print_contents\n");
#endif
#ifdef DEBUG
	printf("zl_prt_contnt filename = %s\n", filename);
#endif
	out = fopen(filename, "w");
	fprintf(out,
			"This file contains the contents of a test run for reading in\n");
	fprintf(out, "a core file in the MPS format\n");
	fprintf(out, "The problem name is: %s\n", test->name);
	fprintf(out, "The objective name is: %s\n", test->objname);
	fprintf(out, "The objective sense is: %d\n", test->objsen);
	fprintf(out, "The value of mar is %d\n", test->mar);
	fprintf(out, "The value of marsz is %d\n", test->marsz);
	fprintf(out, "The value of mac is %d\n", test->mac);
	fprintf(out, "The value of macsz is %d\n", test->macsz);
	fprintf(out, "The value of matsz is %d\n", test->matsz);
	fprintf(out, "The value of cstorsz is %d\n", test->cstorsz);
	fprintf(out, "The value of rstorsz is %d\n\n", test->rstorsz);

#ifdef TIM
	printf("printing matbeg and matcnt\n");
#endif
	/*
	 fprintf(out,"The contents of matbeg and matcnt are : \n");
	 for(i=0; i< test->macsz; i++) 
	 fprintf(out,"\t%d\t%d\n",test->matbeg[i], test->matcnt[i]);
	 fprintf(out,"\n\n");
	 */
#ifdef TIM
	printf("printing matval and matind\n");
#endif
	/*
	 fprintf(out,"The contents of matval and matind are : \n");
	 for(i=0; i< test->matsz; i++) 
	 fprintf(out,"\t%f\t%d\n", test->matval[i], test->matind[i]);
	 fprintf(out,"\n\n");
	 */
#ifdef TIM
	printf("printing cname objx, bdl, and bdu\n");
#endif
	/* zl
	 fprintf(out,"The contents of cname objx, bdl, and bdu are :\n");
	 for(i=0; i< test->macsz; i++) 
	 fprintf(out,"\t%s\t%f\t%f\t%f\n", 
	 test->cname[i], test->objx[i], test->bdl[i], test->bdu[i]);
	 fprintf(out,"\n\n");

	 fprintf(out,"The contents of cstore are:\n");
	 i = 0;
	 while(i< test->cstorsz)
	 {
	 for(;i< test->cstorsz && i< (n*block); i++)
	 {
	 if(test->cstore[i] == '\0')
	 {
	 putc('\\',out);
	 putc('0',out);
	 }
	 else
	 putc(test->cstore[i], out);
	 }
	 n++;
	 fprintf(out,"\n");
	 }  
	 n = 1;
	 fprintf(out,"\n\n");
	 */
#ifdef DEBUG
	printf("zl_prt_contnt ~1\n");
#endif

#ifdef TIM
	printf("printing sense, rname, and rhs\n");
#endif
	/* zl
	 fprintf(out,"The contents of sense, rname, and rhs are : \n");
	 for(i=0; i< test->marsz; i++) 
	 fprintf(out,"\t%c\t%s\t%f\n", test->senx[i], test->rname[i], test->rhsx[i]);
	 fprintf(out,"\n\n");
	 */
#ifdef TIM
	printf("printing rstore\n");
#endif
	/* zl
	 fprintf(out,"The contents of rstore are:\n");
	 i = 0;
	 while(i< test->rstorsz)
	 {
	 for(;i< test->rstorsz && i< (n*block); i++)
	 {
	 if(test->rstore[i] == '\0')
	 {
	 putc('\\',out);
	 putc('0',out);
	 }
	 else
	 putc(test->rstore[i], out);
	 }
	 n++;
	 fprintf(out,"\n");
	 }  
	 n = 1;
	 fprintf(out,"\n\n");
	 */
	fclose(out);

}

void cplex_err_msg(sdglobal_type* sd_global, char *string, prob_type *p,
		cell_type *c, soln_type *s)
{
	int cnt;
	double height, eta_coef;

	printf("CPLEX error occured in %s\n", string);
	print_problem(c->master, "master_err.mps");
	print_problem(c->subprob, "subprob_err.mps");

	/* Print the height of all the cuts at the incumbent */
	for (cnt = 0; cnt < c->cuts->cnt; cnt++)
	{
		height = cut_height(sd_global, c->cuts->val[cnt], s->incumb_x, c,
				p->num);
		eta_coef = c->cuts->val[cnt]->cut_obs / c->k;
		printf("Cut #%d in row#%d: height=%lf, eta-coef=%lf", cnt,
				c->cuts->val[cnt]->row_num, height, eta_coef);
	}
}
