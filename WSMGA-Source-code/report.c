/*This file is used to print the reports*/

# include "moga.h"

extern int no_gener, /*number of generations*/
    no_obj,  /*no. of objective functions*/
	no_cons, /*no. of constraints*/
	no_pro, /*no. of property*/
    no_option[maxchrom1],  /*list of no. of options for each integer decision variable*/
    no_intevar, /*no. of integer decision variables = size of chromosome1*/
	no_realvar, /*no. of real decision variables = size of chromosome2*/
	no_mut1, /*no. of mutations happened to chrom1*/
	no_mut2, /*no. of mutations happened to chrom2*/
	no_cross, /*no. of crossovers happened*/
	no_cross_real, /*no. of real variables crossed*/
	ans;
	
extern double seed, /*random seed*/
	pc, /*crossover probability for both integer and real chromosomes*/
	pm, /*mutation probability for both interger and real chromosomes*/
	di, /*distribution index for SBX*/
	dim, /*distribution index for polynomial mutation*/
	xoption[maxchrom1][maxno_option], //lookout table of 1st property of no_intevar
	yoption[maxchrom1][maxno_option], //lookout table of 2nd property of no_intevar
	zoption[maxchrom1][maxno_option], //lookout table of t3rd property of no_intevar
	uoption[maxchrom1][maxno_option],
	voption[maxchrom1][maxno_option],
	reallim[maxchrom2][2]; /*the lower and upper limts for real avriables*/
	
extern int popsize; /*population size*/

void report(int t,population *pop_ptr1, population *pop_ptr2, FILE *lastit)
{ //t=generation no. i in moga.c; best_pop_ptr; mate_pop_ptr;  lastit;  
	//                                                          plot.out;
	int i, j; //counters
	double *fit_ptr1, //pointer of fitness of inds in pop1
		*fit_ptr2, //pointer of fitness of inds in pop2
		*chrom2_ptr1, //pointer of chrom2(real var) of inds in pop1
	    *chrom2_ptr2; //pointer of chrom2 of inds in pop2
	double *cons_ptr1, //pointer of constraints of inds in pop1(oldpop)
	    *cons_ptr2; //pointer of constraints of inds in pop2(matepop)
	double *pro_ptr1, //pointer of properties of inds in pop1(oldpop)
	    *pro_ptr2; //pointer of properties of inds in pop2(matepop)
		//*err2; //pointer of overallconstraints of inds in pop2(matepop)
	int *chrom1_ptr1, //pointer of chrom1(integer var) of inds in pop1
	    *chrom1_ptr2; //pointer of chrom1 of inds, in pop2
	
		
	if (t==no_gener-1) {
		pop_ptr1->ind_ptr=&(pop_ptr1->ind[0]); 
		pop_ptr2->ind_ptr=&(pop_ptr2->ind[0]); //ind[0] of both oldpop and matepop
		fprintf(lastit, " variables(real %d; integer %d)  fitness(%d)  constraints(%d) overallcons  properties(%d) "
			"rank  cub_len  evaluation_number  generation_number\n", no_realvar,no_intevar,no_obj,no_cons, no_pro);
		for (i=0;i<popsize;i++) {
			chrom1_ptr2=&(pop_ptr2->ind_ptr->chrom1[0]);
			chrom2_ptr2=&(pop_ptr2->ind_ptr->chrom2[0]);
			fit_ptr2=&(pop_ptr2->ind_ptr->fit[0]);
			cons_ptr2=&(pop_ptr2->ind_ptr->cons[0]);
			pro_ptr2=&(pop_ptr2->ind_ptr->property[0]);
			chrom1_ptr1=&(pop_ptr1->ind_ptr->chrom1[0]);
			chrom2_ptr1=&(pop_ptr1->ind_ptr->chrom2[0]);
			fit_ptr1=&(pop_ptr1->ind_ptr->fit[0]);
			cons_ptr1=&(pop_ptr1->ind_ptr->cons[0]);
			pro_ptr1=&(pop_ptr1->ind_ptr->property[0]);
			for (j=0;j<no_realvar;j++) {
				fprintf(lastit, "%.4f ", *chrom2_ptr1++);
			}
			for (j=0;j<no_intevar;j++) {
				fprintf(lastit, "%d ", *chrom1_ptr1++);
			}
			for (j=0;j<no_obj;j++) {
				fprintf(lastit, "  %.4f ", *fit_ptr1++);
			}
			if (no_cons!=0) {
				for (j=0;j<no_cons;j++) {
					fprintf(lastit, "  %g", *cons_ptr1++);
				}
				fprintf(lastit, "  %.2e", pop_ptr1->ind_ptr->overallcons);
			}
			if (no_pro!=0) {
				for (j=0;j<no_pro;j++) {
					fprintf(lastit, "  %g", *pro_ptr1++);
				}
			}
			fprintf(lastit, " %d", pop_ptr1->ind_ptr->rank);
			fprintf(lastit, " %.2e ", pop_ptr1->ind_ptr->cub_len);
			fprintf(lastit, " %d", pop_ptr1->ind_ptr->eval);
			fprintf(lastit, " %d", pop_ptr1->ind_ptr->gen);
			fprintf(lastit, "\n");
			pop_ptr1->ind_ptr=&(pop_ptr1->ind[i+1]);
			pop_ptr2->ind_ptr=&(pop_ptr2->ind[i+1]);
		}
		
	}
	
	pop_ptr1->ind_ptr=&(pop_ptr1->ind[0]); //locate the start of the pointers to 
	pop_ptr2->ind_ptr=&(pop_ptr2->ind[0]); //ind[0] of both bestpop and matepop
	
	return;
}
