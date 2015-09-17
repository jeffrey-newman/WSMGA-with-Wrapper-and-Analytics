/*this file decode the integer chromosome in each individual to
get real values*/

#include "moga.h"

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

void decode(population *pop_ptr)
{
	int i,j,k;
	
	pop_ptr->ind_ptr=&(pop_ptr->ind[0]);
	
	for (i=0;i<popsize;i++) {//for i individuals
		if (no_intevar>0) { //there is integer variables
			for (j=0;j<no_intevar;j++) {//for jth integer DVs in the ith individual
				for (k=0;k<no_option[j];k++) {//for each option of jth DV in ith individual
					if (pop_ptr->ind[i].chrom1[j]==k+1) {
						pop_ptr->ind[i].xchrom1[j]=xoption[j][k];
						//printf("\nxchrom1 for %dth ind %d DV is %g", i+1, j+1,pop_ptr->ind[i].xchrom1[j]);
						pop_ptr->ind[i].ychrom1[j]=yoption[j][k];
						//printf("\nychrom1 for %dth ind %d DV is %g", i+1, j+1,pop_ptr->ind[i].ychrom1[j]);
						pop_ptr->ind[i].zchrom1[j]=zoption[j][k];
						//printf("\nzchrom1 for %dth ind %d DV is %g", i+1, j+1,pop_ptr->ind[i].zchrom1[j]);
						pop_ptr->ind[i].uchrom1[j]=uoption[j][k];
						//printf("\nzchrom1 for %dth ind %d DV is %g", i+1, j+1,pop_ptr->ind[i].uchrom1[j]);
						pop_ptr->ind[i].vchrom1[j]=voption[j][k];
						//printf("\nzchrom1 for %dth ind %d DV is %g", i+1, j+1,pop_ptr->ind[i].vchrom1[j]);
						//printf("\n---------------------------------------------------\n");
					}
				}
			}
		}
	}
	
	return;
}
