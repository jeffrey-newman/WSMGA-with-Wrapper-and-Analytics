/*this file copy the first population to the second population, 
which makes a exact copy of the first population*/

#include "moga.h"

extern int no_gener, /*number of generations*/
    no_obj,  /*no. of objective functions*/
	no_cons, /*no. of constraints*/
	no_pro, /*no. of property*/
    no_option[maxchrom1],  /*list of no. of options for each integer decision variable*/
    no_intevar, /*no. of integer decision variables = size of chromosome1*/
	no_realvar, /*no. of real decision variables = size of chromosome2*/
	no_mut, /*no. of mutations happened*/
	no_cross, /*no. of crossovers happened*/
	no_cross_real, /*no. of real variables crossed*/
	ans;
	
extern double seed, /*random seed*/
	pc, /*crossover probability for both integer and real chromosomes*/
	pm, /*mutation probability for both interger and real chromosomes*/
	di, /*distribution index for SBX*/
	xoption[maxchrom1][maxno_option], //lookout table of 1st property of no_intevar
	yoption[maxchrom1][maxno_option], //lookout table of 2nd property of no_intevar
	zoption[maxchrom1][maxno_option], //lookout table of t3rd property of no_intevar
	uoption[maxchrom1][maxno_option],
	voption[maxchrom1][maxno_option],
	reallim[maxchrom2][2]; /*the lower and upper limts for real avriables*/

extern int popsize; /*population size*/

extern population oldpop, *old_pop_ptr,
		newpop, *new_pop_ptr,
	    matepop,  *mate_pop_ptr,
	    midpop, *mid_pop_ptr; //will be used in select.c



/*function prototype*/
individual copyind(individual ind1);

void copypop(population *pop_ptr1, population *pop_ptr2)
{
	int i;
	
	/**************generate the temperate population - midpop*********************/
	pop_ptr2->maxrank=pop_ptr1->maxrank;
	for (i=0;i<popsize;i++) { //copy the properties of each individual
		pop_ptr2->rankno[i]=pop_ptr1->rankno[i];
		pop_ptr2->ind[i]=copyind(pop_ptr1->ind[i]);
	}
		
	return;
}
