/************This file is the TOURNAMENT SELECTION OPERATOR.*********************** 
 It will copy the current population, which is the old population and 
 shuffles it to generate a temperate population, which is called midpop. 
 Then it will pair each individual in the old population with the corresponding 
 individual in the midpop and select the better one to put into another 
 population called matepop. This matepop is the parent population (the mating pool)
 on which crossover and mutation operators will perform.*/

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
	reallim[maxchrom2][2]; /*the lower and upper limts for real avriables*/
	
extern int popsize; /*population size*/


//function prototype
individual copyind(individual ind1);
void shuffle(population *pop_ptr);

void selectt(population *pop_ptr1, population *pop_ptr2, population *pop_ptr3)
{ //they are oldpop, midpop and matepop
	int i;

	
	/**************shuffle the temperate population - midpop(pop_ptr2)*****************/
	shuffle(pop_ptr2);
	/******************************midpop is shuffled*********************************/
	
	pop_ptr1->ind_ptr=&(pop_ptr1->ind[0]);
	pop_ptr2->ind_ptr=&(pop_ptr2->ind[0]);
	pop_ptr3->ind_ptr=&(pop_ptr3->ind[0]);
	
	/*Tournament Selection - pair oldpop with midpop and generate matepop (pop_ptr3)*/
	//Compare indivisuals by rank first:the smaller the rank, the better the individual
	for (i=0;i<popsize;i++) {
		if ((pop_ptr1->ind[i].rank)<(pop_ptr2->ind[i].rank)) {
			pop_ptr3->ind[i]=copyind(pop_ptr1->ind[i]);
		}
		else if ((pop_ptr1->ind[i].rank)>(pop_ptr2->ind[i].rank)) {
			pop_ptr3->ind[i]=copyind(pop_ptr2->ind[i]);
		}
		else { //when the ranks are the same, compare cub_len; larger cub_len better
			if ((pop_ptr1->ind[i].cub_len)<(pop_ptr2->ind[i].cub_len)) {
				pop_ptr3->ind[i]=copyind(pop_ptr2->ind[i]); 
			}
			else {
				pop_ptr3->ind[i]=copyind(pop_ptr1->ind[i]);
			}
		} //loop over two inds having the same rank ends
	} //loop over the tournament selection over population ends

	return;
}

