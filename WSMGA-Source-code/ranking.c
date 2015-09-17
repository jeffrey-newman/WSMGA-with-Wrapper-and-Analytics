/*this file ranks the individuals without constraints in each generation
 and also demarcate different Pareto Fronts*/

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

/*function prototype*/
int indcmp(double *ptr1, double *ptr2);

void ranking(population *pop_ptr)
{
	int i,j,k, /*counters*/
		rnk, /*rank*/
		val, /*value obtained after comparing two individuals*/
		nondom, /*no. of non-domonated members*/
		maxrank1, /*max rank of the population*/
		rankarr[maxpop], /*array sorting the  number of individuals at a rank*/
		q;//the  number of individuals at a rank
	double *ptr1,*ptr2;
	
	/********************************ranking****************************************/
	
	//initialize the ranks to zero
	rnk=0;
	nondom=0;
	maxrank1=0;
	
	//initialize all of the flags to 2
	for (i=0;i<popsize;i++)
		pop_ptr->ind[i].flag=2;
	q=0; 
	for (k=0;k<popsize;k++,q=0) {
		for (j=0;j<popsize;j++) {
			if (pop_ptr->ind[j].flag!=1)
				break; 
		}
		if (j==popsize)
			break; //break, when all individuals are assigned ranks
		rnk=rnk+1;
		for (j=0;j<popsize;j++) { //set the flags of domonated ind to 2
			if (pop_ptr->ind[j].flag==0)
				pop_ptr->ind[j].flag=2;
		}
		//select two individuals to compare in order to assign ranks
		for (i=0;i<popsize;i++) { 
			pop_ptr->ind_ptr=&(pop_ptr->ind[i]);//select the first individual
			if (pop_ptr->ind_ptr->flag!=1 && pop_ptr->ind_ptr->flag!=0) {//flag check
				ptr1=&(pop_ptr->ind_ptr->fit[0]);
				for (j=0;j<popsize;j++) { //select the second individual
					if (i!=j) {
						if (pop_ptr->ind[j].flag!=1) {
							pop_ptr->ind_ptr=&(pop_ptr->ind[j]);
							ptr2=&(pop_ptr->ind_ptr->fit[0]);
							val=indcmp(ptr1,ptr2); //compare the two individuals
							switch (val) {
								case 2: //first individual is dominated
									pop_ptr->ind[i].flag=0;
									break;
								case 1: //second individual is dominated
									pop_ptr->ind[j].flag=0;
									break;
								case 3:
									nondom++;
									if (pop_ptr->ind[j].flag!=0)
										pop_ptr->ind[j].flag=3;//i & j are non-dominated
									break;
								default:
									printf("\n!ranking has problems in comparing inds!");
									break;
							}
						}
					}
				} //loop over the second individual j ends
				if (j==popsize) {
					pop_ptr->ind[i].rank=rnk;
					if (pop_ptr->ind[i].flag!=0) {
						pop_ptr->ind[i].flag=1;//individual i is domonating
					rankarr[q]=i;
					q++;
					}
				}
			}//loop over flag check ends
		}//loop over first individual i ends
		pop_ptr->rankno[rnk-1]=q;
	}//loop over k ends
	maxrank1=rnk;
	
	/*find max rank of the population*/
	for (i=0;i<popsize;i++) {
		rnk=pop_ptr->ind[i].rank;
		if (rnk>maxrank1)
			maxrank1=rnk;
	}
	pop_ptr->maxrank=maxrank1; //assign the max rank value of the population

	return;
}
