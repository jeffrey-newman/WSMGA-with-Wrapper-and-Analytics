/*this file ranks the individuals with constraints in each generation
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

void rankcon(population *pop_ptr)
{
	int i,j,k, /*counters*/
		rnk, /*rank*/
		val, /*value obtained after comparing two individuals*/
		nondom, /*no. of non-domonated members*/
		maxrank1, /*max rank of the population*/
		rankarr[maxpop], /*array sorting the  number of individuals at a rank*/
		q;//the  number of individuals at a rank
	double *ptr1,*ptr2,*err_ptr1,*err_ptr2;
	
	//double min_fit; /*the min fitness of a better rank*/
      //delta_fit; /*min difference of the crowdign distance between two ranks*/
	
	/********************************ranking****************************************/
	
	//initialize the ranks to zero
	rnk=0;
	nondom=0;
	maxrank1=0;
	
	//initialize all of the flags to 2
	for (i=0;i<popsize;i++)
		pop_ptr->ind[i].flag=2;
	q=0; 
	for (k=0;k<popsize;k++,q=0) {//start ranking each ind in the population
		for (j=0;j<popsize;j++) {
			if (pop_ptr->ind[j].flag!=1)
				break; //break, when all individuals are assigned ranks
		}
		if (j==popsize)
			break;
		rnk=rnk+1;
		for (j=0;j<popsize;j++) { //set the flags of domonated ind to 2
			if (pop_ptr->ind[j].flag==0)
				pop_ptr->ind[j].flag=2;
		}
		//select two individuals to compare in order to assign ranks
		for (i=0;i<popsize;i++) {//select the first individual i
			pop_ptr->ind_ptr=&(pop_ptr->ind[i]);
			if (pop_ptr->ind_ptr->flag!=1 && pop_ptr->ind_ptr->flag!=0) {//flag check
				ptr1=&(pop_ptr->ind_ptr->fit[0]);
				err_ptr1=&(pop_ptr->ind_ptr->overallcons);
				for (j=0;j<popsize;j++) {//select the second individual j
					if (i!=j) {
						if (pop_ptr->ind[j].flag!=1) {//if ind j is not dominating
							pop_ptr->ind_ptr=&(pop_ptr->ind[j]);
							ptr2=&(pop_ptr->ind_ptr->fit[0]);
							err_ptr2=&(pop_ptr->ind_ptr->overallcons);
							if (*err_ptr1>0.0 && *err_ptr2>0.0) {//both infeasible
								if (*err_ptr1<*err_ptr2) //when i is better
									pop_ptr->ind[j].flag=0; //j is dominated
								else {
									if (*err_ptr1>*err_ptr2) {//j is better
										pop_ptr->ind[i].flag=0; //i is dominated
									}
									else {
										nondom++;
										if (pop_ptr->ind[j].flag!=0)
											pop_ptr->ind[j].flag=3;
									}
								}
							}
							else { //not both infeasible
								if (*err_ptr1>0.0 && *err_ptr2<=0.0) {
									pop_ptr->ind[i].flag=0; //i is infeasible & j is feasible
									break;
								}
								else {
									if (*err_ptr1<=0.0 && *err_ptr2>0.0) {
										pop_ptr->ind[j].flag=0;//i is feasible & j is infeasible
									}
									else { //both feasible
										val=indcmp(ptr1,ptr2);
										switch (val) {
											case 2:
												pop_ptr->ind[i].flag=0;
												break;//i is dominated
											case 1:
												pop_ptr->ind[j].flag=0;
												break;//j is dominated
											case 3://i and j are not dominated by each other
												nondom++;
												if (pop_ptr->ind[j].flag!=0)
													pop_ptr->ind[j].flag=3;
												break;
											default:
												printf("!rankcon has problems in comparing inds!\n");
												break;
										}
									} //both feasible
								}
							} //not both feasible
						} //j is not dominating
					}//i!=j loop ends
				}//loop over j ends
				if (j==popsize) {
					pop_ptr->ind[i].rank=rnk;
					if (pop_ptr->ind[i].flag!=0) {
						pop_ptr->ind[i].flag=1;//i is dominating
						rankarr[q]=i;
						q++;
					}
				}
			}//loop over flag check ends
		}//loop over i ends
		pop_ptr->rankno[rnk-1]=q;
	}
	maxrank1=rnk;
	
	//find the max rank of the population
	for (i=0;i<popsize;i++) {
		rnk=pop_ptr->ind[i].rank;
		if (rnk>maxrank1)
			maxrank1=rnk;
	}
	
	pop_ptr->maxrank=maxrank1;

	return;
}

