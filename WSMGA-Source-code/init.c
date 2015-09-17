/*initialize population*/

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
	ans; //if limits on real variables are rigid, ans=1
	
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


/*declare function prototypes used in this file*/
double randomperc(void);
int rnd(int low, int high);

void init(population *pop_ptr)
{
	int i,j,ran1;
	double d, d1;
	
	for (i=0;i<popsize;i++) { /*for the whole population*/
		pop_ptr->ind[i].tag=0;
		pop_ptr->ind[i].eval=0;
		pop_ptr->ind[i].gen=0;
		if (no_intevar!=0) { /*there is integer variables*/
			for (j=0;j<no_intevar;j++) { /*initialize the first chromosome*/
				ran1=rnd(1, no_option[j]);
				//printf("\n rand()= %d", ran1);
				pop_ptr->ind[i].chrom1[j]=ran1;//form 1 to no_option
			}
		}  /*chromosome1 is finished*/
		else { //when there is not integer DV
			for (j=0;j<maxchrom1;j++) {
				pop_ptr->ind[i].chrom1[j]=0;
				xoption[i][j]=0;
				yoption[i][j]=0;
				zoption[i][j]=0;
				uoption[i][j]=0;
				voption[i][j]=0;
			}
		}
		if (no_realvar!=0) {  /*initialize the second chromosome*/
			for (j=0;j<no_realvar;j++) {
				d=randomperc();
				//printf("\nd for j=%d is %f\n", j, d);
				d1=2*d-1; 
				//printf("\nd1 for j=%d is %f\n", j, d1);
				/*if no limits, generates any no. between 0 and infinity*/
				if (ans!=1)
					pop_ptr->ind[i].chrom2[j]=1/d1;
				/*if limits specified, generates no. between lower and upper limits
				as x = random no.*(upper_limit - lower_limit)+lower_limit*/
				else
					pop_ptr->ind[i].chrom2[j]=d*(reallim[j][1]-reallim[j][0])+reallim[j][0];
			}  /*chromosome2 is finished*/
		}  /*one individual is finished*/
		else { //when there is no real DV
			for (j=0;j<maxchrom2;j++)
				pop_ptr->ind[i].chrom2[j]=0.0;
		}
	}  /*the population is finished*/
	
	return;
}
