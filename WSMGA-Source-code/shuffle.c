/*this file contains a operator, which shuffles a population randomly*/

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

extern population oldpop, *old_pop_ptr,
		newpop, *new_pop_ptr,
	    matepop,  *mate_pop_ptr,
	    midpop, *mid_pop_ptr; //will be used in select.c

/*function prototype*/
int rnd(int low, int high);

void shuffle(population *pop_ptr)
{
	int i,j,m,n;
	int temp_chrom1[maxchrom1],
	    temp_rank,
		temp_flag,
		temp_tag,
		temp_eval,
		temp_gen;
	double temp_cub_len,
	    temp_overallcons,
		temp_fit[maxobj],
		temp_cons[maxcons],
		temp_property[maxpro];
	double temp_chrom2[maxchrom2];
	double temp_xchrom1[maxchrom1],
		temp_ychrom1[maxchrom1],
		temp_zchrom1[maxchrom1],
		temp_uchrom1[maxchrom1],
		temp_vchrom1[maxchrom1];
		
	pop_ptr->ind_ptr=&(pop_ptr->ind[0]); 
	for (i=0;i<popsize;i++) { /*shuffle midpop by swaping ind[n]&ind[m] for popsize times*/
		n=rnd(1,popsize)-1;
		m=rnd(1,popsize)-1;
		for (j=0;j<no_intevar;j++) { //shuffle chrom1, xchrom1, ychrom1, zchrom1
			temp_chrom1[j]=pop_ptr->ind[n].chrom1[j];
			pop_ptr->ind[n].chrom1[j]=pop_ptr->ind[m].chrom1[j];
			pop_ptr->ind[m].chrom1[j]=temp_chrom1[j];
			temp_xchrom1[j]=pop_ptr->ind[n].xchrom1[j];
			pop_ptr->ind[n].xchrom1[j]=pop_ptr->ind[m].xchrom1[j];
			pop_ptr->ind[m].xchrom1[j]=temp_xchrom1[j];
			temp_ychrom1[j]=pop_ptr->ind[n].ychrom1[j];
			pop_ptr->ind[n].ychrom1[j]=pop_ptr->ind[m].ychrom1[j];
			pop_ptr->ind[m].ychrom1[j]=temp_ychrom1[j];
			temp_zchrom1[j]=pop_ptr->ind[n].zchrom1[j];
			pop_ptr->ind[n].zchrom1[j]=pop_ptr->ind[m].zchrom1[j];
			pop_ptr->ind[m].zchrom1[j]=temp_zchrom1[j];
			temp_uchrom1[j]=pop_ptr->ind[n].uchrom1[j];
			pop_ptr->ind[n].uchrom1[j]=pop_ptr->ind[m].uchrom1[j];
			pop_ptr->ind[m].uchrom1[j]=temp_uchrom1[j];
			temp_vchrom1[j]=pop_ptr->ind[n].vchrom1[j];
			pop_ptr->ind[n].vchrom1[j]=pop_ptr->ind[m].vchrom1[j];
			pop_ptr->ind[m].vchrom1[j]=temp_vchrom1[j];
		}
		for (j=0;j<no_obj;j++) { //shuffle fitness
			temp_fit[j]=pop_ptr->ind[n].fit[j];
			pop_ptr->ind[n].fit[j]=pop_ptr->ind[m].fit[j];
			pop_ptr->ind[m].fit[j]=temp_fit[j];
		}
		for (j=0;j<no_cons;j++) { //shuffle cons
			temp_cons[j]=pop_ptr->ind[n].cons[j];
			pop_ptr->ind[n].cons[j]=pop_ptr->ind[m].cons[j];
			pop_ptr->ind[m].cons[j]=temp_cons[j];
		}
		for (j=0;j<no_pro;j++) { //shuffle properties
			temp_property[j]=pop_ptr->ind[n].property[j];
			pop_ptr->ind[n].property[j]=pop_ptr->ind[m].property[j];
			pop_ptr->ind[m].property[j]=temp_property[j];
		}
		for (j=0;j<no_realvar;j++) { //shuffle chrom2
			temp_chrom2[j]=pop_ptr->ind[n].chrom2[j];
			pop_ptr->ind[n].chrom2[j]=pop_ptr->ind[m].chrom2[j];
			pop_ptr->ind[m].chrom2[j]=temp_chrom2[j];
		}
		//shuffle the tag,eval,gen,rank, flag, cub_len & overallcons
		temp_tag=pop_ptr->ind[n].tag;
		pop_ptr->ind[n].tag=pop_ptr->ind[m].tag;
		pop_ptr->ind[m].tag=temp_tag;
		temp_eval=pop_ptr->ind[n].eval;
		pop_ptr->ind[n].eval=pop_ptr->ind[m].eval;
		pop_ptr->ind[m].eval=temp_eval;
		temp_gen=pop_ptr->ind[n].gen;
		pop_ptr->ind[n].gen=pop_ptr->ind[m].gen;
		pop_ptr->ind[m].gen=temp_gen;
		temp_rank=pop_ptr->ind[n].rank;
		pop_ptr->ind[n].rank=pop_ptr->ind[m].rank;
		pop_ptr->ind[m].rank=temp_rank;
		temp_flag=pop_ptr->ind[n].flag;
		pop_ptr->ind[n].flag=pop_ptr->ind[m].flag;
		pop_ptr->ind[m].flag=temp_flag;
		temp_cub_len=pop_ptr->ind[n].cub_len;
		pop_ptr->ind[n].cub_len=pop_ptr->ind[m].cub_len;
		pop_ptr->ind[m].cub_len=temp_cub_len;
		temp_overallcons=pop_ptr->ind[n].overallcons;
		pop_ptr->ind[n].overallcons=pop_ptr->ind[m].overallcons;
		pop_ptr->ind[m].overallcons=temp_overallcons;
	} 
	return;
}
