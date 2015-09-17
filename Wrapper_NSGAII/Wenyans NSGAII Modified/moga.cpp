/*this is a multi-objective GA, which can take integer and real decision variable values
 this program is based on NSGA-II developed by Prof. Deb and his students.

The program generates a numbwer of output files:
1. output.out
   This file has the detailed record for all variables, fitness values, constraint
violation values,overall constraint violation

*/

#include "moga.h"

int no_gener, /*number of generations*/
    no_obj,  /*no. of objective functions*/
	no_cons, /*no. of constraints*/
	no_pro, /*no. of property*/
	tot_option, //no. of total options
    no_option[maxchrom1],  /*list of no. of options for each integer decision variable*/
    no_intevar, /*no. of integer decision variables = size of chromosome1*/
	no_realvar, /*no. of real decision variables = size of chromosome2*/
	no_mut1, /*no. of mutations happened to chrom1*/
	no_mut2, /*no. of mutations happened to chrom2*/
	no_cross, /*no. of crossovers happened*/
	no_cross_real, /*no. of real variables crossed*/
	ans; //if limits on real variables are rigid, ans=1
	
double seed, /*random seed*/
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

int popsize; /*population size*/

int noeval; //no. of evaluations

population oldpop, *old_pop_ptr,
		newpop, *new_pop_ptr,
	    matepop,  *mate_pop_ptr,
		bestpop, *best_pop_ptr,
	    midpop, *mid_pop_ptr; //will be used in select.c


/*define function prototypes*/	
void input(FILE *rep_ptr);	
void randomize(void);
void init(population *pop_prt);
void decode(population *pop_ptr);
void selectt(population *pop_ptr1, population *pop_ptr2, population *pop_ptr3);
void crossover(population *pop_ptr1, population *pop_ptr2);
void mutation(population *new_pop_ptr);
void evaluatepop(population *pop_ptr, int gen); //update fit and cons of changed inds in a pop
void ranking(population *pop_ptr);
void rankcon(population *pop_ptr);
void keepalive(population *pop_ptr1, population *pop_ptr2, population *pop_ptr3, int gen);
void copypop(population *pop_ptr1, population *pop_ptr2);
individual copyind(individual ind1);
void selectbest(population *pop_ptr1, population *pop_ptr2, int t);
void mogareport(int t,population *pop_ptr1, population *pop_ptr2, FILE *lastit);


/*define other small functions*/
void errors(int value);

int main(int argc, char *argv[])
{
	/*define local variables*/
	clock_t start, end, t1, t2, t3, t4;
	int i,g; //counters
	    //maxrank1; //the larger maxrank between oldpop and matepop
	//double tot; //sum of no. of inds in a rank in both oldpop and newpop
	FILE *rep_ptr, *lastit;/*File Pointers*/
	
    if( argc == 2 )
		printf("The argument supplied is %s\n", argv[1]);
	else if( argc > 2 )
		printf("Too many arguments supplied.\n");
	else
		printf("One argument expected.\n");		
	seed = (double)atof(argv[1]);
	if (seed<=0.0 || seed>=1.0)
    {
        printf("\n Entered seed value is wrong, seed value must be in (0,1) \n");
        exit(11);
    }
  /*open files*/
	rep_ptr=fopen ("output.out","w");
	lastit = fopen("plot.out","w");
	
	old_pop_ptr=&(oldpop);
	
	no_mut1=0;
	no_mut2=0;
	no_cross = 0;
	
	input(rep_ptr);  /*obtain inputs*/
	randomize(); /*initialize random no. generator*/
	
  /*initialize population*/
	init(old_pop_ptr);
	//debug
	/*printf("\n The initial population are:");
	printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
	for (i=0;i<popsize;i++) {
		printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			old_pop_ptr->ind[i].chrom1[0],old_pop_ptr->ind[i].chrom1[1],
			old_pop_ptr->ind[i].chrom2[0],old_pop_ptr->ind[i].chrom2[1],
			old_pop_ptr->ind[i].fit[0],old_pop_ptr->ind[i].fit[1],
			old_pop_ptr->ind[i].cons[0],old_pop_ptr->ind[i].cons[1],
			old_pop_ptr->ind[i].overallcons,old_pop_ptr->ind[i].cub_len,
			old_pop_ptr->ind[i].rank,old_pop_ptr->ind[i].tag,
			old_pop_ptr->ind[i].eval,old_pop_ptr->ind[i].gen);
	}*/

  /*decode integer chromosomes*/
	if (no_intevar>0) {
		decode(old_pop_ptr);
	}
	//debug
	/*printf("\n The decoded initial population are:");
	printf("\n DV1 DV2   D_DV1x   D_DV1y   D_DV1z   D_DV1u   D_DV1v   D_DV1x   D_DV1y   D_DV1z   D_DV1u   D_DV1v");
	for (i=0;i<popsize;i++) {
		printf("\n %d   %d  %.4f    %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f ",
			old_pop_ptr->ind[i].chrom1[0],old_pop_ptr->ind[i].chrom1[1],
			old_pop_ptr->ind[i].xchrom1[0],old_pop_ptr->ind[i].ychrom1[0],
			old_pop_ptr->ind[i].zchrom1[0],old_pop_ptr->ind[i].uchrom1[0],
			old_pop_ptr->ind[i].vchrom1[0],old_pop_ptr->ind[i].xchrom1[1],
			old_pop_ptr->ind[i].ychrom1[1],old_pop_ptr->ind[i].zchrom1[1],
			old_pop_ptr->ind[i].uchrom1[1],old_pop_ptr->ind[i].vchrom1[1]);
	}*/
  
  /*initialize the rank array having different individuals at
  a particular rank to zero*/
	new_pop_ptr=&(newpop);
	for (i=0;i<popsize;i++) {
		old_pop_ptr->rankno[i]=0;
		new_pop_ptr->rankno[i]=0;
		old_pop_ptr->ind[i].tag=1;
	}
  
  /*evaluation the objective and constraint functions of the current population*/
	noeval=0;
	//assig generation number to each ind in matepop
	evaluatepop(old_pop_ptr, -1);
	for (i=0;i<popsize;i++) {
		old_pop_ptr->ind[i].tag=0;
	}
	//debug
	/*printf("\n After evaluation the initial oldpop are:");
	printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
	for (i=0;i<popsize;i++) {
		printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			old_pop_ptr->ind[i].chrom1[0],old_pop_ptr->ind[i].chrom1[1],
			old_pop_ptr->ind[i].chrom2[0],old_pop_ptr->ind[i].chrom2[1],
			old_pop_ptr->ind[i].fit[0],old_pop_ptr->ind[i].fit[1],
			old_pop_ptr->ind[i].cons[0],old_pop_ptr->ind[i].cons[1],
			old_pop_ptr->ind[i].overallcons,old_pop_ptr->ind[i].cub_len,
			old_pop_ptr->ind[i].rank,old_pop_ptr->ind[i].tag,
			old_pop_ptr->ind[i].eval,old_pop_ptr->ind[i].gen);
	}*/
	//rank the oldpop before copy it into midpop and selection
	if (no_cons==0) {
		old_pop_ptr->ind_ptr->overallcons=0.0;
		ranking(old_pop_ptr);
	}
	else {
		rankcon(old_pop_ptr);
	}
	//debug
	/*printf("\n After ranking/rankcon the initial oldpop are:");
	printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
	for (i=0;i<popsize;i++) {
		printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			old_pop_ptr->ind[i].chrom1[0],old_pop_ptr->ind[i].chrom1[1],
			old_pop_ptr->ind[i].chrom2[0],old_pop_ptr->ind[i].chrom2[1],
			old_pop_ptr->ind[i].fit[0],old_pop_ptr->ind[i].fit[1],
			old_pop_ptr->ind[i].cons[0],old_pop_ptr->ind[i].cons[1],
			old_pop_ptr->ind[i].overallcons,old_pop_ptr->ind[i].cub_len,
			old_pop_ptr->ind[i].rank,old_pop_ptr->ind[i].tag,
			old_pop_ptr->ind[i].eval,old_pop_ptr->ind[i].gen);
	}*/
	
	
  /***********************************************************************************/
  /****************************GNENRATION LOOP STARTS*********************************/
	
	for (g=0;g<no_gener;g++) {
		printf("\nGeneration number = %d\n", g+1);
		start=clock();
		//fprintf(rep_ptr,"Population at Generation No. -->%d\n",g+1);
		//fprintf(rep_ptr,"#Generation No. -->%d\n",g+1);
		//fprintf(rep_ptr,"#Variable_vector  Fitness_vector  Constraint_violation  Overall_violation\n");
		
		old_pop_ptr=&(oldpop);
		mid_pop_ptr=&(midpop);
		mate_pop_ptr=&(matepop);
		new_pop_ptr=&(newpop);
		best_pop_ptr=&(bestpop);
		
	 /**********************Two-member Tournament Selection***************************/
		//copy oldpop into midpop-midpop is the exact copy of oldpop
		copypop(old_pop_ptr, mid_pop_ptr);
		
		//debug
		/*printf("\n ===================Copypop==========================\n");
		printf("\n In generation loop:After copypop the initial oldpop are:");
		printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
		for (i=0;i<popsize;i++) {
			printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			old_pop_ptr->ind[i].chrom1[0],old_pop_ptr->ind[i].chrom1[1],
			old_pop_ptr->ind[i].chrom2[0],old_pop_ptr->ind[i].chrom2[1],
			old_pop_ptr->ind[i].fit[0],old_pop_ptr->ind[i].fit[1],
			old_pop_ptr->ind[i].cons[0],old_pop_ptr->ind[i].cons[1],
			old_pop_ptr->ind[i].overallcons,old_pop_ptr->ind[i].cub_len,
			old_pop_ptr->ind[i].rank,old_pop_ptr->ind[i].tag,
			old_pop_ptr->ind[i].eval,old_pop_ptr->ind[i].gen);
		}
		printf("\n In generation loop:After copypop the midpop are:");
		printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
		for (i=0;i<popsize;i++) {
			printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			mid_pop_ptr->ind[i].chrom1[0],mid_pop_ptr->ind[i].chrom1[1],
			mid_pop_ptr->ind[i].chrom2[0],mid_pop_ptr->ind[i].chrom2[1],
			mid_pop_ptr->ind[i].fit[0],mid_pop_ptr->ind[i].fit[1],
			mid_pop_ptr->ind[i].cons[0],mid_pop_ptr->ind[i].cons[1],
			mid_pop_ptr->ind[i].overallcons,mid_pop_ptr->ind[i].cub_len,
			mid_pop_ptr->ind[i].rank,mid_pop_ptr->ind[i].tag,
			mid_pop_ptr->ind[i].eval,mid_pop_ptr->ind[i].gen);
		}*/
		
		//shuffle midpop and select from oldpop+midpop to form matepop
		selectt(old_pop_ptr, mid_pop_ptr, mate_pop_ptr);
		
		//debug
		/*printf("\n===================Select===========================\n");
		printf("\n In generation loop:After selection the initial oldpop are:");
		printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
		for (i=0;i<popsize;i++) {
			printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			old_pop_ptr->ind[i].chrom1[0],old_pop_ptr->ind[i].chrom1[1],
			old_pop_ptr->ind[i].chrom2[0],old_pop_ptr->ind[i].chrom2[1],
			old_pop_ptr->ind[i].fit[0],old_pop_ptr->ind[i].fit[1],
			old_pop_ptr->ind[i].cons[0],old_pop_ptr->ind[i].cons[1],
			old_pop_ptr->ind[i].overallcons,old_pop_ptr->ind[i].cub_len,
			old_pop_ptr->ind[i].rank,old_pop_ptr->ind[i].tag,
			old_pop_ptr->ind[i].eval,old_pop_ptr->ind[i].gen);
		}
		printf("\n In generation loop:After selection the midpop are:");
		printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
		for (i=0;i<popsize;i++) {
			printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			mid_pop_ptr->ind[i].chrom1[0],mid_pop_ptr->ind[i].chrom1[1],
			mid_pop_ptr->ind[i].chrom2[0],mid_pop_ptr->ind[i].chrom2[1],
			mid_pop_ptr->ind[i].fit[0],mid_pop_ptr->ind[i].fit[1],
			mid_pop_ptr->ind[i].cons[0],mid_pop_ptr->ind[i].cons[1],
			mid_pop_ptr->ind[i].overallcons,mid_pop_ptr->ind[i].cub_len,
			mid_pop_ptr->ind[i].rank,mid_pop_ptr->ind[i].tag,
			mid_pop_ptr->ind[i].eval,mid_pop_ptr->ind[i].gen);
		}
		printf("\n In generation loop:After selection the matepop are:");
		printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
		for (i=0;i<popsize;i++) {
			printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			mate_pop_ptr->ind[i].chrom1[0],mate_pop_ptr->ind[i].chrom1[1],
			mate_pop_ptr->ind[i].chrom2[0],mate_pop_ptr->ind[i].chrom2[1],
			mate_pop_ptr->ind[i].fit[0],mate_pop_ptr->ind[i].fit[1],
			mate_pop_ptr->ind[i].cons[0],mate_pop_ptr->ind[i].cons[1],
			mate_pop_ptr->ind[i].overallcons,mate_pop_ptr->ind[i].cub_len,
			mate_pop_ptr->ind[i].rank,mate_pop_ptr->ind[i].tag,
			mate_pop_ptr->ind[i].eval,mate_pop_ptr->ind[i].gen);
		}*/
	  
	 /***********************************Crossover************************************/
	 //one-point crossover for chrom1 and SBX for chrom2
	 //before crossover and mutation set tag to zero
		for (i=0;i<popsize;i++) {
			new_pop_ptr->ind[i].tag=0;
		}
		crossover(new_pop_ptr, mate_pop_ptr);
		
		//debug
		/*printf("\n===================Crossover===========================\n");
		printf("\n In generation loop:After crossover matepop, the initial oldpop are:");
		printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
		for (i=0;i<popsize;i++) {
			printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			old_pop_ptr->ind[i].chrom1[0],old_pop_ptr->ind[i].chrom1[1],
			old_pop_ptr->ind[i].chrom2[0],old_pop_ptr->ind[i].chrom2[1],
			old_pop_ptr->ind[i].fit[0],old_pop_ptr->ind[i].fit[1],
			old_pop_ptr->ind[i].cons[0],old_pop_ptr->ind[i].cons[1],
			old_pop_ptr->ind[i].overallcons,old_pop_ptr->ind[i].cub_len,
			old_pop_ptr->ind[i].rank,old_pop_ptr->ind[i].tag,
			old_pop_ptr->ind[i].eval,old_pop_ptr->ind[i].gen);
		}
		printf("\n In generation loop:After crossover matepop, the midpop are:");
		printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
		for (i=0;i<popsize;i++) {
			printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			mid_pop_ptr->ind[i].chrom1[0],mid_pop_ptr->ind[i].chrom1[1],
			mid_pop_ptr->ind[i].chrom2[0],mid_pop_ptr->ind[i].chrom2[1],
			mid_pop_ptr->ind[i].fit[0],mid_pop_ptr->ind[i].fit[1],
			mid_pop_ptr->ind[i].cons[0],mid_pop_ptr->ind[i].cons[1],
			mid_pop_ptr->ind[i].overallcons,mid_pop_ptr->ind[i].cub_len,
			mid_pop_ptr->ind[i].rank,mid_pop_ptr->ind[i].tag,
			mid_pop_ptr->ind[i].eval,mid_pop_ptr->ind[i].gen);
		}
		printf("\n In generation loop:After crossover matepop, the matepop are:");
		printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
		for (i=0;i<popsize;i++) {
			printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			mate_pop_ptr->ind[i].chrom1[0],mate_pop_ptr->ind[i].chrom1[1],
			mate_pop_ptr->ind[i].chrom2[0],mate_pop_ptr->ind[i].chrom2[1],
			mate_pop_ptr->ind[i].fit[0],mate_pop_ptr->ind[i].fit[1],
			mate_pop_ptr->ind[i].cons[0],mate_pop_ptr->ind[i].cons[1],
			mate_pop_ptr->ind[i].overallcons,mate_pop_ptr->ind[i].cub_len,
			mate_pop_ptr->ind[i].rank,mate_pop_ptr->ind[i].tag,
			mate_pop_ptr->ind[i].eval,mate_pop_ptr->ind[i].gen);
		}
		printf("\n In generation loop:After crossover matepop, the newpop are:");
		printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
		for (i=0;i<popsize;i++) {
			printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			new_pop_ptr->ind[i].chrom1[0],new_pop_ptr->ind[i].chrom1[1],
			new_pop_ptr->ind[i].chrom2[0],new_pop_ptr->ind[i].chrom2[1],
			new_pop_ptr->ind[i].fit[0],new_pop_ptr->ind[i].fit[1],
			new_pop_ptr->ind[i].cons[0],new_pop_ptr->ind[i].cons[1],
			new_pop_ptr->ind[i].overallcons,new_pop_ptr->ind[i].cub_len,
			new_pop_ptr->ind[i].rank,new_pop_ptr->ind[i].tag,
			new_pop_ptr->ind[i].eval,new_pop_ptr->ind[i].gen);
		}*/
		
	 /***********************************Mutation************************************/
	 //adjacent mutation for chrom1 and polynomial mutation for chrom2
		//debug to check if tag assigning is working in mutation.c
		/*for(i=0;i<popsize;i++) {
			new_pop_ptr->ind[i].tag=0;
		}*/
		mutation(new_pop_ptr);
		
		//debug
		/*printf("\n===================Mutation===========================\n");
		printf("\n In generation loop:After mutating newpop, the initial oldpop are:");
		printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
		for (i=0;i<popsize;i++) {
			printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			old_pop_ptr->ind[i].chrom1[0],old_pop_ptr->ind[i].chrom1[1],
			old_pop_ptr->ind[i].chrom2[0],old_pop_ptr->ind[i].chrom2[1],
			old_pop_ptr->ind[i].fit[0],old_pop_ptr->ind[i].fit[1],
			old_pop_ptr->ind[i].cons[0],old_pop_ptr->ind[i].cons[1],
			old_pop_ptr->ind[i].overallcons,old_pop_ptr->ind[i].cub_len,
			old_pop_ptr->ind[i].rank,old_pop_ptr->ind[i].tag,
			old_pop_ptr->ind[i].eval,old_pop_ptr->ind[i].gen);
		}
		printf("\n In generation loop:After mutating newpop, the midpop are:");
		printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
		for (i=0;i<popsize;i++) {
			printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			mid_pop_ptr->ind[i].chrom1[0],mid_pop_ptr->ind[i].chrom1[1],
			mid_pop_ptr->ind[i].chrom2[0],mid_pop_ptr->ind[i].chrom2[1],
			mid_pop_ptr->ind[i].fit[0],mid_pop_ptr->ind[i].fit[1],
			mid_pop_ptr->ind[i].cons[0],mid_pop_ptr->ind[i].cons[1],
			mid_pop_ptr->ind[i].overallcons,mid_pop_ptr->ind[i].cub_len,
			mid_pop_ptr->ind[i].rank,mid_pop_ptr->ind[i].tag,
			mid_pop_ptr->ind[i].eval,mid_pop_ptr->ind[i].gen);
		}
		printf("\n In generation loop:After mutating newpop, the matepop are:");
		printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
		for (i=0;i<popsize;i++) {
			printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			mate_pop_ptr->ind[i].chrom1[0],mate_pop_ptr->ind[i].chrom1[1],
			mate_pop_ptr->ind[i].chrom2[0],mate_pop_ptr->ind[i].chrom2[1],
			mate_pop_ptr->ind[i].fit[0],mate_pop_ptr->ind[i].fit[1],
			mate_pop_ptr->ind[i].cons[0],mate_pop_ptr->ind[i].cons[1],
			mate_pop_ptr->ind[i].overallcons,mate_pop_ptr->ind[i].cub_len,
			mate_pop_ptr->ind[i].rank,mate_pop_ptr->ind[i].tag,
			mate_pop_ptr->ind[i].eval,mate_pop_ptr->ind[i].gen);
		}
		printf("\n In generation loop:After mutating newpop, the newpop are:");
		printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
		for (i=0;i<popsize;i++) {
			printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			new_pop_ptr->ind[i].chrom1[0],new_pop_ptr->ind[i].chrom1[1],
			new_pop_ptr->ind[i].chrom2[0],new_pop_ptr->ind[i].chrom2[1],
			new_pop_ptr->ind[i].fit[0],new_pop_ptr->ind[i].fit[1],
			new_pop_ptr->ind[i].cons[0],new_pop_ptr->ind[i].cons[1],
			new_pop_ptr->ind[i].overallcons,new_pop_ptr->ind[i].cub_len,
			new_pop_ptr->ind[i].rank,new_pop_ptr->ind[i].tag,
			new_pop_ptr->ind[i].eval,new_pop_ptr->ind[i].gen);
		}*/
		
		
	 /*************Decode the integer chromosome of the child population*************/ 
		if (no_intevar>0) 
			decode(new_pop_ptr);
		//debug
		/*printf("\n The decoded newpop after mutation are:");
		printf("\n DV1 DV2   D_DV1x   D_DV1y   D_DV1z   D_DV1u   D_DV1v   D_DV1x   D_DV1y   D_DV1z   D_DV1u   D_DV1v");
		for (i=0;i<popsize;i++) {
			printf("\n %d   %d  %.4f    %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f ",
			new_pop_ptr->ind[i].chrom1[0],new_pop_ptr->ind[i].chrom1[1],
			new_pop_ptr->ind[i].xchrom1[0],new_pop_ptr->ind[i].ychrom1[0],
			new_pop_ptr->ind[i].zchrom1[0],new_pop_ptr->ind[i].uchrom1[0],
			new_pop_ptr->ind[i].vchrom1[0],new_pop_ptr->ind[i].xchrom1[1],
			new_pop_ptr->ind[i].ychrom1[1],new_pop_ptr->ind[i].zchrom1[1],
			new_pop_ptr->ind[i].uchrom1[1],new_pop_ptr->ind[i].vchrom1[1]);
		}
		printf("\n In generation loop:After decoding mutated newpop, the newpop are:");
		printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
		for (i=0;i<popsize;i++) {
			printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			new_pop_ptr->ind[i].chrom1[0],new_pop_ptr->ind[i].chrom1[1],
			new_pop_ptr->ind[i].chrom2[0],new_pop_ptr->ind[i].chrom2[1],
			new_pop_ptr->ind[i].fit[0],new_pop_ptr->ind[i].fit[1],
			new_pop_ptr->ind[i].cons[0],new_pop_ptr->ind[i].cons[1],
			new_pop_ptr->ind[i].overallcons,new_pop_ptr->ind[i].cub_len,
			new_pop_ptr->ind[i].rank,new_pop_ptr->ind[i].tag,
			new_pop_ptr->ind[i].eval,new_pop_ptr->ind[i].gen);
		}*/
		
	 /*Evaluation the objective and constraint functions of the child population*/
		t1=clock();
		evaluatepop(new_pop_ptr, g);
		t2=clock();
		//debug
		/*printf("\n================Evaluatepop in generation loop=================\n");
		printf("\n In generation loop:After evaluatepop  newpop, the newpop are:");
		printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
		for (i=0;i<popsize;i++) {
			printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			new_pop_ptr->ind[i].chrom1[0],new_pop_ptr->ind[i].chrom1[1],
			new_pop_ptr->ind[i].chrom2[0],new_pop_ptr->ind[i].chrom2[1],
			new_pop_ptr->ind[i].fit[0],new_pop_ptr->ind[i].fit[1],
			new_pop_ptr->ind[i].cons[0],new_pop_ptr->ind[i].cons[1],
			new_pop_ptr->ind[i].overallcons,new_pop_ptr->ind[i].cub_len,
			new_pop_ptr->ind[i].rank,new_pop_ptr->ind[i].tag,
			new_pop_ptr->ind[i].eval,new_pop_ptr->ind[i].gen);
		}*/
		
		/*if (no_cons==0) {
			new_pop_ptr->ind_ptr->overallcons=0.0;
			ranking(new_pop_ptr);
		}
		else {
			rankcon(new_pop_ptr);
		}
		printf("\n In generation loop:After ranking/rankcon newpop, the newpop are:");
		printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
		for (i=0;i<popsize;i++) {
			printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			new_pop_ptr->ind[i].chrom1[0],new_pop_ptr->ind[i].chrom1[1],
			new_pop_ptr->ind[i].chrom2[0],new_pop_ptr->ind[i].chrom2[1],
			new_pop_ptr->ind[i].fit[0],new_pop_ptr->ind[i].fit[1],
			new_pop_ptr->ind[i].cons[0],new_pop_ptr->ind[i].cons[1],
			new_pop_ptr->ind[i].overallcons,new_pop_ptr->ind[i].cub_len,
			new_pop_ptr->ind[i].rank,new_pop_ptr->ind[i].tag,
			new_pop_ptr->ind[i].eval,new_pop_ptr->ind[i].gen);
		}*/
	  
	 /************************Selection keeping fronts alive************************/
		 //Elitism and sharing implementation
		t3=clock();
		keepalive(old_pop_ptr, new_pop_ptr, mate_pop_ptr, g+1);
		t4=clock();
		//debug
		/*printf("\n\n=====================Keepalive==========================\n");
		printf("\n In generation loop:After keepalive, the initial oldpop are:");
		printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
		for (i=0;i<popsize;i++) {
			printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			old_pop_ptr->ind[i].chrom1[0],old_pop_ptr->ind[i].chrom1[1],
			old_pop_ptr->ind[i].chrom2[0],old_pop_ptr->ind[i].chrom2[1],
			old_pop_ptr->ind[i].fit[0],old_pop_ptr->ind[i].fit[1],
			old_pop_ptr->ind[i].cons[0],old_pop_ptr->ind[i].cons[1],
			old_pop_ptr->ind[i].overallcons,old_pop_ptr->ind[i].cub_len,
			old_pop_ptr->ind[i].rank,old_pop_ptr->ind[i].tag,
			old_pop_ptr->ind[i].eval,old_pop_ptr->ind[i].gen);
		}
		printf("\n In generation loop:After keepalive, the newpop are:");
		printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
		for (i=0;i<popsize;i++) {
			printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			new_pop_ptr->ind[i].chrom1[0],new_pop_ptr->ind[i].chrom1[1],
			new_pop_ptr->ind[i].chrom2[0],new_pop_ptr->ind[i].chrom2[1],
			new_pop_ptr->ind[i].fit[0],new_pop_ptr->ind[i].fit[1],
			new_pop_ptr->ind[i].cons[0],new_pop_ptr->ind[i].cons[1],
			new_pop_ptr->ind[i].overallcons,new_pop_ptr->ind[i].cub_len,
			new_pop_ptr->ind[i].rank,new_pop_ptr->ind[i].tag,
			new_pop_ptr->ind[i].eval,new_pop_ptr->ind[i].gen);
		}
		printf("\n In generation loop:After keepalive, the newly formed matepop are:");
		printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
		for (i=0;i<popsize;i++) {
			printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			mate_pop_ptr->ind[i].chrom1[0],mate_pop_ptr->ind[i].chrom1[1],
			mate_pop_ptr->ind[i].chrom2[0],mate_pop_ptr->ind[i].chrom2[1],
			mate_pop_ptr->ind[i].fit[0],mate_pop_ptr->ind[i].fit[1],
			mate_pop_ptr->ind[i].cons[0],mate_pop_ptr->ind[i].cons[1],
			mate_pop_ptr->ind[i].overallcons,mate_pop_ptr->ind[i].cub_len,
			mate_pop_ptr->ind[i].rank,mate_pop_ptr->ind[i].tag,
			mate_pop_ptr->ind[i].eval,mate_pop_ptr->ind[i].gen);
		}*/
		
		/*if (no_intevar>0)
			decode(mate_pop_ptr);*/
		//debug
		/*printf("\n The decoded matepop after keepalive are:");
		printf("\n DV1 DV2   D_DV1x   D_DV1y   D_DV1z   D_DV1u   D_DV1v   D_DV1x   D_DV1y   D_DV1z   D_DV1u   D_DV1v");
		for (i=0;i<popsize;i++) {
			printf("\n %d   %d  %.4f    %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.4f ",
			mate_pop_ptr->ind[i].chrom1[0],mate_pop_ptr->ind[i].chrom1[1],
			mate_pop_ptr->ind[i].xchrom1[0],mate_pop_ptr->ind[i].ychrom1[0],
			mate_pop_ptr->ind[i].zchrom1[0],mate_pop_ptr->ind[i].uchrom1[0],
			mate_pop_ptr->ind[i].vchrom1[0],mate_pop_ptr->ind[i].xchrom1[1],
			mate_pop_ptr->ind[i].ychrom1[1],mate_pop_ptr->ind[i].zchrom1[1],
			mate_pop_ptr->ind[i].uchrom1[1],mate_pop_ptr->ind[i].vchrom1[1]);
		}*/
		
	 /*********************************Report printing*****************************/
		
		//debug
		/*printf("\n In generation loop:Before selectbest, the old bestpop are:");
		printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
		for (i=0;i<popsize;i++) {
			printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			best_pop_ptr->ind[i].chrom1[0],best_pop_ptr->ind[i].chrom1[1],
			best_pop_ptr->ind[i].chrom2[0],best_pop_ptr->ind[i].chrom2[1],
			best_pop_ptr->ind[i].fit[0],best_pop_ptr->ind[i].fit[1],
			best_pop_ptr->ind[i].cons[0],best_pop_ptr->ind[i].cons[1],
			best_pop_ptr->ind[i].overallcons,best_pop_ptr->ind[i].cub_len,
			best_pop_ptr->ind[i].rank,best_pop_ptr->ind[i].tag,
			best_pop_ptr->ind[i].eval,best_pop_ptr->ind[i].gen);
		}*/
		//select best pop
		if (g==0) { //for the first generation form the first bestpop from matepop
			copypop(mate_pop_ptr, best_pop_ptr);
		}
		else {
			selectbest(mate_pop_ptr, best_pop_ptr, g);
		}
		//debug
		/*printf("\n In generation loop:After selectbest, the matepop are:");
		printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
		for (i=0;i<popsize;i++) {
			printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			mate_pop_ptr->ind[i].chrom1[0],mate_pop_ptr->ind[i].chrom1[1],
			mate_pop_ptr->ind[i].chrom2[0],mate_pop_ptr->ind[i].chrom2[1],
			mate_pop_ptr->ind[i].fit[0],mate_pop_ptr->ind[i].fit[1],
			mate_pop_ptr->ind[i].cons[0],mate_pop_ptr->ind[i].cons[1],
			mate_pop_ptr->ind[i].overallcons,mate_pop_ptr->ind[i].cub_len,
			mate_pop_ptr->ind[i].rank,mate_pop_ptr->ind[i].tag,
			mate_pop_ptr->ind[i].eval,mate_pop_ptr->ind[i].gen);
		}
		printf("\n In generation loop:After selectbest, the bestpop are:");
		printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen");
		for (i=0;i<popsize;i++) {
			printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d",
			best_pop_ptr->ind[i].chrom1[0],best_pop_ptr->ind[i].chrom1[1],
			best_pop_ptr->ind[i].chrom2[0],best_pop_ptr->ind[i].chrom2[1],
			best_pop_ptr->ind[i].fit[0],best_pop_ptr->ind[i].fit[1],
			best_pop_ptr->ind[i].cons[0],best_pop_ptr->ind[i].cons[1],
			best_pop_ptr->ind[i].overallcons,best_pop_ptr->ind[i].cub_len,
			best_pop_ptr->ind[i].rank,best_pop_ptr->ind[i].tag,
			best_pop_ptr->ind[i].eval,best_pop_ptr->ind[i].gen);
		}*/
		
		mogareport(g, best_pop_ptr, mate_pop_ptr, lastit);
		//                                 output.out,
	 /***************************Copy the newpop to oldpop**************************/
	 /***Assign oldpop the value of matepop(cross, mutate, keepalive_selection, 
		so the original oldpop is updated***/
		copypop(mate_pop_ptr, old_pop_ptr);
		for (i=0;i<popsize;i++) {
			old_pop_ptr->ind[i].tag=0;
		}
	  
	 /***************Print the fitness record for the last generation***************/
	   /* if (g==no_gener-1) { //for the last generation
			old_pop_ptr=&(matepop);
			for (f=0;f<popsize;f++) { //printing loop
				old_pop_ptr->ind_ptr=&(old_pop_ptr->ind[f]);
			    //for all feasible and non-dominating solutions
				if ((old_pop_ptr->ind_ptr->overallcons<=0.0) && 
					(old_pop_ptr->ind_ptr->rank==1)) { //if solution is feasible
					for (l=0;l<no_obj;l++) //for each objective
						fprintf(end_ptr, "%f\t", old_pop_ptr->ind_ptr->fit[l]);
					for (l=0;l<no_cons;l++) //for each constraint
						fprintf(end_ptr, "%f\t", old_pop_ptr->ind_ptr->cons[l]);
					if (no_cons>0)
						fprintf(end_ptr, "%f\t", old_pop_ptr->ind_ptr->overallcons);
					fprintf(end_ptr, "\n");
					if (no_realvar>0) { //for real variables
						for (l=0;l<no_realvar;l++)
							fprintf(g_var, "%f\t", old_pop_ptr->ind_ptr->chrom2[l]);
						fprintf(g_var, " ");
					} //loop over realvar ends
					if (no_intevar>0) { //for integer variables
						for (l=0;l<no_intevar;l++)
							fprintf(g_var, "%d\t", old_pop_ptr->ind_ptr->chrom1[l]);
					} //loop over integer var ends
					fprintf(g_var,"\n");
				} //feasibility check ends
			} //loop over f ends (printing)
	    } //end of the last generation */
		end=clock();
		//printf("\n The time to do evaluatepop of this generation is %f seconds", ((double) t2-t1)/CLOCKS_PER_SEC);
		//printf("\n The time to doNS of this generation is %f seconds", ((double) t4-t3)/CLOCKS_PER_SEC);
		//printf("\n The time to do this generation is %f seconds", ((double) end-start)/CLOCKS_PER_SEC);
		//system(" rm -f en* ");
		printf("\n The %d th Generation has finished. ", g+1);
    } /*end of the ith generation*/
  
  /****************************GNENRATION LOOP FINISHES*********************************/
  /*************************************************************************************/
  fprintf(rep_ptr, "No. OF CROSSOVER = %d\n", no_cross);
  fprintf(rep_ptr, "NO. OF MUTATION FOR INTEGER CHROMOSOME= %d\n", no_mut1);
  fprintf(rep_ptr, "NO. OF MUTATION FOR REAL-CODED CHROMOSOME= %d\n", no_mut2);
  fprintf(rep_ptr,
  "-----------------------------------------------------------------------\n");
  fprintf(rep_ptr,
  "----------------Now you can look in the ourput files-------------------\n");
  
  //close files
  fclose(rep_ptr);
  fclose(lastit);
  return(0);
}

