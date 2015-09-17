/*this file taks inputs from input files or the screen for the GA*/

#include "moga.h"

extern int no_gener, /*number of generations*/
    no_obj,  /*no. of objective functions*/
	no_cons, /*no. of constraints*/
	no_pro, /*no. of property*/
	//tot_option, //no. of total options
    no_option[maxchrom1],  /*list of no. of options for each integer decision variable*/
    no_intevar, /*no. of integer decision variables = size of chromosome1*/
	no_realvar, /*no. of real decision variables = size of chromosome2*/
	no_mut1, /*no. of mutations happened to chrom1*/
	no_mut2, /*no. of mutations happened to chrom2*/
	no_cross, /*no. of crossovers happened*/
	no_cross_real, /*no. of real variables crossed*/
	ans; //if limits on real variables are rigid, ans=1
	
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

extern population oldpop, *old_pop_ptr,
		newpop, *new_pop_ptr,
	    matepop,  *mate_pop_ptr,
	    midpop, *mid_pop_ptr; //will be used in select.c

void input(FILE *rep_ptr)
{
	int i, j, k; //counters
	int defaulto, no1;
	int nogroup=11, nd[11], no[11], DVcount;
	//int nd1, nd2, nd3, nd4, nd5, nd6, nd7, nd8, nd9, nd10, nd11;
	//int no1, no2, no3, no4, no5, no6, no7, no8, no9, no10, no11;
	double xtemp[maxno_option],
		ytemp[maxno_option],
		ztemp[maxno_option],
		utemp[maxno_option],
		vtemp[maxno_option];
	
	printf("---------------------------------------------------------------------\n\n\n");
	printf("This is a Multi-objective Genetic Algorithm with constraints handling\n\n\n");
	printf("---------------------------------------------------------------------\n\n\n");
	printf("Give the values of problem specified variables\n\n\n");
	printf("---------------------------------------------------------------------\n");
	printf("---------------------------------------------------------------------\n");
	
	printf("\nGive number of generations: \n");
	scanf("%d", &no_gener);
	printf("\nGive population size: \n");
	scanf("%d", &popsize);
	if (popsize>maxpop) {
		printf("\nIncerase the maximum population size!\n");
		exit(1);
	}
	printf("\nGive number of objective functions: \n");
	scanf("%d",&no_obj);
	if (no_obj>maxobj) {
		printf("\nIncerase the maximum number of objective functions!\n");
		exit(1);
	}
	printf("\nGive number of constraints: \n");
	scanf("%d",&no_cons);
	if (no_cons>maxcons) {
		printf("\nIncrease the maximum number of constraints!\n");
		exit(1);
	}
	printf("\nGive number of properties: \n");
	scanf("%d",&no_pro);
	if (no_pro>maxpro) {
		printf("\nIncrease the maximum number of properties!\n");
		exit(1);
	}
	printf("\nGive the probability of crossover (0 to 1): \n");
	scanf("%lf", &pc);
	printf("\nGive the probability of mutation (0 to 1): \n");
	scanf("%lf", &pm);
	/*printf("\nGive the initial random seed (0 to 1): ");
	scanf("%f", &seed);*/
	//=============================for integer variables===========================//
	printf("\nGive the number of integer decision variables: \n");
	scanf("%d",&no_intevar);
	if (no_intevar>maxchrom1) {
		printf("\nIncrease the number of maximum integer variables: maxchrom1!\n");
		exit(1);
	}
	//input no. of options for each Integer DV
	if (no_intevar>0) {
		printf("\nAre the input options the same for all integer decision variables(1 for yes;0 for no)?\n");
		scanf("%d", &defaulto);
		if (defaulto==1) {
			printf("\nThe number of options for integer decision variables are: \n");
			scanf("%d", &no1);
			for (i=0;i<no_intevar;i++)
				no_option[i]=no1;
			printf("\nNow input the lookout array of 5 properties of each option.\n");
			for (j=0;j<no1;j++) {
				printf("\nInput 5 properties(x,y,z,u,v) of %d th option.\n",j+1);
				scanf("%lf %lf %lf %lf %lf", &xtemp[j], &ytemp[j], &ztemp[j], &utemp[j], &vtemp[j]);
			}
			for (i=0;i<no_intevar;i++) {
				for (j=0;j<no_option[i];j++) {
					xoption[i][j]=xtemp[j];
					yoption[i][j]=ytemp[j];
					zoption[i][j]=ztemp[j];
					uoption[i][j]=utemp[j];
					voption[i][j]=vtemp[j];
				}
			}
		}
		else { //when the options for different Integer DV are different
			/*printf("\nGive number of options for each integer variable: ");
			for (i=0;i<no_intevar;i++) {
				printf("\nGive number of options of the %d th DV: ", i+1);
				scanf("%d", &no_option[i]);
				if (no_option[i]>maxno_option) {
					printf("Increase the maximum number of options!");
					exit(1);
				}
			}
			printf("\nNow input the look out array of 3 properties of each option.\n");
			for (i=0;i<no_intevar;i++) { //for each integer DV
				for (j=0;j<no_option[i];j++) {
					printf("\nInput 5 properties(x,y,z,u,v) of %d th option, %d th integer DV: \n",
					j+1, i+1);
					scanf("%lf %lf %lf %lf %lf", &xoption[i][j], &yoption[i][j], &zoption[i][j],&uoption[i][j],&voption[i][j]);
				}
			}*/
			
			/*========================WDSA 2012 BNWII inputs below=================================*/
			printf("\nGive number of groups of integer DV: \n");
			scanf("%d", &nogroup);
			DVcount=0;
			for (i=0;i<nogroup;i++) {
				printf("\nGive number of integer DV in group %d: \n", i+1);
				scanf("%d", &nd[i]);
				printf("\nGive number of options for integer DV in group %d: \n", i+1);
				scanf("%d", &no[i]);
				if (no[i]>maxno_option) {
					printf("\nIncrease the maximum number of options!\n");
					exit(1);
				}
				else {
					for (k=DVcount;k<DVcount+nd[i];k++) { //assign number of options to each DV
						no_option[k]=no[i];
					}
				}
				for (j=0;j<no[i];j++) { //input the no[i] options for DV in group i+1
					printf("\nInput 5 properties(x,y,z,u,v) of option %d for group %d: \n",j+1, i+1);
					scanf("%lf %lf %lf %lf %lf", &xtemp[j], &ytemp[j], &ztemp[j],&utemp[j],&vtemp[j]);
				}
				for (k=DVcount;k<DVcount+nd[i];k++) { //assign the options to each DV in group i+1 
					for (j=0;j<no[i];j++) { //input the 25 ptions
						xoption[k][j]=xtemp[j];
						yoption[k][j]=ytemp[j];
						zoption[k][j]=ztemp[j];
						uoption[k][j]=utemp[j];
						voption[k][j]=vtemp[j];
					}
				}
				DVcount+=nd[i];	
			}
			
			
			
			
			/*//GROUP 1
			printf("\nGive number of DVs in the 1st group: ");
			scanf("%d", &nd1);
			printf("\nGive number of options of the 1st group: ");
			scanf("%d", &no1);
			if (no2>maxno_option) {
				printf("Increase the maximum number of options!");
				exit(1);
			}
			for (j=0;j<no1;j++) { //input the 25 options for DV in group 1
				printf("\nInput 5 properties(x,y,z,u,v) of %d th option for group 1: \n",j+1);
				scanf("%lf %lf %lf %lf %lf", &xtemp[j], &ytemp[j], &ztemp[j],&utemp[j],&vtemp[j]);
			}
			for (i=0;i<nd1;i++) { //asign the options to each DV in group 1 
				for (j=0;j<no1;j++) { //input the 25 ptions
					xoption[i][j]=xtemp[j];
					yoption[i][j]=ytemp[j];
					zoption[i][j]=ztemp[j];
					uoption[i][j]=utemp[j];
					voption[i][j]=vtemp[j];
				}
			}
			//GROUP 2
			printf("\nGive number of DVs in the 2nd group: ");
			scanf("%d", &nd2);
			printf("\nGive number of options of the 2nd group: ");
			scanf("%d", &no2);
			if (no2>maxno_option) {
				printf("Increase the maximum number of options!");
				exit(1);
			}
			for (j=0;j<no2;j++) { //input the 24 options for DV in group 2
				printf("\nInput 5 properties(x,y,z,u,v) of %d th option for group 2: \n",j+1);
				scanf("%lf %lf %lf %lf %lf", &xtemp[j], &ytemp[j], &ztemp[j],&utemp[j],&vtemp[j]);
			}
			for (i=nd1;i<nd1+nd2;i++) { //asign the options to each DV in group 1 
				for (j=0;j<no2;j++) { //input the 25 ptions
					xoption[i][j]=xtemp[j];
					yoption[i][j]=ytemp[j];
					zoption[i][j]=ztemp[j];
					uoption[i][j]=utemp[j];
					voption[i][j]=vtemp[j];
				}
			}*/
			
			
			/*========================WDSA 2012 BNWII inputs above=================================*/
		}
		//debug
		/*printf("\nNow print out the lookout table:\n");
		printf("\n----------------------------------------------\n");
		printf("Decision variable   xoption   yoption   zoption  uoption  voption");
		for (i=0;i<no_intevar;i++) {
			for (j=0;j<no_option[i];j++) {
				printf("\nDV index %d, Option%d:\t\t\t %g\t%g\t%g\t%g\t%g\t\n",i,j+1,
				xoption[i][j],yoption[i][j],zoption[i][j],uoption[i][j],voption[i][j]);
			}
		}
		printf("\n----------------------------------------------\n");
		exit(1);*/
	}
	//=================================================================================
	//==============================for real value variables=========================//
	printf("\nGive number of real value decision variables: \n");
	scanf("%d", &no_realvar);
	if (no_realvar>maxchrom2) {
		printf("\nIncrease the maximum number of real value variables!\n");
		exit(1);
	}
	if (no_realvar>0) {
		for (i=0;i<no_realvar;i++) { 
			printf("\nGive lower & upper limits of the %d th real variable: \n", i+1);
			scanf("%lf %lf", &reallim[i][0], &reallim[i][1]);
		}
		printf("\nIf limits of real value variables are rigid (1 for yes): \n");
		scanf("%d", &ans);
		printf("\nGive the distribution index for SBX (5 to 20): \n");
		scanf("%lf", &di);
		printf("\nGive the distribution index for polynomial mutation (5 to 50): \n");
		scanf("%lf", &dim);
	}
	//=============================================================================//
	//=============================================================================//
	
	//===========Print GA parameters and problem parameters in output.dat==========//
	fprintf(rep_ptr, "GA Parameters: \n");
	fprintf(rep_ptr, "---------------------------------------------------------------\n");
	fprintf(rep_ptr, "Population Size ---> %d\n", popsize);
	fprintf(rep_ptr, "Number of Generations ---> %d\n", no_gener);
	fprintf(rep_ptr, "Number of Objective Functions ---> %d\n", no_obj);
	fprintf(rep_ptr, "Number of constraints ---> %d\n", no_cons);
	fprintf(rep_ptr, "Number of properties ---> %d\n", no_pro);
	fprintf(rep_ptr, "Probability of crossover ---> %f\n", pc);
	fprintf(rep_ptr, "Probability of mutation ---> %f\n", pm);
	fprintf(rep_ptr, "Random seed ---. %f\n", seed);
	if (no_intevar>0) {
		fprintf(rep_ptr, "Number of integer variables (pipes)---> %d\n", no_intevar);
	}
	if (no_realvar>0) {
		fprintf(rep_ptr, "Number of real variables ---> %d\n", no_realvar);
		for (i=0;i<no_realvar;i++) {
			fprintf(rep_ptr, "\nLower and Upper limits for %d th real variable are : %f %f",
			i+1, reallim[i][0], reallim[i][1]);
		}
		if (ans==1)
			fprintf(rep_ptr, "\nReal variable bounds are rigid");
		else
			fprintf(rep_ptr, "\nReal variable bounds are not rigid");
	}
	return;
}
