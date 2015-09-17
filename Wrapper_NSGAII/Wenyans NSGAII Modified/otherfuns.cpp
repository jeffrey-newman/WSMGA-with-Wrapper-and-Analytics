/*this file stores other small functions*/

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

/********************function that report errors related to EPANET*******************/
void errors(int value)
{
	/*if(value!=0 && value<=100) {//negative pressure warnings
		printf("\nERROR CODE: %i\n",value);
		ENclose();
	}*/
	if (value>100) {
		printf("\nERROR CODE: %i\n",value);
		//ENclose();
		//exit(8);//when there is an error, it will exit
	}
}  

/********************function that find sum of an array*******************/
double sum(double value[], int n)
{
	int i;
	double sum=0.0;
	
	for (i=0;i<n;i++) 
		sum+=value[i];
	return(sum);
}

/********************function that find minimum value in an array*******************/
double min(double value[], int n)
{
	int i;
	double min=value[0];
	
	for (i=0;i<n;i++) 
		if (min>value[i])
			min=value[i];
	return(min);
}

/********************function that find minimum value in an array of float numbers*******************/
double minf(float value[], int n)
{
	int i;
	double min=value[0];
	
	for (i=0;i<n;i++) 
		if (min>value[i])
			min=value[i];
	return(min);
}

/********************function that find maximum value in an array*******************/
double max(double value[], int n)
{
	int i;
	double max=value[0];
	
	for (i=0;i<n;i++) 
		if (max<value[i])
			max=value[i];
	return(max);
}

/*************function that compare two individuals using their fitness************/
int indcmp(double *ptr1, double *ptr2)
{
	double fit1[maxobj], fit2[maxobj];
	int   flag1=0, flag2=0;
	int i, value;
	
	for (i=0;i<no_obj;i++) {
		fit1[i]=*ptr1++;
		fit2[i]=*ptr2++;
	}
	for (i=0;i<no_obj;i++) {
		if (fit1[i]<fit2[i])
			flag1=1;
		else {
			if (fit1[i]>fit2[i])
				flag2=1;
		}
	}
	if (flag1==1 && flag2==0)
		value=1; //i is better
	else {
		if (flag1==0 && flag2==1)
			value=2; //j is better
		else
			value=3; //i and j are nondominated
	}
	//printf("\nthe value is : %d", value);
	return(value);
}
