/*this file copy the first individual to the second individual, 
which makes a exact copy of the first individual in the second one*/

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

individual copyind(individual ind1)
{
	int i;
	individual ind2;
	
	ind2.flag=ind1.flag;
	ind2.rank=ind1.rank;
	ind2.overallcons=ind1.overallcons;
	ind2.cub_len=ind1.cub_len;
	ind2.tag=ind1.tag;
	ind2.eval=ind1.eval;
	ind2.gen=ind1.gen;
	
	for (i=0;i<no_intevar;i++) {
		ind2.chrom1[i]=ind1.chrom1[i];
		ind2.xchrom1[i]=ind1.xchrom1[i];
		ind2.ychrom1[i]=ind1.ychrom1[i];
		ind2.zchrom1[i]=ind1.zchrom1[i];
		ind2.uchrom1[i]=ind1.uchrom1[i];
		ind2.vchrom1[i]=ind1.vchrom1[i];
	}
	for (i=0;i<no_realvar;i++)
		ind2.chrom2[i]=ind1.chrom2[i];
	for (i=0;i<no_obj;i++)
		ind2.fit[i]=ind1.fit[i];
	for(i=0;i<no_cons;i++)
		ind2.cons[i]=ind1.cons[i];
	for(i=0;i<no_pro;i++)
		ind2.property[i]=ind1.property[i];
	return ind2;
}
