/*this file contains the mutation operator, which applies adjcent mutation 
to the integer string and polynomial mutation to the real-coded string. In the adjcent 
mutation, once the mutation is decided to take place, it has 50% chance to change 
the value to the next upper value and 50% chance to the next lower value; if the
value is the boundary value, it will be mutated the value next to it*/

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
	uoption[maxchrom1][maxno_option],
	voption[maxchrom1][maxno_option],
	reallim[maxchrom2][2]; /*the lower and upper limts for real avriables*/

extern int popsize; /*population size*/

extern population oldpop, *old_pop_ptr,
		newpop, *new_pop_ptr,
	    matepop,  *mate_pop_ptr,
	    midpop, *mid_pop_ptr; //will be used in select.c

/*function prototype*/
double randomperc(void); 

void mutation(population *pop_ptr)
{
	int //r, //random number for random mutation
	    i, //individual counter in the population
	    j, //bit counter in chrom1
		k, //bit counter in chrom2
		x; //the value of bits in chrom1
	double rnd,
	    rnd1,
		rnd2,
		yl, //the lower bound of ith realvar y
		yu, //the upper bound of ith realvar y
		y, //the value of bits in chrom2
		u, val,
		expp,
		delta,
		deltaq;
	
    for (i=0;i<popsize;i++) { //the loop over the population
		//printf("\n-----------------------------------\n");
		//printf("\nFor the %dth individual:", i+1);
		for (j=0;j<no_intevar;j++) {//mutation loop for chrom1
			x=pop_ptr->ind[i].chrom1[j]; //x is jth bit in chrom1 of ith ind
			rnd=randomperc();
			//printf("\nFor the %dth Integer DV: rnd= %f; pm= %f ", j+1,rnd,pm);
			//printf("\nBefore mutation the values is: %d", x);
			if (rnd<=pm) {//perform mutation to jth integer DV in ith ind
				no_mut1++;
				rnd1=randomperc(); 
				x=floor(rnd1*no_option[j])+1; //random mutation
				/*if (x==no_option[j]) {//if value of x is the last option of jth bit
					x=no_option[j]-1;
					//if (rnd1<=0.5)
						//x=no_option[j]-1;//assign x the value of second last option
					//else
						//x=no_option[j];
				}
				else if (x==1) {//if x selected the first option of jth bit
					x=2;
					//if (rnd1<=0.5)
						//x=2;//assign x the value of second variable
					//else
						//x=1;
				}
				else { //for the inside values
					if (rnd1<=0.5) 
						x=x-1;//mutate x to the next upper option
					else
						x=x+1;//mutate x to the next lower option
				}*/  //adjcent mutation finished
				if (x==0) {
					printf("ERROR:Mutated value of %dth ind %dth Integer DV is 0!", i+1,j+1);
					exit(-1);
				}
			}
			if (x!=pop_ptr->ind[i].chrom1[j]) {
				pop_ptr->ind[i].tag=1; //set tag for evaluation
			}
			//printf("\nAfter mutation the values is: %d\n", x);
			pop_ptr->ind[i].chrom1[j]=x;//assign mutated value back to the array
		} //integer mutation for chrom1 ends
		for (k=0;k<no_realvar;k++) {//mutation loop over chrom2
			y=pop_ptr->ind[i].chrom2[k];
			rnd=randomperc();
			//printf("\nFor the %dth R-DV: rnd= %f; pm= %f ", k+1,rnd,pm);
			//printf("\nBefore mutation the values is: %f", y);
			if (rnd<=pm) {//perform mutation to kth real bit in ind i
				no_mut2++;
				yl=reallim[k][0];
				yu=reallim[k][1];
				if (y>yl) { //another posibility is y=yl
					//calculate delta
					if ((y-yl)<=(yu-y)) //y is closer to yl
						delta=(y-yl)/(yu-y);
					else //y is closer to yu
						delta=(yu-y)/(y-yl);
					expp=1.0/(dim+1.0);
					
					//calculate deltaq
					rnd2=randomperc();
					if (rnd2<=0.5) {
						u=1.0-delta;
						val=2.0*rnd2+(1-2.0*rnd2)*(pow(u,(dim+1)));
						deltaq=pow(val,expp)-1.0;
					}
					else {//if rnd2>0.5
						u=1.0-delta;
						val=2.0*(1.0-rnd2)+2.0*(rnd2-0.5)*(pow(u,(dim+1)));
						deltaq=1.0-pow(val,expp);
					}
					
					y=y+deltaq*(yu-yl);
					if (y<yl) 
						y=yl;
					if (y>yu) 
						y=yu;
				}
				else {//if y==yl
					u=randomperc();
					y=u*(yu-yl)+yl;
				}
				if (y<yl || y>yu) {
					printf("\nERROR:Mutated value of %dth ind %dth R-DV out of range!",
					i+1,k+1);
					exit(-2);
				}
				//printf("\nAfter mutation the values is: %f", y);
			} //mutation loop for kth bit ends
			else { //rnd>pm
				y=y;
			}
			if (y!=pop_ptr->ind[i].chrom2[k]) {
				pop_ptr->ind[i].tag=1;
			}
			pop_ptr->ind[i].chrom2[k]=y;
		} //mutation loop over chrom2 ends
	} //mutation loop over population ends
	return;
}
