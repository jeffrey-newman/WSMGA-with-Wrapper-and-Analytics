/*this file evaluate the objective function values, the constraint 
violations and overall constraint violation of the newly formed inds in newpop
Note:
1. The maximum number of objective and properties that can be used in this file 
are defined in the moga.h file
2. All constraints need to be constructed in the form of ConstraintFunction>=0*/

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

/*function prototype*/
double sum(double *, int);
double min(double *, int);
double max(double *, int);
void errors(int value);


void evaluatepop(population *pop_ptr, int gen)
{
	int i,j;
	int k;
	double error;
	double *x_ptr,  /*pointer to xchrom1*/
		*y_ptr, //pointer to ychrom1
		*z_ptr, //pointer to zchrom1
		*u_ptr,
		*v_ptr,
		x[maxchrom1], //array of xchrom1 values of one ind
		y[maxchrom1], //array of ychrom1 values of one ind
		z[maxchrom1], //array of zchrom1 values of one ind
		u[maxchrom1],
		v[maxchrom1],
		r[maxchrom2], //array of chrom2 values of one ind
	    *chrom2_ptr,  /*pointer to the array of real-coded variables*/
		//*fit_ptr,  //pointer to the array of objective function values
		f[maxobj],  /*array of fitness values for each individual*/
		//*err_ptr,  //pointer to the overall constraints violation
		cstr[maxcons], /*list of constraint functions*/
		pro[maxpro]; /*list of network properties*/
		
	 /*********************declare problem dependent variables********************/
	double x1, x2; 
	double obj1, obj2;

	
	

	//--------------------------Please do not change the section below	
	pop_ptr->ind_ptr=&(pop_ptr->ind[0]);
	pop_ptr->maxrank=0;
	
	for (i=0;i<popsize;i++) { //this loop asigns fit, cons & overalcons to each ind
		pop_ptr->ind_ptr=&(pop_ptr->ind[i]);
		if (pop_ptr->ind_ptr->tag==1) {
			x_ptr=&(pop_ptr->ind_ptr->xchrom1[0]);
			y_ptr=&(pop_ptr->ind_ptr->ychrom1[0]);
			z_ptr=&(pop_ptr->ind_ptr->zchrom1[0]);
			u_ptr=&(pop_ptr->ind_ptr->uchrom1[0]);
			v_ptr=&(pop_ptr->ind_ptr->vchrom1[0]);
			chrom2_ptr=&(pop_ptr->ind_ptr->chrom2[0]);
			if (no_intevar>0) {//when there are integer variables
				for (j=0;j<no_intevar;j++) {//for each integer DV in the ith ind
					x[j]=*x_ptr++; //array of xchrom1 values of ith ind
					y[j]=*y_ptr++; //array of ychrom1 values of ith ind
					z[j]=*z_ptr++; //array of zchrom1 values of ith ind
					u[j]=*u_ptr++;
					v[j]=*v_ptr++;
				}
			}
			if (no_realvar>0) {//when there are real variables
				for (j=0;j<no_realvar;j++)//for each real DV in the ith ind
					r[j]=*chrom2_ptr++; ///array of chrom2 values of ith ind
			}
	//----------------------------Please do not change the section above
			
			
			/*=============fitness and constraint functions are coded below============*/
			
			
			f[0]=obj1;
			f[1]=obj2;  
			cstr[0]=pow(x1,2)+pow(x2,2)-1-0.1*cos(16*atan(x1/x2));
			cstr[1]=-1*(pow((x1-0.5),2)+pow((x2-0.5),2));
			pro[0]=x1;	
			pro[1]=x2;  
			pro[2]=obj1;	
			pro[3]=obj2;	
			
			
		/*=====assign fitness and constraints to individuals in this generation====*/
			for (k=0;k<no_obj;k++)
				pop_ptr->ind_ptr->fit[k]=f[k]; /*assign fitness to fitness array fit[no_obj] of ith indi*/
			if (no_cons>0) {
				for (k=0;k<no_cons;k++)
					pop_ptr->ind_ptr->cons[k]=cstr[k]; 
				error=0.0;
				for (k=0;k<no_cons;k++) { /*calculate the overall constraints of ith individual*/
					if(cstr[k]<0.0) /*test if the constraint k is violated*/
						error-=cstr[k]; /*if it is violated then add it to error*/
				}
				pop_ptr->ind_ptr->overallcons=error; /*assign the overall constarints to the individual*/	
			}
			if (no_pro>0) { //copy network properties to current individual
				for (k=0;k<no_pro;k++)
					pop_ptr->ind_ptr->property[k]=pro[k];
			}
			noeval+=1;
			pop_ptr->ind_ptr->eval=noeval;
			pop_ptr->ind_ptr->gen=gen+1;
		} //if the tag==1
	} //check each ind in the population
	return;
}


