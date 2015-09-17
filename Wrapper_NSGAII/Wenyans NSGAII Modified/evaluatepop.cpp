/*this file evaluate the objective function values, the constraint 
violations and overall constraint violation of the newly formed inds in newpop*/



#include "moga.h"
#include <stdlib.h>
#include <vector>
#include "decision_variables_2_ascii.h"

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
extern int noeval; //no. of evaluations

/*function prototype*/
double sum(double *, int);
double min(double *, int);
double minf(float *, int);
double max(double *, int);
void errors(int value);



void evaluatepop(population *pop_ptr, int gen)
{
    

    
    
    
    /*===========================start population evaluation==================================*/
    pop_ptr->ind_ptr=&(pop_ptr->ind[0]);
    pop_ptr->maxrank=0;
    
    for (int i=0;i<popsize;i++)
    { //this loop asigns fit, cons & overalcons to each ind
        pop_ptr->ind_ptr=&(pop_ptr->ind[i]);
        if (pop_ptr->ind_ptr->tag==1)          ///////// JEFF: What is this. Perhaps a tag to
            ///////// specify whether the obj functs and
            ///////// constraints have been evaluated yet?
        {
            
            /****************************************************************************
             *                                                                           *
             *          Jeff's interpretation of variables in this function:             *
             *                                                                           *
             *   1. It looks like Wenyan's NSGAII can consider decision variable vectors *
             *       consisting of both float and integer types.                         *
             *   2. To make matters more confusing, it looks like the integer decision   *
             *       variables are stored in pop_ptr->ind_ptr->chrom1. Then, in a funct  *
             *       decode, which is called before evaluatepop, these are variables are *
             *       given a physical meaning through mapping the int. values to things  *
             *       such as pipe size etc, which are stored in the arrays               *
             *          xchrom1[maxchrom1],                                              *
             *          ychrom1[maxchrom1],                                              *
             *          zchrom1[maxchrom1],                                              *
             *          uchrom1[maxchrom1],                                              *
             *          vchrom1[maxchrom1];                                              *
             *       These were then copied into local array variables defined in this   *
             *       function, in Wenyans original code.                                 *
             *       NOTE: Copying takes time and memory. There is no reason not to use  *
             *       the variable values as stored in the pop struct which is passed     *
             *           to this function as a pointer. [As far as I can tell            *
             *                                                                           *
             *    MY CHANGES                                                             *
             *       a. Integers dvs copied into the vector x,                           *
             *                                     copied from pop_ptr->ind_ptr->chrom1  *
             *                                                                           *
             *       b. Floats dvs copied into the vector r,                             *
             *                                     copied from pop_ptr->ind_ptr->chrom2  *
             *                                                                           *
             *                                                                           *
             *                                                                           *
             *****************************************************************************/
            
            
            std::vector<int> x(pop_ptr->ind_ptr->chrom1, pop_ptr->ind_ptr->chrom1 + no_intevar);
            std::vector<double> r(pop_ptr->ind_ptr->chrom2, pop_ptr->ind_ptr->chrom2 + no_realvar);
            
            /****************************************************************************
             *                                                                           *
             *          Jeff's interpretation of what Michael needs                      *
             *    1. It looks as though Michael is optimising a problem with only        *
             *        integer decision variables, and we print these to an ascii file    *
             *    2. The name of this ascii file has not been specified. We will call it *
             *          "dv.txt"                                                         *
             ****************************************************************************/
            
            DecisionVariable2Ascii<int> myWriter("dv.txt");
            myWriter(x);
            
            
            /****************************************************************************
             *                                                                           *
             *          Evaluate the objectives                                          *
             *              1. Call the system excutable                                 *
             *              2. Read in the obj. funct. values and constraints            *
             ****************************************************************************/
            
            //Call the system executable
            /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             *                                                                           *
             *      REPLACE "A.OUT" WITH THE COMMAND LINE CALL, TO EXECUTE YOUR CODE     *
             *                                                                           *
             !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
            system("a.out");
            
            //Read in the obj funct. values and constraints
            bool is_obj = no_obj > 0 ? true : false;
            bool is_cnstr = no_cons > 0 ? true : false;
            Ascii2ObjectiveValues<double, double> myReader("testObj.txt", is_obj, is_cnstr);
            std::shared_ptr<ObjectivesAndConstraints<double, double> > obj_and_cnstr = myReader();
            
            
            /****************************************************************************
             *                                                                           *
             *          Jeff's interpretation of outputs required from this function:    *
             *                                                                           *
             *   1. The objective functions get put into the fit[] array of each ind in  *
             *           as accessed through the pop_ptr pointer. MAKE SURE TO SPECIFY   *
             *           NUMBER OF OBJECTIVE FUNCTIONS, CONSTRAINTS ETC.                 *
             *       a. It looks to me that inputs for the GA are specified in moga.c or *
             *          are asked of the user from the command line when run             *
             *                                                                           *
             ****************************************************************************/
            
            // Sanity checks
            if (obj_and_cnstr->objectives.size() != no_obj)
            {
                std::cerr << "ERROR: Wrong number of objectives\n";
            }
            
            if (obj_and_cnstr->constraints.size() != no_cons)
            {
                std::cerr << "ERROR: Wrong number of constraints\n";
            }
            
            
            
            /****************************************************************************
             *                                                                           *
             *          NOTE ON OVERALLCONS                                              *
             *                                                                           *
             *   Wenyans code uses a variable overallcons (associated with an individual *
             *   to specify the overall constraint.  Mike - your input file defined this *
             *   as a vector of threshold values for constrained MUSIC outputs -         *
             *           constraint met if = 0.0,                                        *
             *           constraint violated if > 0.0                                    *
             *                                                                           *
             *      WHAT I'VE DONE.                                                      *
             *   If the "cons" value is > then the "overallcons" value, then I work out  *
             *   how much it is over (cons - overallcons).                               *
             *      from Wenyan's code, I think overallcons should be positive. With     *
             *      greater values indicating greater constraint violation               *
             ****************************************************************************/
            
            
            for (int k=0;k<no_obj;k++)
                pop_ptr->ind_ptr->fit[k]=obj_and_cnstr->objectives[k]; /*assign fitness to fitness array fit[no_obj] of ith indi */
            if (no_cons>0)
            {
                double error = 0.0;
                for (int k=0;k<no_cons;k++)
                {
                    double constraint_violation = obj_and_cnstr->constraints[k]
                                                    - obj_and_cnstr->constraint_thresholds[k];
                    pop_ptr->ind_ptr->cons[k] = (constraint_violation > 0) ?
                                                    constraint_violation : 0;
                    error += pop_ptr->ind_ptr->cons[k];
                }
                pop_ptr->ind_ptr->overallcons=error; /*assign the overall constarints to the individual*/
            }
            
            noeval+=1;
            pop_ptr->ind_ptr->eval=noeval;
            pop_ptr->ind_ptr->gen=gen+1;
        } //if the tag==1
    } //check each ind in the population
    
}


