/*this file contains the crossover operator, which apply one point crossover to 
the integer chromosome chrom1 and SBX to the real-coded chromsome chrom2.*/

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
double randomperc(void);
void shuffle(population *pop_ptr);

void crossover(population *pop_ptr1,population *pop_ptr2)
{ //they are newpop(children) and matepop(parents)
	int i, 
	    y, //counter for parents in matepop
		n, //counter for children in newpop
		k, //counter for bits in chrom1
		j, //counter for bits in chrom2
	    mating_site;
    double rnd, rnd2;
	
	//variables for SBX
	double par1,par2,child1,child2, temp,
		y1, //the parent with smaller value
		y2, //the parent with larger value
		alpha, //
		beta, //spead factor for SBX
		betaq, //
		yl, //upper limit of jth variable in ith individual reallim[j][0]
		yu; //lower limit of jth variable in ith individual reallim[j][1]
	
	
	
	/**************shuffle the matepop(pop_ptr2) before applying crossover*************/
    shuffle(pop_ptr2);
	/******************************matepop is shuffled*********************************/
	
	rnd=randomperc();
	/*initialize population pointer*/
	pop_ptr1->ind_ptr=&(pop_ptr1->ind[0]);
	pop_ptr2->ind_ptr=&(pop_ptr2->ind[0]); 
	for (i=0,n=0,y=0;i<popsize/2;i++) {
		/*initialize individual and chromosome pointers*/
		//printf("\n--------------------------------------------\n");
		//printf("\ni=%d; n=%d; y=%d;", i,n,y);
		pop_ptr1->ind_ptr=&(pop_ptr1->ind[n]); //children in children pop
		pop_ptr2->ind_ptr=&(pop_ptr2->ind[y]); //parents in parents pop
		
		rnd=randomperc();
		//printf("\nrnd= %f; pc= %f ", rnd, pc);
		if (rnd<=pc) { //crossover will take place
			no_cross++;
			//printf("\nno_cross= %d ", no_cross);
			//apply one-point crossover to chrom1
			//printf("\nOne-Point Crossover For Integer DVs:");
			rnd=randomperc();
			mating_site=floor(rnd*no_intevar);
			if (mating_site==0)
				mating_site=1;
			//printf("\nmating_site= %d ", mating_site);
			for (k=0;k<no_intevar;k++) {//for each bit in chrom1
				if (k>(mating_site-1)) {//on the right hand side of mating_site
					//swap the bits of the parents to form childen
					pop_ptr1->ind[n].chrom1[k]=pop_ptr2->ind[y+1].chrom1[k];
					pop_ptr1->ind[n+1].chrom1[k]=pop_ptr2->ind[y].chrom1[k];
				}
				else {//on the left hand side of mating_site, chrom1 keeps the same
					pop_ptr1->ind[n].chrom1[k]=pop_ptr2->ind[y].chrom1[k];
					pop_ptr1->ind[n+1].chrom1[k]=pop_ptr2->ind[y+1].chrom1[k];
				}
				//change tag for evaluation
				pop_ptr1->ind[n].tag=1;
				pop_ptr1->ind[n+1].tag=1;
			}//one-point crossover to chrom1 ends
			
			//apply SBX to chrom2
			//printf("\n\nSBX For Real Value DVs:");
			for (j=0;j<no_realvar;j++) {//for each bit in chrom2
				//select two parents
				par1=pop_ptr2->ind[y].chrom2[j];
				par2=pop_ptr2->ind[y+1].chrom2[j];
				yl=reallim[j][0];//lower limt of variable j in chrom2 of ind i
				yu=reallim[j][1];//upper limt of variable j in chrom2 of ind i
				//choose a random number
				rnd=randomperc();
				//printf("\nrnd= %f ", rnd);
				if (rnd<=0.5) {//only half of the bits in chrom2 will be crossed
					no_cross_real++;
					if (fabs(par1-par2)>SIGMA) {//values of two parents are different
						if (par2>par1) {
							y2=par2; //y2 is the parent that has larger value
							y1=par1;
						}
						else {//if par2<par1
							y2=par1;
							y1=par2;
						}
						rnd2=randomperc();
						//for the first child
						beta=1.0+(2.0*(y1-yl)/(y2-y1));
						alpha=2.0-pow(beta,-(di+1.0));
						if (rnd2<=(1.0/alpha)) 
							betaq=pow((rnd2*alpha), (1.0/(di+1.0)));
						else
							betaq=pow((1.0/(2.0-rnd2*alpha)), (1.0/(di+1.0)));
						child1=0.5*((y1+y2)-betaq*(y2-y1));
						//for the second child
						beta=1.0+(2.0*(yu-y2)/(y2-y1));
						alpha=2.0-pow(beta, -(di+1.0));
						if (rnd2<=(1.0/alpha))
							betaq=pow((rnd2*alpha), (1.0/(di+1.0)));
						else
							betaq=pow((1.0/(2.0-rnd2*alpha)), (1.0/(di+1.0)));
						child2=0.5*((y1+y2)+betaq*(y2-y1));
						if (child1<yl) child1=yl;
						if (child1>yu) child1=yu;
						if (child2<yl) child2=yl;
						if (child2>yu) child2=y2;
						if (randomperc()<=0.5) {
							temp=child1;
							child1=child2;
							child2=temp;
						}
					}//loop over two different parents ends
					else {//two parents have similar values; fabs(par2-par1)<0.0
						child1=par1;
						child2=par2;
					}
				}//crossover to half bits in chrom2 ends
				else {//to the bits that will not perform crossover
					child1=par1;
					child2=par2;
				}
				
				//copy the values of children into the newpop
				pop_ptr1->ind[n].chrom2[j]=child1;
				pop_ptr1->ind[n+1].chrom2[j]=child2;
				//change tag for evaluation
				pop_ptr1->ind[n].tag=1;
				pop_ptr1->ind[n+1].tag=1;
			}//SBX to chrom2 ends (j loop ends)
			
		}//loop over crossover application ends
		else {//rnd>pc, crossover won't take place
			for (k=0;k<no_intevar;k++) {//for chrom1
				pop_ptr1->ind[n].chrom1[k]=pop_ptr2->ind[y].chrom1[k];
				pop_ptr1->ind[n+1].chrom1[k]=pop_ptr2->ind[y+1].chrom1[k];
			}
			for (j=0;j<no_realvar;j++) {//for chrom2
				pop_ptr1->ind[n].chrom2[j]=pop_ptr2->ind[y].chrom2[j];
				pop_ptr1->ind[n+1].chrom2[j]=pop_ptr2->ind[y+1].chrom2[j];
			}
			pop_ptr1->ind[n].tag=pop_ptr2->ind[y].tag;
			pop_ptr1->ind[n+1].tag=pop_ptr2->ind[y+1].tag;
			//if they are not crossed, copy the fit, cons and overallcons
			for (j=0;j<no_obj;j++) {
				pop_ptr1->ind[n].fit[j]=pop_ptr2->ind[y].fit[j];
				pop_ptr1->ind[n+1].fit[j]=pop_ptr2->ind[y+1].fit[j];
			}
			for (j=0;j<no_cons;j++) {
				pop_ptr1->ind[n].cons[j]=pop_ptr2->ind[y].cons[j];
				pop_ptr1->ind[n+1].cons[j]=pop_ptr2->ind[y+1].cons[j];
			}
			for (j=0;j<no_pro;j++) {
				pop_ptr1->ind[n].property[j]=pop_ptr2->ind[y].property[j];
				pop_ptr1->ind[n+1].property[j]=pop_ptr2->ind[y+1].property[j];
			}
			pop_ptr1->ind[n].overallcons=pop_ptr2->ind[y].overallcons;
			pop_ptr1->ind[n+1].overallcons=pop_ptr2->ind[y+1].overallcons;
		}//loop over noncrossover application ends
		pop_ptr1->ind[n].eval=pop_ptr2->ind[y].eval;
		pop_ptr1->ind[n+1].eval=pop_ptr2->ind[y+1].eval;
		pop_ptr1->ind[n].gen=pop_ptr2->ind[y].gen;
		pop_ptr1->ind[n+1].gen=pop_ptr2->ind[y+1].gen;
		n=n+2;
		y=y+2;
		//printf("\ncounter for children n= %d; counter for parents y= %d", n, y);
	}//loop over crossover ends
	return;
}

