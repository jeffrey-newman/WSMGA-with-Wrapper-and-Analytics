/*this funcation keeps the top popsize rank one solution*/

#include "moga.h"

extern int no_gener, //number of generations
    no_obj,  //no. of objective functions
	no_cons, //no. of constraints
	no_pro, /*no. of property*/
    no_option[maxchrom1],  //list of no. of options for each integer decision variable
    no_intevar, //no. of integer decision variables = size of chromosome1
	no_realvar, //no. of real decision variables = size of chromosome2
	no_mut1, //no. of mutations happened to chrom1
	no_mut2, //no. of mutations happened to chrom2
	no_cross, //no. of crossovers happened
	no_cross_real, //no. of real variables crossed
	ans;
	
extern double seed, //random seed
	pc, //crossover probability for both integer and real chromosomes
	pm, //mutation probability for both interger and real chromosomes
	di, //distribution index for SBX
	dim, //distribution index for polynomial mutation
	xoption[maxchrom1][maxno_option], //lookout table of 1st property of no_intevar
	yoption[maxchrom1][maxno_option], //lookout table of 2nd property of no_intevar
	zoption[maxchrom1][maxno_option], //lookout table of t3rd property of no_intevar
	uoption[maxchrom1][maxno_option],
	voption[maxchrom1][maxno_option],
	reallim[maxchrom2][2]; //the lower and upper limts for real avriables
	
extern int popsize; //population size

//global structure defination
typedef struct {
	int maxrank, //max ranks of the global population
		rankarray[2*maxpop][2*maxpop], //array of ind's no. at a rank
        rankno[2*maxpop]; //no. of inds at a rank
	int chrom1[2*maxpop][maxchrom1], //integer variables
	    rank[2*maxpop], //rank of different inds
		tag[2*maxpop],
		eval[2*maxpop],
		gen[2*maxpop],
		flag[2*maxpop];
	double chrom2[2*maxpop][maxchrom2], //real value variables
		fit[2*maxpop][maxobj], //no_obj fitnesses of inds
		cub_len[2*maxpop], //crowding distance
		overallcons[2*maxpop], //overallcons of each ind 
		cons[2*maxpop][maxcons], // constraints of each ind in globalpop
		property[2*maxpop][maxpro]; // properties of each ind in globalpop
} globpop;

globpop globalpop, *global_pop_ptr;

void shuffle(population *pop_ptr);
void grank(int gen); //ranking the globalpop when there is no constraints
int indcmp(double *ptr1, double *ptr2); //compare two inds
void grankcon(int gen); //ranking the globalpop when there are constraints
void gshare(int rnk); //Calculating Crowding distance cub_len
void gbsort(int rnk, int sel); //Sort the arrays of cub_len
void sort(int m1); //sort array fpara1[i][1] in ascending order of the fitness
int rnd (int low, int high); //Fetch a single random integer between low and high including

void selectbest(population *pop_ptr1, population *pop_ptr2, int t)
{//        mate_pop_ptr, best_pop_ptr, generation number
	int i,j,k,m,jj;
	int pool,st,sel,lastrank;
	int *chrom1_ptr1, //pointer of chrom1 of globalpop
	    *chrom1_ptr2; //pointer of chrom1 of pop2 (bestpop)
	double *chrom2_ptr1, //pointer of chrom2 of globalpop
		*chrom2_ptr2; //pointer of chrom2 of pop2 (bestpop)
	
	//shuffle(pop_ptr1);
	//Form globalpop using both bestpop and matepop
	for (i=0;i<popsize;i++) {
		if (no_intevar>0) { //copy chrom1 into the global pool
			for (j=0;j<no_intevar;j++) { 
				globalpop.chrom1[i][j]=pop_ptr2->ind[i].chrom1[j];
				globalpop.chrom1[i+popsize][j]=pop_ptr1->ind[i].chrom1[j];
			}
		} //integer variables are copied
		if (no_realvar>0) { //copy chrom2 into the global pool
			for (j=0;j<no_realvar;j++) {
				globalpop.chrom2[i][j]=pop_ptr2->ind[i].chrom2[j];
				globalpop.chrom2[i+popsize][j]=pop_ptr1->ind[i].chrom2[j];
			}
		}// real value variables are copied
		for (j=0;j<no_obj;j++) { //copy fitnesses into the global pool
			globalpop.fit[i][j]=pop_ptr2->ind[i].fit[j];
			globalpop.fit[i+popsize][j]=pop_ptr1->ind[i].fit[j];
		} //fitnesses are copied
		for (j=0;j<no_cons;j++) { //copy constaints into the global pool
			globalpop.cons[i][j]=pop_ptr2->ind[i].cons[j];
			globalpop.cons[i+popsize][j]=pop_ptr1->ind[i].cons[j];
		} //constraints are copied
		for (j=0;j<no_pro;j++) { //copy properties into the global pool
			globalpop.property[i][j]=pop_ptr2->ind[i].property[j];
			globalpop.property[i+popsize][j]=pop_ptr1->ind[i].property[j];
		} //properties are copied
		//copy overall constraints into the global pool
		globalpop.overallcons[i]=pop_ptr2->ind[i].overallcons;
		globalpop.overallcons[i+popsize]=pop_ptr1->ind[i].overallcons;
		globalpop.tag[i]=pop_ptr2->ind[i].tag;
		globalpop.tag[i+popsize]=pop_ptr1->ind[i].tag;
		globalpop.eval[i]=pop_ptr2->ind[i].eval;
		globalpop.eval[i+popsize]=pop_ptr1->ind[i].eval;
		globalpop.gen[i]=pop_ptr2->ind[i].gen;
		globalpop.gen[i+popsize]=pop_ptr1->ind[i].gen;
		//initialize the crowding distance to zero
		globalpop.cub_len[i]=0.0;
		globalpop.cub_len[i+popsize]=0.0;
		globalpop.rank[i]=0;
		globalpop.rank[i+popsize]=0;
		globalpop.flag[i]=0;
		globalpop.flag[i+popsize]=0;
	}
	global_pop_ptr=&(globalpop);
	//debug
	/*printf("\nIn selectbest: Print out globalpop after it is formed from bestpop+newpop");
	printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen flag");
	for (i=0;i<2*popsize;i++) {
		printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d   %d",
			global_pop_ptr->chrom1[i][0],global_pop_ptr->chrom1[i][1],
			global_pop_ptr->chrom2[i][0],global_pop_ptr->chrom2[i][1],
			global_pop_ptr->fit[i][0],global_pop_ptr->fit[i][1],
			global_pop_ptr->cons[i][0],global_pop_ptr->cons[i][1],
			global_pop_ptr->overallcons[i],global_pop_ptr->cub_len[i],
			global_pop_ptr->rank[i],global_pop_ptr->tag[i],
			global_pop_ptr->eval[i],global_pop_ptr->gen[i],global_pop_ptr->flag[i]);
	}*/
	//========================= globalpop is formed =================================	
	
	if (no_cons==0)
		grank(t);
	else
		grankcon(t);
	//printf("\n***********bgrank/bgrankcon finishes here in selectbast***********");
	m=globalpop.maxrank;
	for (i=0;i<m;i++)
		gshare(i+1);
	
	//printf("\n In selectbest, gloablpop.rankno[0]= %d", globalpop.rankno[0]);
	st=0; //no. of inds considered by now
	pool=0; //no. of inds in rank 1 to rank i+1
	if(globalpop.rankno[0]>=popsize) { //if the first rank no.>popsize
		gbsort(1,popsize);
	}
	else {
		for (i=0;i<m;i++) { // m is the global maxrank
		pool+=globalpop.rankno[i];//no. of inds in rank 1 to rank i+1
			if(pool>popsize) {
				st=pool-globalpop.rankno[i]; //no. of inds considered by the last rank
				for (j=0;j<i;j++) { //for the ranks before current rank
					for (k=0;k<2*popsize;k++) {
						if(globalpop.rank[k]==j+1)
							globalpop.flag[k]=1;
					}
				}
				sel=popsize-st; //no. of ind in curernt rank need to sort according cub_len
				lastrank=i+1;
				pop_ptr2->rankno[i]=sel;
				gbsort(i+1, sel);
				/*printf("\nlastrank= %d", lastrank);
				printf("\nPool= %d", pool);
				printf("\nglobalpop.rankno[%d]= %d",i,globalpop.rankno[i]);
				printf("\nst=pool-globalpop.rankno[i]=%d", st);
				printf("\nNo. of ind in last rank can be selected are sel=%d", sel);*/
				goto next2;
			}
		} // now go to next rank i+1
	} //loop over i ends, all inds in globalpop are considered
	next2:
	//debug
	/*printf("\nIn selectbest: Print out globalpop after rank, share, assign flag");
	printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen flag");
	for (i=0;i<2*popsize;i++) {
		printf("\n %d   %d  %.4f     %.4f   %.4f   %.4f   %.4f   %.4f   %.4f   %.2e     %d   %d  %d   %d   %d",
			global_pop_ptr->chrom1[i][0],global_pop_ptr->chrom1[i][1],
			global_pop_ptr->chrom2[i][0],global_pop_ptr->chrom2[i][1],
			global_pop_ptr->fit[i][0],global_pop_ptr->fit[i][1],
			global_pop_ptr->cons[i][0],global_pop_ptr->cons[i][1],
			global_pop_ptr->overallcons[i],global_pop_ptr->cub_len[i],
			global_pop_ptr->rank[i],global_pop_ptr->tag[i],
			global_pop_ptr->eval[i],global_pop_ptr->gen[i],global_pop_ptr->flag[i]);
	}*/
	/**********************select inds from globalpop to update pop2(bestpop)*****/
	k=0; //counter of pop2 (bestpop)
	for (i=0,k=0;i<2*popsize && k<popsize;i++) {
		if (globalpop.flag[i]==1) { //when i is dominating
			if (no_intevar>0) { //for integer variables
				chrom1_ptr1=&(globalpop.chrom1[i][0]);
				pop_ptr2->ind_ptr=&(pop_ptr2->ind[k]);
				chrom1_ptr2=&(pop_ptr2->ind_ptr->chrom1[0]);
				for (j=0;j<no_intevar;j++) //assign values of integer vars of pop3
					*chrom1_ptr2++=*chrom1_ptr1++;
				}//loop over inter variables ends
			if (no_realvar>0) { //for real variables
				chrom2_ptr1=&(globalpop.chrom2[i][0]);
				pop_ptr2->ind_ptr=&(pop_ptr2->ind[k]);
				chrom2_ptr2=&(pop_ptr2->ind_ptr->chrom2[0]);
				for (j=0;j<no_realvar;j++)
					*chrom2_ptr2++=*chrom2_ptr1++;
			}//loop over real variables ends
			//copy other values
			for (j=0;j<no_obj;j++)
				pop_ptr2->ind[k].fit[j]=globalpop.fit[i][j];
			pop_ptr2->ind[k].cub_len=globalpop.cub_len[i];
			if (no_cons!=0) 
				pop_ptr2->ind[k].overallcons=globalpop.overallcons[i];
			for (jj=0;jj<no_cons;jj++)
				pop_ptr2->ind[k].cons[jj]=globalpop.cons[i][jj];
			for (jj=0;jj<no_pro;jj++)
				pop_ptr2->ind[k].property[jj]=globalpop.property[i][jj];
			pop_ptr2->ind[k].rank=globalpop.rank[i];
			pop_ptr2->ind[k].tag=globalpop.tag[i];
			pop_ptr2->ind[k].eval=globalpop.eval[i];
			pop_ptr2->ind[k].gen=globalpop.gen[i];
			k++; //increment of pop3 counter
		} //loop over i dominating ends
		
		//debug
		//printf("\n In keepalive form matepop i=%d, k=%d ", i,k);
	} //pops(matepop) is formed (loop over i ends)
	pop_ptr2->maxrank=lastrank;
	
	return;
}

/********************Sort the arrays of cub_len in descending order****************/
void gbsort(int rnk, int sel)
{
	int i,j,a,q;
	double array[2*maxpop][3],temp,temp1;
	
	q=globalpop.rankno[rnk-1];
	
	for (i=0;i<q;i++) { //for each ind at a rank
		array[i][0]=globalpop.rankarray[rnk-1][i];
		a=globalpop.rankarray[rnk-1][i];
		array[i][1]=globalpop.cub_len[a];
		array[i][2]=globalpop.eval[a];
	}
	for (i=0;i<q;i++) {
		for (j=i+1;j<q;j++) {
			if (array[i][1]<array[j][1]) {
				temp=array[i][1];
				temp1=array[i][0];
				array[i][1]=array[j][1];
				array[i][0]=array[j][0];
				array[j][1]=temp;
				array[j][0]=temp1;
			}
			if (array[i][1]==0.0 && array[j][1]==0.0) {
				if(array[i][2]>array[j][2]) {
					temp=array[i][2];
					temp1=array[i][0];
					array[i][2]=array[j][2];
					array[i][0]=array[j][0];
					array[j][2]=temp;
					array[j][0]=temp1;
				}
			}
		}
	}
	for (i=0;i<sel;i++) {
		a=array[i][0];
		globalpop.flag[a]=1;
	}
	return;
}

