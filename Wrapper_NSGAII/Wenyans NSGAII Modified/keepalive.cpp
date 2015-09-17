/*This file adds the first population (oldpop) and the second population (newpop)
into a single population pool, which is the globalpop. Then the individuals in this
global pool are sorted firstly non-dominatedly into different ranks; then according
to their crowding distances (cub_len). The top half individuals will be selected into 
the third population, which is the matepop. Then the individuals in this mating 
pool will be used for selection, crossover and mutation; thus non-dominated sorting is
applied and diversity (density estimation: cub_len) and elitism (oldpop+newpop)are kept.*/

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

globpop my_globpop, *my_global_pop_ptr;

//variables used in this file
int lastrank;
double fpara1[2*maxpop][2];

//function prototy
void grank(int gen); //ranking the globalpop when there is no constraints
int indcmp(double *ptr1, double *ptr2); //compare two inds
void grankcon(int gen); //ranking the globalpop when there are constraints
void gshare(int rnk); //Calculating Crowding distance cub_len
void gsort(int rnk, int sel); //Sort the arrays of cub_len
void sort(int m1); //sort array fpara1[i][1] in ascending order of the fitness
int rnd (int low, int high); //Fetch a single random integer between low and high including the bounds


void keepalive(population *pop_ptr1, population *pop_ptr2, population *pop_ptr3, int gen)
{ //the 3 populations are oldpop, newpop(through crossover & mutation) and matepop
	int i,j,k,l,ii,jj, //counters
		m, //maxrank of globalpop
		st, //no. of inds considered by now
		sel, //no. of inds that haven't been considered by now; sel=popsize-st
		pool; //no. of accumulated inds in rank 1 to rank i+1
		
	int *chrom1_ptr1, //pointer of chrom1 of globalpop
	    *chrom1_ptr2; //pointer of chrom1 of pop3
		
	double *chrom2_ptr1, //pointer of chrom2 of globalpop
		*chrom2_ptr2; //pointer of chrom2 of pop3
	
	
	/*****************************Forming the Global Pool****************************/
	for (i=0;i<popsize;i++) {
		if (no_intevar>0) { //copy chrom1 into the global pool
			for (k=0;k<no_intevar;k++) { 
				my_globpop.chrom1[i][k]=pop_ptr1->ind[i].chrom1[k];
				my_globpop.chrom1[i+popsize][k]=pop_ptr2->ind[i].chrom1[k];
			}
		} //integer variables are copied
		if (no_realvar>0) { //copy chrom2 into the global pool
			for (k=0;k<no_realvar;k++) {
				my_globpop.chrom2[i][k]=pop_ptr1->ind[i].chrom2[k];
				my_globpop.chrom2[i+popsize][k]=pop_ptr2->ind[i].chrom2[k];
			}
		}// real value variables are copied
		for (l=0;l<no_obj;l++) { //copy fitnesses into the global pool
			my_globpop.fit[i][l]=pop_ptr1->ind[i].fit[l];
			my_globpop.fit[i+popsize][l]=pop_ptr2->ind[i].fit[l];
		} //fitnesses are copied
		for (ii=0;ii<no_cons;ii++) { //copy constaints into the global pool
			my_globpop.cons[i][ii]=pop_ptr1->ind[i].cons[ii];
			my_globpop.cons[i+popsize][ii]=pop_ptr2->ind[i].cons[ii];
		} //constraints are copied
		for (ii=0;ii<no_pro;ii++) { //copy properties into the global pool
			my_globpop.property[i][ii]=pop_ptr1->ind[i].property[ii];
			my_globpop.property[i+popsize][ii]=pop_ptr2->ind[i].property[ii];
		} //properties are copied
		//copy overall constraints into the global pool
		my_globpop.overallcons[i]=pop_ptr1->ind[i].overallcons;
		my_globpop.overallcons[i+popsize]=pop_ptr2->ind[i].overallcons;
		my_globpop.tag[i]=pop_ptr1->ind[i].tag;
		my_globpop.tag[i+popsize]=pop_ptr2->ind[i].tag;
		my_globpop.eval[i]=pop_ptr1->ind[i].eval;
		my_globpop.eval[i+popsize]=pop_ptr2->ind[i].eval;
		my_globpop.gen[i]=pop_ptr1->ind[i].gen;
		my_globpop.gen[i+popsize]=pop_ptr2->ind[i].gen;
		//initialize the crowding distance to zero
		my_globpop.cub_len[i]=0.0;
		my_globpop.cub_len[i+popsize]=0.0;
	}
	my_global_pop_ptr=&(my_globpop);
	//debug
	/*printf("\nPrint out globalpop right after it is formed from oldpop+newpop");
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
	/*************************Global Pool is Formed******************************/
	
	//Finding global ranks
	//printf("\n***********grank/grankcon starts here in keepalive************");
	if (no_cons==0)
		grank(gen);
	else
		grankcon(gen);
	//printf("\n***********grank/grankcon finishes here in keepalive***********");
	m=my_globpop.maxrank;
	//debug
	/*printf("\nPrint out globalpop after grankcon.c");
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
	
	//Calculate Crowding distance cub_len by sharing the fitnesses
	//printf("\nThe maxrank m is : %d", m);
	for (i=0;i<m;i++) {
		gshare(i+1);
	}
	//Initialize the falgs of the population to zero
	for (i=0;i<2*popsize;i++)
		my_globpop.flag[i]=0;
	//debug
	/*printf("\nPrint out globalpop after gshare.c and initializing flag");
	printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen  flag");
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
	
	//decide which solutions belong to pop3 -the global-matepop
	st=0; //no. of inds considered by now
	pool=0; //no. of inds in rank 1 to rank i+1
	//Elitism applied here--rank globalpop to form global-matepop
	//printf("\nglobalpop.rankno[0]= %d", globalpop.rankno[0]);
	if(my_globpop.rankno[0]>=popsize) { //if the first rank no.>popsize
		gsort(1,popsize);
	}
	else {
		for (i=0;i<m;i++) { // m is the global maxrank
		pool+=my_globpop.rankno[i];//no. of inds in rank 1 to rank i+1
			if(pool>=popsize) {
				st=pool-my_globpop.rankno[i]; //no. of inds considered by the last rank
				for (j=0;j<i;j++) { //for the ranks before current rank
					for (k=0;k<2*popsize;k++) {
						if(my_globpop.rank[k]==j+1)
							my_globpop.flag[k]=1;
					}
				}
				sel=popsize-st; //no. of ind in curernt rank need to sort according cub_len
				lastrank=i+1;
				pop_ptr3->rankno[i]=sel;
				gsort(i+1, sel);
				/*printf("\nlastrank= %d", lastrank);
				printf("\nPool= %d", pool);
				printf("\nglobalpop.rankno[%d]= %d",i,globalpop.rankno[i]);
				printf("\nst=pool-globalpop.rankno[i]=%d", st);
				printf("\nNo. of ind in last rank can be selected are sel=%d", sel);*/
				goto next1;
			}
		} // now go to next rank i+1
	} //loop over i ends, all inds in globalpop are considered
	next1:
	//debug
	/*printf("\nPrint out globalpop after assigning flag");
	printf("\n DV1 DV2   RDV1   RDV2   fit1   fit2   cons1  cons2  Ocons  cub_len  rank tag eval gen  flag");
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
	
	/*****select inds from globalpop to form pop3(matepop)*****/
	k=0; //counter of pop3 (matepop)
	for (i=0,k=0;i<2*popsize && k<popsize;i++) {
		if (my_globpop.flag[i]==1) { //when i is dominating
			if (no_intevar>0) { //for integer variables
				chrom1_ptr1=&(my_globpop.chrom1[i][0]);
				pop_ptr3->ind_ptr=&(pop_ptr3->ind[k]);
				chrom1_ptr2=&(pop_ptr3->ind_ptr->chrom1[0]);
				for (j=0;j<no_intevar;j++) //assign values of integer vars of pop3
					*chrom1_ptr2++=*chrom1_ptr1++;
				}//loop over inter variables ends
			if (no_realvar>0) { //for real variables
				chrom2_ptr1=&(my_globpop.chrom2[i][0]);
				pop_ptr3->ind_ptr=&(pop_ptr3->ind[k]);
				chrom2_ptr2=&(pop_ptr3->ind_ptr->chrom2[0]);
				for (j=0;j<no_realvar;j++)
					*chrom2_ptr2++=*chrom2_ptr1++;
			}//loop over real variables ends
			//copy other values
			for (j=0;j<no_obj;j++)
				pop_ptr3->ind[k].fit[j]=my_globpop.fit[i][j];
			pop_ptr3->ind[k].cub_len=my_globpop.cub_len[i];
			if (no_cons!=0) 
				pop_ptr3->ind[k].overallcons=my_globpop.overallcons[i];
			for (jj=0;jj<no_cons;jj++)
				pop_ptr3->ind[k].cons[jj]=my_globpop.cons[i][jj];
			for (jj=0;jj<no_pro;jj++)
				pop_ptr3->ind[k].property[jj]=my_globpop.property[i][jj];
			pop_ptr3->ind[k].rank=my_globpop.rank[i];
			pop_ptr3->ind[k].tag=my_globpop.tag[i];
			pop_ptr3->ind[k].eval=my_globpop.eval[i];
			pop_ptr3->ind[k].gen=my_globpop.gen[i];
			k++; //increment of pop3 counter
		} //loop over i dominating ends
		
		//debug
		//printf("\n In keepalive form matepop i=%d, k=%d ", i,k);
	} //pops(matepop) is formed (loop over i ends)
	pop_ptr3->maxrank=lastrank;
	//debug
	//assign cub_len from globalpop to oldpop
	/*for (i=0;i<popsize;i++)
		pop_ptr1->ind[i].cub_len=globalpop.cub_len[i];
	printf("\nprint out the cub_len and rank of both oldpop and just formed matepop\n");
	printf("\nIND    fit,cub_len, rank in oldpop    fit, cub_len, rank in matepop");
	for (i=0;i<popsize;i++)
		printf("\n#%d  \t%g\t%g\t%g\t%d \t\t\t%g\t%g\t%g \t%d",i+1,
		pop_ptr1->ind[i].fit[0],pop_ptr1->ind[i].fit[1],pop_ptr1->ind[i].cub_len,pop_ptr1->ind[i].rank,
		pop_ptr3->ind[i].fit[0],pop_ptr3->ind[i].fit[1],pop_ptr3->ind[i].cub_len,pop_ptr3->ind[i].rank);
	*/
	return;
}

/*****************ranking the globalpop when there is no constraints*************/
void grank(int gen)
{
	int i,j,k, /*counters*/
		rnk, /*rank*/
		val, /*value obtained after comparing two individuals*/
		nondom, /*no. of non-domonated members*/
		gpopsize, //popsize of globalpop
		gflag[2*maxpop], //array of flags in globalpop
		q;//the  number of individuals at a rank
	double *ptr1,*ptr2;
	//FILE *gr;
	
	/*gr=fopen("g_rank_record.out", "a");
	fprintf(gr, "Generation no. = %d\n", gen);*/
	
	/****************Ranking starts here**********************/
	//initialize ranks to zero
	rnk=0;
	nondom=0;
	gpopsize=2*popsize;
	//initialize all flags to 2, which means dominated
	for (i=0;i<gpopsize;i++)
		gflag[i]=2;
	for (k=0;k<gpopsize;k++) { //loop for each ind in the globalpop
		q=0; //initialize the no. of individual at the rank to zero
		for (j=0;j<gpopsize;j++) {
			if (gflag[j]!=1) 
				break; 
		}
		if (j==gpopsize)
			break; //break, when all inds are assigned ranks
		rnk=rnk+1;
		for (j=0;j<gpopsize;j++) { //set the flags of dominated inds to 2
			if (gflag[j]==0)
				gflag[j]=2;
		} 
		//select two individuals to compare in order to assign ranks
		for (i=0;i<gpopsize;i++) { //select the first ind i
			if (gflag[i]!=0 && gflag[i]!=1) {
				ptr1=&(my_global_pop_ptr->fit[i][0]);
				for (j=0;j<gpopsize;j++) { //select the second ind j
					if (i!=j) {
						if (gflag[j]!=1) { //if j is not dominating
							ptr2=&(my_global_pop_ptr->fit[j][0]);
							val=indcmp(ptr1,ptr2); //compare i and j by fit
							switch (val) {
								case 2: //j is dominating
									gflag[i]=0; //i is dominated
									break;
								case 1: //i is dominating
									gflag[j]=0;
									break;
								case 3: //i and j are nondominated
									nondom++;
									if (gflag[j]!=0)
										gflag[j]=3;
									break;
								default:
									printf("\n!grank has problems in comparing inds!");
									break;
							}
							//printf("\n---------------------------------------------\n");
							//printf("\nfor the comparation in grank.c, val= %d", val);
							//printf("\nflag of first ind(%d): %d", i+1, gflag[i]);
							//printf("\nflag of second ind(%d): %d", j+1, gflag[j]);
						}
					}
				} //loop over j ends
				if (j==gpopsize) {
					my_global_pop_ptr->rank[i]=rnk;
					if (gflag[i]!=0) {
						gflag[i]=1; //i is dominating
						my_global_pop_ptr->rankarray[rnk-1][q]=i;
						q++;
					}
				}
			}
		} //loop over i ends
		my_global_pop_ptr->rankno[rnk-1]=q; //there are q inds at rank rnk
	} //loop over globalpop ends
	my_global_pop_ptr->maxrank=rnk;
	
	/*fprintf(gr, "       RANK      No. of Individuals\n");
	for (i=0;i<rnk;i++)
		fprintf(gr, "\t%d\t%d\n", i+1, globalpop.rankno[i]);
	fclose(gr);*/
	//debug
	/*printf("\n====================In grank.c=======================\n");
	for (i=0;i<gpopsize;i++) {
		printf("\nThe rank for ind %d is: %d", i+1, globalpop.rank[i]);
		printf("\nThe flag for ind %d is: %d", i+1, gflag[i]);
	}
	printf("\nMaxrank of the population is: %d", globalpop.maxrank);
	for (i=0;i<global_pop_ptr->maxrank;i++) {
		printf("\nthe rankno at rank %d is : %d", i+1, globalpop.rankno[i]);
	}
	printf("\n========================================================\n");*/
	return;
}
/*=============================================================================*/

/******************ranking the globalpop when there are constraints***************/
void grankcon(int gen)
{
	int i,j,k, /*counters*/
		rnk, /*rank*/
		val, /*value obtained after comparing two individuals*/
		nondom, /*no. of non-domonated members*/
		gpopsize, //popsize of globalpop
		gflag[2*maxpop], //array of flags in globalpop
		q;//the  number of individuals at a rank
	double *ptr1,*ptr2, //fitness pointers
		*err_ptr1, *err_ptr2; //overall constraints pointers
	//double min_fit; //the min fitness of a better rank
      //delta_fit; //min difference of the crowdign distance between two ranks
	//FILE *gr;
	
	//gr=fopen("g_rank_record.out", "a");
	//fprintf(gr, "Generation no. = %d\n", gen);
	
	/****************Ranking starts here**********************/
	
	//initialize ranks to zero
	rnk = 0;
	nondom = 0;
    gpopsize = 2*popsize;
	
	//initialize all flags to 2
	for (i=0;i<gpopsize;i++)
		gflag[i]=2;
	for (k=0;k<gpopsize;k++) { //strat ranking each ind in the globalpop
		q=0; //initialize the no. of ind in the rank to zero
		for (j=0;j<gpopsize;j++) {
			if (gflag[j]!=1) {
				//printf("\nbreak1, in grankcon");
				break; //break, when all inds are assigned ranks
			}
		}
		if (j==gpopsize) {
			//printf("\nbreak2, in grankcon");
			break;
		}
		for (j=0;j<gpopsize;j++) { //set the flags of dominated inds to 2
			if (gflag[j]==0)
				gflag[j]=2;
		}
		rnk=rnk+1;
		//select two individuals to compare in order to assign ranks
		for (i=0;i<gpopsize;i++) { //select the first ind i
			if (gflag[i]!=1 && gflag[i]!=0) { //flag check
				ptr1=&(my_global_pop_ptr->fit[i][0]);
				err_ptr1=&(my_global_pop_ptr->overallcons[i]);
				for (j=0;j<gpopsize;j++) { //select the second ind j
					if (i!=j) {
						if (gflag[j]!=1) { //j is not dominating
							ptr2=&(my_global_pop_ptr->fit[j][0]);
							err_ptr2=&(my_global_pop_ptr->overallcons[j]);
							if (*err_ptr1>0.0 && *err_ptr2>0.0) {//both infeasible
								if (*err_ptr1<*err_ptr2) {//i is better
									gflag[j]=0; //j is dominated
								}
								else {
									if (*err_ptr1>*err_ptr2) {//j is better
										gflag[i]=0;
										//printf("\nbreak3, in grankcon");
										break;
									}
									else {//*err_ptr1==*err_ptr2
										nondom++;
										if (gflag[j]!=0)
											gflag[j]=3;
									}
								}
							}
							else {//not both infeasible 
								if (*err_ptr1>0.0 && *err_ptr2<=0.0) {
									gflag[i]=0; //i is infeasible & j is feasible
									//printf("\nbreak4, in grankcon");
									break;
								}
								else {
									if (*err_ptr1<=0.0 && *err_ptr2>0.0) {
										gflag[j]=0;//i is feasible & j is infeasible
									}
									else { //both feasible
										val=indcmp(ptr1,ptr2);
										switch (val) {
										case 2:
											gflag[i]=0; //i is dominated
											//printf("\nbreak5, in grankcon");
											break;
										case 1:
											gflag[j]=0; //j is dominated
											//printf("\nbreak6, in grankcon");
											break;
										case 3:
											nondom++;
											if (gflag[j]!=0)
												gflag[j]=3; //i & j are nondominating
											//printf("\nbreak7, in grankcon");
											break;
										default:
											printf("!grankcon has problems in comparing inds!\n");
											break;
										}
									}//both feasible
								}
							}
						} //loop when j is not dominating ends
					} //loop over i!=j ends
				} //loop over j ends
				if (j==gpopsize) {
					my_global_pop_ptr->rank[i]=rnk;
					if (gflag[i]!=0) {
						gflag[i]=1; //i is dominating
						my_global_pop_ptr->rankarray[rnk-1][q]=i;
						q++;
					}
				}
			} //flag chenk ends
		} //loop over i ends
		my_global_pop_ptr->rankno[rnk-1]=q;
		//printf("\nIn grankcon:Glogalpop counter k= %d; The rank rnk = %d; q= %d", k, rnk, q);
	}  //loop over globalpop ends
	my_global_pop_ptr->maxrank=rnk; //due to break statement
	/*fprintf(gr,"    RANK     No of Individuals\n");
	for (i=0;i<rnk;i++)
		fprintf(gr, "\t%d\t%d\n", i+1,globalpop.rankno[i]);
	fclose(gr);*/
	//debug
	/*printf("\n====================In grankcon.c=======================\n");
	for (i=0;i<gpopsize;i++) {
		printf("\nThe rank for ind %d is: %d", i+1, globalpop.rank[i]);
		printf("\nThe flag for ind %d is: %d", i+1, gflag[i]);
	}
	printf("\nMaxrank of the population is: %d", globalpop.maxrank);
	for (i=0;i<globalpop.maxrank;i++) {
		printf("\nthe rankno at rank %d is : %d", i+1, globalpop.rankno[i]);
	}
	printf("\n========================================================\n");*/
	return;
}
/*=============================================================================*/

/******************Calculating Crowding distance cub_len at each rank******************/
void gshare(int rnk)
{
	double length[2*maxpop][no_obj+1], //
		max, //max fitness in the array fparal[i][1]
		min, //min fitness in the array fparal[i][1]
		diff; //difference between min and max
	int p,j, //counters
	    m1, //no. of inds at the rnkth rank 
		a; //the ind at the rnkth rank
	
	m1=my_globpop.rankno[rnk-1]; //no. of ind at rank rnk
	for (j=0;j<no_obj;j++) { //for each objective
		for (p=0;p<m1;p++) { //for each ind at the rnkth rank
			fpara1[p][0]=0.0; //initialize to zero
			fpara1[p][1]=0.0;
		}
		for (p=0;p<m1;p++) {
			a=my_globpop.rankarray[rnk-1][p];
			fpara1[p][0]=(double)a; //index of ind
			fpara1[p][1]=my_globpop.fit[a][j]; //corrispondinf fit of ind a
		}
		sort(m1); //sort the array fpara1[i][1] in ascending order of the fitness
		max=fpara1[m1-1][1];
		min=fpara1[0][1];
		diff=max-min;
		/*printf("\nfor the %dth rank,%dth obj, maxfit, minfit and diff are: %g %g %g",
		rnk,j+1,max,min,diff);*/
		if (diff<0.0) {
			printf("\nERROR:Something wrong in keepalive.gshare (max-min<0)\n");
			exit(-1);
		}
		for (p=0;p<m1;p++) {
			if (p==0 || p==(m1-1)) { //for boundary solutions
				length[p][0]=fpara1[p][0];
				length[p][j+1]=INF;
			}
		}
		for (p=1;p<m1-1;p++) {//for intermediate solutions
			if (diff==0.0) { //fit of all ind on the rank the same
				length[p][0]=fpara1[p][0];
				length[p][j+1]=0.0;
			}
			else {
				length[p][0]=fpara1[p][0];
				length[p][j+1]=fabs(fpara1[p+1][1]-fpara1[p-1][1])/diff;
				//normalize the crowding distance
			}
		}
	} //loop over each obj ends
	//debug print out the matrix length
	/*printf("\n===========> The %dth rank:", rnk);
	printf("\n(Gshare)The matrix 'length': 1st row is index+1(ind); others are distance for each obj");
	for (p=0;p<m1;p++) {
		printf("\n%d \t%.2e \t%.2e", (int)length[p][0]+1,length[p][1],length[p][2]);
	}*/
	//assign the length to globalpop
	for (p=0;p<m1;p++) {
		a=length[p][0];
		if (p==0 || p==m1-1) {
			my_globpop.cub_len[a]=INF;
		}
		else {
			for (j=0;j<no_obj;j++)
				my_globpop.cub_len[a]=+length[p][j+1];
		my_globpop.cub_len[a]=my_globpop.cub_len[a]/no_obj;
		}
	}
	
	return;
}
/*=============================================================================*/

/********************Sort the arrays of cub_len in descending order****************/
void gsort(int rnk, int sel)
{
	int i,j,a,q;
	double array[2*maxpop][2],temp,temp1;
	
	q=my_globpop.rankno[rnk-1];
	
	for (i=0;i<q;i++) { //for each ind at a rank
		array[i][0]=my_globpop.rankarray[rnk-1][i];
		a=my_globpop.rankarray[rnk-1][i];
		array[i][1]=my_globpop.cub_len[a];
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
		}
	}
	for (i=0;i<sel;i++) {
		a=array[i][0];
		my_globpop.flag[a]=1;
	}
	return;
}
/*=============================================================================*/

/**********Sort the array fpara1[i][1] in ascending order of the fitness**********/
void sort(int m1)
{
	double temp,temp1;
	int i1,k1;
	
	for (k1=0;k1<m1-1;k1++) {
		for (i1=k1+1;i1<m1;i1++) {
			if (fpara1[k1][1]>fpara1[i1][1]) {
				temp=fpara1[k1][1];
				temp1=fpara1[k1][0];
				fpara1[k1][1]=fpara1[i1][1];
				fpara1[k1][0]=fpara1[i1][0];
				fpara1[i1][1]=temp;
				fpara1[i1][0]=temp1;
			}
		}
	}
	return;
}
/*=============================================================================*/
