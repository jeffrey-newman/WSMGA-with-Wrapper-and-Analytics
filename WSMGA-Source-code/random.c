/* Declaration for random number related variables and routines */

# ifndef _RAND_H_
# define _RAND_H_

/* Variable declarations for the random number generator */
extern double seed;
extern double oldrand[55];
extern int jrand;

/* Function declarations for the random number generator */
void randomize(void);
void warmup_random (double seed);
void advance_random (void);
double randomperc(void);
int rnd (int low, int high);
//double rndreal (double low, double high);

# endif

/*this file stores functions related to random number generator*/

#include "math.h"
#include "stdio.h"
 
/* variables are declared static so that they cannot conflict 
with names of   */ 
/* other global variables in other files.  See K&R, p 80, for 
scope of static */ 
double seed;
double oldrand[55]; /* Array of 55 random numbers */ 
int jrand;  /* current random number */ 
//static double rndx2;  /* used with random normal deviate */ 
//static int rndcalcflag;  /* used with random normal deviate */ 

//void advance_random(void);
/* Get seed number for random and start it up */
void randomize(void)
{
    int j1;
    for(j1=0; j1<=54; j1++)
    {
        oldrand[j1] = 0.0;
    }
    jrand=0;
    warmup_random (seed);
    return;
}

void warmup_random(double random_seed)
{
	int j1, ii; 
    double new_random, prev_random; 
 
    oldrand[54] = random_seed; 
    new_random = 0.000000001; 
    prev_random = random_seed; 
    for(j1 = 1 ; j1 <= 54; j1++) 
    { 
        ii = (21*j1)%54; 
        oldrand[ii] = new_random; 
        new_random = prev_random-new_random; 
        if(new_random<0.0) new_random = new_random + 1.0; 
        prev_random = oldrand[ii]; 
    } 
    advance_random(); 
    advance_random(); 
    advance_random(); 
    jrand = 0;
	return;
}

void advance_random(void) 
/* Create next batch of 55 random numbers */ 
{ 
    int j1; 
    double new_random; 
 
    for(j1 = 0; j1 < 24; j1++) 
    { 
        new_random = oldrand[j1] - oldrand[j1+31]; 
        if(new_random < 0.0) new_random = new_random + 1.0; 
        oldrand[j1] = new_random; 
    } 
    for(j1 = 24; j1 < 55; j1++) 
    { 
        new_random = oldrand [j1] - oldrand [j1-24]; 
        if(new_random < 0.0) new_random = new_random + 1.0; 
        oldrand[j1] = new_random; 
    } 
} 

double randomperc(void) 
/* Fetch a single random number between 0.0 and 1.0 - Subtractive Method */ 
/* See Knuth, D. (1969), v. 2 for details */ 
/* name changed from random() to avoid library conflicts on some machines*/ 
{ 
    jrand++; 
    if(jrand >= 55) 
    { 
        jrand = 1; 
        advance_random(); 
    } 
    return((double) oldrand[jrand]); 
}

/* Fetch a single random integer between low and high including the bounds */
int rnd (int low, int high)
{
    int res;
    if (low >= high)
    {
        res = low;
    }
    else
    {
        res = low + (randomperc()*(high-low+1));
        if (res > high)
        {
            res = high;
        }
    }
    return (res);
}
