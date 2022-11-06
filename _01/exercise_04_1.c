/* 
 * parallelized trapezoidal rule with collective communication
 * Wilfredo Colmenares (2130541)
 * David Teran (2132546)
 */
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include <stdlib.h>

/* we define a value for pi/4 */
const double pi_q = 22.0 / 7.0 / 4.0;

/* function to calculate the boundaries of each processor */
void find_my_interval(int p, int my_rank, float a, float b, float *boundaries)
{
	boundaries[0] = a + my_rank*((b - a) / p);
	boundaries[1] = a + (my_rank + 1)*((b - a) / p);
}


/* function to calculate the trapezoidal integral */
float trapezoidal(float a, float b, int n) {
    float  integral;   /* Store result in integral   */
    float  h;          /* Trapezoid base width       */
    float  x;
    int    i;

    float f(float x);  /* Function we're integrating */

    h = (b-a)/n;
    integral = (f(a) + f(b))/2.0;
    x = a;
    for (i = 1; i <= n-1; i++) {
        x = x + h;
        integral = integral + f(x);
    }
    integral = integral*h;
    
	return integral;
}


/* the function we want to integrate */
float f(float x) {
    float return_val;

	/* Calculate f(x).  Store calculation in return_val. */
    return_val = 1/(1+x*x);
    return return_val;
}


/* main function */
int main(int argc, char* argv[]) {
   /* 
 	*  arc : # of arguments from the console
 	*  argv[] : array of strings coming from the console
 	*
 	*/
    int         my_rank;       /* rank of process      */
    int         p;             /* number of processes  */
    int         source;        /* rank of sender       */
    int         dest;          /* rank of receiver     */
    int         tag = 0;       /* tag for messages     */
    char        message[100];  /* storage for message  */
    MPI_Status  status;        /* return status for receive */
	char		temp[100]; 
	int 		n; 
	float 		a;
	float 		b;
	float		ans;
	float		my_boundaries[] = {0, 0};
	double init_time, final_time, total_time;
    /* Start up MPI */
    MPI_Init(&argc, &argv);

    /* Find out process rank  */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Find out number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &p);

	// we get the arguments and send them to processors
	if(my_rank==0)
	{
		//printf("Enter the 'a' 'b' points and the 'n' trapezoids\n");
		// when using scanf the first argument specify the type and the second
		// where do I want to store my result (in term of memory address)
		scanf("%f %f %d", &a, &b, &n);
	}
	
	// here starts the parallel execution
	init_time = MPI_Wtime();

	// here we distribute info to the processors
	MPI_Bcast(&a, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&b, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	//printf("hi from %d these are my a: %f b: %f and n: %d\n", my_rank, a, b, n);
	// first we assign the subintervals
	find_my_interval(p, my_rank, a, b, my_boundaries);	
    
	/* we calculate the integral with my new intervals */
	float my_ans = trapezoidal(my_boundaries[0], my_boundaries[1], n);	

	/* we send all the results to zero and print */
	MPI_Reduce(&my_ans, &ans, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
	
	if(my_rank==0)
	{
		//printf("final answer is %f\n", ans);
		// here we determine the time
		final_time = MPI_Wtime();
		total_time = final_time - init_time;	/* the result for the total execution time */
		printf("%d;%f\n", p, total_time);
	}
	/* Shut down MPI */
    MPI_Finalize();
	return 0;
} 


