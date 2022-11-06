/* 
 * Timing in MPI
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
	float		my_boundaries[] = {0, 0};
	double		init_time, final_time, total_time, init_sum_time, final_sum_time, total_sum_time; /* to record the perfomance */
    
	/* Start up MPI */
    MPI_Init(&argc, &argv);
	
    /* Find out process rank  */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Find out number of processes */
    MPI_Comm_size(MPI_COMM_WORLD, &p);
	
	/* from here we begin our parallel execution for all */
	MPI_Barrier(MPI_COMM_WORLD);
	init_time = MPI_Wtime();

	// we get the arguments and convert them to float
	a = atof(argv[1]);
	b = atof(argv[2]);
	n = atoi(argv[3]);

	if(my_rank == 0)
	{
		//printf("Hi! we are integrating between %f and %f with %d partitions for each processor\n", a, b, n);
	}
	
	// first we assign the subintervals
	find_my_interval(p, my_rank, a, b, my_boundaries);	
    
	/* we calculate the integral with my new intervals */
	float my_ans = trapezoidal(my_boundaries[0], my_boundaries[1], n);	

	/* we send all the results to zero */
	if(my_rank != 0)  // even node
	{
		dest = 0;
		MPI_Send(&my_ans,			  /* what am i sending */ 
				 1,				      /* size of the message */
				 MPI_FLOAT,		      /* type of message */ 
				 dest,                /* where am i sending my message */
				 tag,                 /* what's the tag of my message */
				 MPI_COMM_WORLD);     /* the communicator */
	}
	else  // im the 0 process
	{
		float final_ans = my_ans;
		float message;

		for(int q = 1; q < p; q++)
		{
			init_sum_time = MPI_Wtime();	
			source = q;
			MPI_Recv(&message,       /* where am I receiving my message */
					 1,              /* the size of the message */
					 MPI_FLOAT,      /* what type of message am i receiving */
					 source,         /* from where am i receiving */
					 tag,            /* what tag has my message */
					 MPI_COMM_WORLD, /* the communicator */
					 &status);       /* to inform what's the status of my message (output) */

			final_ans += message;
		}
		
		/* we print the output in 0 */
		//printf("the final result is %f\n with %d processes\n", final_ans, p);
		final_sum_time = MPI_Wtime();
		total_sum_time = final_sum_time - init_sum_time;
	}

	
	MPI_Barrier(MPI_COMM_WORLD);
	if(my_rank==0)
	{
		final_time = MPI_Wtime();
		total_time = final_time - init_time;	/* the result for the total execution time */
		printf("%d;%f;%f\n", p, total_time, (total_sum_time/total_time));
	}

	/* Shut down MPI */
    MPI_Finalize();
	return 0;
} 


