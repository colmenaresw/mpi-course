/* 
 * parallelized trapezoidal rule with collective communication
 * only results is printed
 * now the function is passed as argument
 * Wilfredo Colmenares (2130541)
 * David Teran (2132546)
 */
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

/* we define a value for pi/4 */

#define _USE_MATH_DEFINES
const double pi_q = M_PI_4;

/* function to calculate the boundaries of each processor */
void find_my_interval(int p, int my_rank, double a, double b, double *boundaries)
{
	boundaries[0] = a + my_rank*((b - a) / p);
	boundaries[1] = a + (my_rank + 1)*((b - a) / p);
}


/* function to calculate the trapezoidal integral */
double trapezoidal(float a, float b, int n, double (*fun)(float x)) {
    double  integral;   /* Store result in integral   */
    double  h;          /* Trapezoid base width       */
    double  x;
    int    i;

    //double f(float x);  /* Function we're integrating */

    h = (b-a)/n;
    integral = (fun(a) + fun(b))/2.0;
    x = a;
    for (i = 1; i <= n-1; i++) {
        x = x + h;
        integral = integral + fun(x);
    }
    integral = integral*h;
    
	return integral;
}


/* the function we want to integrate */
double f(float x) {
    double return_val;

	/* Calculate f(x).  Store calculation in return_val. */
    return_val = 1/(1+x*x);
    return return_val;
}


double f2(float x) {
    double return_val;

	/* Calculate f(x).  Store calculation in return_val. */
    return_val = x*x;
    return return_val;
}

/* function to estimate T2N from TN
 * here a is lower limit b is upper tn is the previous integral and n is 
 * how many sub-intervals we are going to calculate
 */
double trapez_2n(float a, float b, float tn, int n, double (*fun)(float x))
{

	double t2n = 0;				// final result
	double h = (b-a) / (2 * n);	// step size

	// we calculate f(x) for each new sub interval
	for(double i = a + h; i < b; i += 2 * h)
	{
		t2n += fun(i);
	}

	// we calculate the final result
	t2n = t2n * h + tn * 0.5;
	
	return t2n;
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
	int 		n; 
	double 		a;
	double 		b;
	double		total_s2n;
	double		my_boundaries[] = {0, 0};
	double init_time, final_time, total_time;
	double temp = 999;  // to store temporary values
	double final_ans = 0;

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
		scanf("%lf %lf %d", &a, &b, &n);
	}
	
	// here starts the parallel execution
	init_time = MPI_Wtime();

	// here we distribute info to the processors
	MPI_Bcast(&a, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&b, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	//printf("hi from %d these are my a: %f b: %f and n: %d\n", my_rank, a, b, n);
	// first we assign the subintervals
	find_my_interval(p, my_rank, a, b, my_boundaries);	
    
	/* we calculate the integral with my new intervals */
	double tn = trapezoidal(my_boundaries[0], my_boundaries[1], n, f2);
	//printf("im %d and Im calculating between %f and %f and my result is: %0.8f\n", my_rank, my_boundaries[0], my_boundaries[1], tn); 
	bool stop = false;  // to know when to stop the loop
	
	while(true)
	{
		

		/* calculate the value of t2n */
		double t2n = trapez_2n(my_boundaries[0], my_boundaries[1], tn, 2 * n, f2);
		/* local sub-integral */
		double s2n = (4.0/3.0) * t2n - (1.0/3.0) * tn;
		//printf("im %d and my tn is: %0.8f and my t2n is: %0.8f\n", my_rank, tn, t2n);
		//printf("im %d and my local ans is: %0.8f\n", my_rank,  s2n);

		/* we send all the results to zero and check the accuracy of the result*/
		MPI_Reduce(&s2n, &total_s2n, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

		if(my_rank == 0)
		{	
			//printf("the reduce results is: %0.8f\n", total_s2n);
			double err = 0;
			double tol = 1e-8;
			//err = fabs(total_s2n - pi_q);
			err = fabs(total_s2n - (1.0/3.0));
			//printf("with n: %d err is: %0.10f as we are substracting %0.10f and %0.10f and tol is: %0.8f\n", n, err, total_s2n, pi_q, tol);
			if(err < tol || n == 1048576)
			{
				//printf("im inside the stoping condition\n");
				stop = true;
			}

			// the best posible results is save
			if(err < temp)
			{
				temp = err;
				final_ans = total_s2n;
			}
 
		}

		// Sincronize all processors and update boolean condition
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&stop, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);  // we tell the other processors to stop
		//printf("im %d and the value of stop is: %d\n", my_rank, stop);
		
		if(stop == true)
		{
			break;
		}
		else
		{
			n = 2 * n;  // we double the amount of intervals
			tn = t2n;   // we use the previous result
		}

		
	}
	
	if(my_rank==0)
	{
		printf("final answer is %0.9f with tolerance %0.8f \n", final_ans, temp);
	}
	/* Shut down MPI */
    MPI_Finalize();
	return 0;
} 


