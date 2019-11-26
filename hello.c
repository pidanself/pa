#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main()
{
    int nthreads, tid;
    nthreads=2*omp_get_num_procs()-1;

	double ctemp,stemp,temp;
	int m,n;
	m=;
	n=;
	ctemp=12,stemp=13;
	
	for(int i=0;i<=n;i++){
		

	}

    /* Fork a team of threads giving them their own copies of variables */
    #pragma omp parallel num_threads(3) private(nthreads, tid)
    {

        /* Obtain thread number */
        tid = omp_get_thread_num();
        printf("Hello World from thread = %d\n", tid);

        /* Only master thread does this */
        if (tid == 0)
        {
            nthreads = omp_get_num_threads();
            printf("Number of threads = %d\n", nthreads);
        }

    }  /* All threads join master thread and disband */
    return 0;
}
