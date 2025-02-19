#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#define N 20000

int a[N][N];

// #define SCHED_OPT static
// #define SCHED_OPT dynamic
// #define SCHED_OPT guided
#define SCHED_OPT auto
// #define SCHED_OPT runtime

int main()
{

	int i, j, k = 1;
	double start;
	start = omp_get_wtime();
	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++)
			a[i][j] = rand();

	printf("sync:		%02.3lf s.\n", omp_get_wtime() - start);
	start = omp_get_wtime();
	for (i = 0; i < N; i++) {
         #pragma omp parallel for schedule(SCHED_OPT)
					for (j = 0; j < N; j++)
						a[i][j] = rand();
	
	}

	printf("internal:	%02.3lf s.\n", omp_get_wtime() - start);
	start = omp_get_wtime();
   #pragma omp parallel for private(j) schedule(SCHED_OPT)
		for (i = 0; i < N; i++)
			for (j = 0; j < N; j++)
				a[i][j] = rand();

	printf("external:	%02.3lf s.\n", omp_get_wtime() - start);
	start = omp_get_wtime();

    #pragma omp parallel
	{
    #pragma omp for collapse(2) schedule(SCHED_OPT)
		for (i = 0; i < N; i++)
			for (j = 0; j < N; j++)
				a[i][j] = rand();
	}

	printf("collapse:	%02.3lf s.\n", omp_get_wtime() - start);

	return 0;
}
