
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <omp.h>
#include "locale.h"

using namespace std;


int Calc_ser(double** u, double** f, int N, double eps);  
int Calc_blk(double** u, double** f, int N, double eps);  
int Calc_par(double** u, double** f, int N, double eps); 


void Init(double **u, double **f, int N);
double** new_arr(int N);
void delete_arr(double** arr, int N);

void Output(double** u, int N);


int main(int argc, char **argv)
{
	double **u=NULL, **f=NULL;
	
	const int N =600;       
	const double eps = 0.00001;   
	int icnt;                 
	double stime = -1;        
	
	f = new_arr(N);   
	u = new_arr(N + 2);  



		cout << "\n\t*** Serial version ***\n";
		Init(u, f, N);                 
    stime = omp_get_wtime();
		icnt = Calc_ser(u, f, N, eps); 
		cout << "Solution time = " << omp_get_wtime() -stime << endl;
    cout << "Iterations =    " << icnt << endl;
		cout << "Results:\n";
		Output(u, N);                   

	

		cout << "\n\t*** Block serial version ***\n";
		Init(u, f, N);                  
		stime = omp_get_wtime();
		icnt = Calc_blk(u, f, N, eps);  
		cout << "Solution time = " << omp_get_wtime() -stime << endl;
    cout << "Iterations =    " << icnt << endl;
		cout << "Results:\n";
		Output(u, N);


	
		cout << "\n\t*** Parallel version ***\n";
		Init(u, f, N);                 
		stime = omp_get_wtime();
		icnt = Calc_par(u, f, N, eps);  
		cout << "Solution time = " << omp_get_wtime() -stime << endl;
    cout << "Iterations =    " << icnt << endl;
		cout << "Results:\n";
		Output(u, N);


  delete_arr(f, N);
  delete_arr(u, N + 2);
  
	return 0;
}


int Calc_ser(double** u, double** f, int N, double eps)
{
	double max;                  
	double h = 1.0 / (N + 1);   
	int icnt = 0;             

	do
	{
		icnt++;
		max = 0;
		for (int i = 1; i <= N; i++)
			for (int j = 1; j <= N; j++)
			{
				double u0 = u[i][j];
				u[i][j] = 0.25 * (u[i - 1][j] + u[i + 1][j]   
					         + u[i][j - 1] + u[i][j + 1] - h * h * f[i - 1][j - 1]);
				double d = fabs(u[i][j] - u0);             
				if (d > max)                      
					max = d;
			}
	}
  while (max > eps);

	return icnt;
}

int Calc_blk(double** u, double** f, int N, double eps)
{
	double max_global;
	double h = 1.0 / (N + 1);
	int icnt = 0;
	const int BlockSize = 50;  
	int bcnt = N / BlockSize;

	if (N % BlockSize != 0) {
			cout << "Error: N must be divisible by BlockSize!" << endl;
			return -1;
	}

	do {
					icnt++;
					max_global = 0;
					// printf("\n--- Iteration %d ---\n\n", icnt);
					for (int wave = 0; wave < bcnt; wave++) {
					  
							for (int diag = 0; diag <= wave; diag++) {
									int i_block = diag;
									int j_block = wave - diag;
									int i_start = i_block * BlockSize + 1;
									int i_end = (i_block + 1) * BlockSize;
									int j_start = j_block * BlockSize + 1;
									int j_end = (j_block + 1) * BlockSize;
							// 		 printf("  Obrabotka bloka (%d, %d): uzli ot (%d, %d) do (%d, %d)\n",
							// i_block, j_block, i_start, j_start, i_end, j_end);
									for (int i = i_start; i <= i_end; i++) {
											for (int j = j_start; j <= j_end; j++) {
													double u_old = u[i][j];
													u[i][j] = 0.25 * (u[i-1][j] + u[i+1][j] 
																					+ u[i][j-1] + u[i][j+1] 
																					- h * h * f[i-1][j-1]);
													double diff = fabs(u[i][j] - u_old);
													if (diff > max_global) max_global = diff;
										
											}
									}
							}
					
			}
			// printf("\nSecond cycle:\n");
			for (int wave = bcnt; wave <2 * bcnt-1; wave++) {
				for (int j_block = wave - bcnt + 1; j_block < bcnt; j_block++) {
					  int i_block = wave - j_block;
						int i_start = i_block * BlockSize + 1;
						int i_end = (i_block + 1) * BlockSize;
						int j_start = j_block * BlockSize + 1;
						int j_end = (j_block + 1) * BlockSize;
						// printf("  Obrabotka bloka (%d, %d): uzli ot (%d, %d) do (%d, %d)\n",
						// 	i_block, j_block, i_start, j_start, i_end, j_end);
						for (int i = i_start; i <= i_end; i++) {
								for (int j = j_start; j <= j_end; j++) {
										double u_old = u[i][j];
										u[i][j] = 0.25 * (u[i-1][j] + u[i+1][j] 
																		+ u[i][j-1] + u[i][j+1] 
																		- h * h * f[i-1][j-1]);
										double diff = fabs(u[i][j] - u_old);
										if (diff > max_global) max_global = diff;
						
								}
						}
				}
		
}


	} while (max_global > eps);

	return icnt;
}

int Calc_par(double** u, double** f, int N, double eps) {
	double max_global;
	double h = 1.0 / (N + 1);
	int icnt = 0;
	const int BlockSize = 50;  
	int bcnt = N / BlockSize;

	if (N % BlockSize != 0) {
			cout << "Error: N must be divisible by BlockSize!" << endl;
			return -1;
	}

	do {
					icnt++;
					max_global = 0;
		
					for (int wave = 0; wave < bcnt; wave++) {
							#pragma omp parallel for reduction(max:max_global)
						  
							for (int diag = 0; diag <= wave; diag++) {
									int i_block = diag;
									int j_block = wave - diag;
									if (i_block >= bcnt || j_block < 0 || j_block >= bcnt) continue;
						
									
									int i_start = i_block * BlockSize + 1;
									int i_end = (i_block + 1) * BlockSize;
									int j_start = j_block * BlockSize + 1;
									int j_end = (j_block + 1) * BlockSize;
									double max_local=0;
									for (int i = i_start; i <= i_end; i++) {
											for (int j = j_start; j <= j_end; j++) {
													double u_old = u[i][j];
													u[i][j] = 0.25 * (u[i-1][j] + u[i+1][j] 
																					+ u[i][j-1] + u[i][j+1] 
																					- h * h * f[i-1][j-1]);
													double diff = fabs(u[i][j] - u_old);
													if (diff > max_local) max_local = diff;
											}
									}
									if(max_local  > max_global) max_global = max_local;
							}
							#pragma omp barrier
					
			}

			for (int wave = bcnt; wave <2 * bcnt-1; wave++) {
				#pragma omp parallel for reduction(max:max_global)
			 
				for (int j_block = wave - bcnt + 1; j_block < bcnt; j_block++) {
					int i_block = wave - j_block;
					if (i_block < 0 || i_block >= bcnt || j_block >= bcnt) continue;
					int i_start = i_block * BlockSize + 1;
					int i_end = (i_block + 1) * BlockSize;
					int j_start = j_block * BlockSize + 1;
					int j_end = (j_block + 1) * BlockSize;

						double max_local=0;
						for (int i = i_start; i <= i_end; i++) {
								for (int j = j_start; j <= j_end; j++) {
										double u_old = u[i][j];
										u[i][j] = 0.25 * (u[i-1][j] + u[i+1][j] 
																		+ u[i][j-1] + u[i][j+1] 
																		- h * h * f[i-1][j-1]);
										double diff = fabs(u[i][j] - u_old);
										if (diff > max_local) max_local = diff;
								}
						}
						if(max_local  > max_global) max_global = max_local;
				}
				#pragma omp barrier
	}
	} while (max_global > eps);
	return icnt;
}


double** new_arr(int N)
{
	double** f = new double* [N];
	for (int i = 0; i < N; i++)
	{
		f[i] = new double [N];
	}
	return f;
}


void delete_arr(double** arr, int N)
{
	for (int i = 0; i < N; i++)
	{
		delete[] arr[i];
	}
	delete[] arr;
}


double G(double x, double y)
{
	if (x == 0) return 1 - 2 * y;
	if (x == 1) return -1 + 2 * y;
	if (y == 0) return 1 - 2 * x;
	if (y == 1) return -1 + 2 * x;
}

double F(double x, double y)
{
	return 2.2;
}

   
void Init(double **u, double **f, int N)
{
	double h = 1.0 / (N + 1);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
			f[i][j] = F((i + 1) * h, (j + 1) * h);
	}
	for (int i = 1; i < N + 1; i++)
	{
		for (int j = 1; j < N + 1; j++)
			u[i][j] = 0.2;
		u[i][0] = G(i * h, 0);
		u[i][N + 1] = G(i * h, (N + 1) * h);
	}
	for (int j = 0; j < N + 2; j++)
	{
		u[0][j] = G(0, j * h);
		u[N + 1][j] = G((N + 1) * h, j * h);
	}
}


void Output(double** u, int N)
{
  const int K = 5;
	cout << fixed << setprecision(8);
	for (int i = 0; i <= K; i++)
	{
		for (int j = 0; j <= K; j++)
			cout << setw(12) << u[i * (N + 1) / K][j * (N + 1) / K];
		cout << endl;
	}
}

