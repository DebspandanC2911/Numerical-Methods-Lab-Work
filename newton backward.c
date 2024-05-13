// C program to demonstrate 
// both Forward and Backward 
// Newton's Interpolation 
#include <stdio.h> 
void BackWard(float x[5], float y[5][5], int n); 
int main() 
{ 
    /*Enter value of x in array of x and change the first column of y correspondingly*/
	int i, j; 
	int n = 5; // number of arguments 
	float x[5] = { 0, 1, 2,3,4 }; 
	float y[5][5] = { 
		{ 0, 0, 0, 0,0 }, 
		{ 1, 0, 0, 0 ,0}, 
		{ 16, 0, 0, 0 ,0}, 
		{ 81, 0, 0, 0 ,0}, 
		{ 256, 0, 0, 0,0 }, 
		 
	}; 
	BackWard(x, y, n);
	return 0; 
} 
void BackWard(float x[5], float y[5][5], int n) 
{ 
	int i, j; 
	float a = 1.5; // interpolation point 
	float h, u, sum, p; 
	for (j = 1; j < n; j++) { 
		for (i = j; i < n; i++) { 
			y[i][j] = y[i][j - 1] - y[i - 1][j - 1]; 
		} 
	} 
	printf("\nThe backward difference table is:\n"); 
	for (i = 0; i < n; i++) { 
		printf("\n"); 
		for (j = 0; j <= i; j++) { 
			printf("%f\t", y[i][j]); 
		} 
	} 

	p = 1.0; 
	sum = y[n - 1][0]; 
	h = x[1] - x[0]; 
	u = (a - x[n - 1]) / h; 
	for (j = 1; j < n; j++) { 
		p = p * (u + j - 1) / j; 
		sum = sum + p * y[n - 1][j]; 
	} 

	printf("\nThe value of y at x=%0.1f is %0.4f", a, sum); 
}

