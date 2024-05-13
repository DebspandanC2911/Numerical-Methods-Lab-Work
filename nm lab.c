//1-Bisection Method
#include<stdio.h>
#include<math.h>
/*
 Defining equation to be solved.
 Change this equation to solve another problem.
*/
#define f(x) cos(x) - x * exp(x)

void main()
{
	 float x0, x1, x2, f0, f1, f2, e;
	 int step = 1;
	 clrscr();
	 /* Inputs */
	 up:
	 printf("\nEnter two initial guesses:\n");
	 scanf("%f%f", &x0, &x1);
	 printf("Enter tolerable error:\n");
	 scanf("%f", &e);
	 /* Calculating Functional Value */
	 f0 = f(x0);
	 f1 = f(x1);
	 /* Checking whether given guesses brackets the root or not. */
	 if( f0 * f1 > 0.0)
	 {
		  printf("Incorrect Initial Guesses.\n");
		  goto up;
	 }
   /* Implementing Bisection Method */
	 printf("\nStep\t\tx0\t\tx1\t\tx2\t\tf(x2)\n");
	 do
	 {
		  x2 = (x0 + x1)/2;
		  f2 = f(x2);
		
		  printf("%d\t\t%f\t%f\t%f\t%f\n",step, x0, x1, x2, f2);
		
		  if( f0 * f2 < 0)
		  {
			   x1 = x2;
			   f1 = f2;
		  }
		  else
		  {
			   x0 = x2;
			   f0 = f2;
		  }
		  step = step + 1;
	 }while(fabs(f2)>e);
	 printf("\nRoot is: %f", x2);
	 }


//2- Regula Falsi Method-
#include<stdio.h>
#include<math.h>
/* Defining equation to be solved.
   Change this equation to solve another problem. */
#define   f(x)   x*log10(x) - 1.2

int main()
{
	
	 float x0, x1, x2, f0, f1, f2, e;
	 int step = 1;
	 	 /* Inputs */
	 up:
	 printf("\nEnter two initial guesses:\n");
	 scanf("%f%f", &x0, &x1);
	 printf("Enter tolerable error:\n");
	 scanf("%f", &e);
	 /* Calculating Functional Values */
	 f0 = f(x0);
	 f1 = f(x1);
	 /* Checking whether given guesses brackets the root or not. */
	 if( f0*f1 > 0.0)
	 {
		  printf("Incorrect Initial Guesses.\n");
		  goto up;
	 }
	 /* Implementing Regula Falsi or False Position Method */
	 printf("\nStep\t\tx0\t\tx1\t\tx2\t\tf(x2)\n");
	 do
	 {
		  x2 = x0 - (x0-x1) * f0/(f0-f1);
		  f2 = f(x2);
		  printf("%d\t\t%f\t%f\t%f\t%f\n",step, x0, x1, x2, f2);
		
		  if(f0*f2 < 0)
		  {
			   x1 = x2;
			   f1 = f2;
		  }
		  else
		  {
			   x0 = x2;
			   f0 = f2;
		  }
		  step = step + 1;
	
	 }while(fabs(f2)>e);

	 printf("\nRoot is: %f", x2);
	 	 return 0;
}


//3-Newton Raphson Method
#include<stdio.h>
#include<math.h>
#include<stdlib.h>

/* Defining equation to be solved.
   Change this equation to solve another problem. */
#define    f(x)    3*x - cos(x) -1

/* Defining derivative of g(x).
   As you change f(x), change this function also. */
#define   g(x)   3 + sin(x)

void main()
{
	 float x0, x1, f0, f1, g0, e;
	 int step = 1, N;
	 
     /* Inputs */
	 printf("\nEnter initial guess:\n");
	 scanf("%f", &x0);
	 printf("Enter tolerable error:\n");
	 scanf("%f", &e);
	 printf("Enter maximum iteration:\n");
	 scanf("%d", &N);
	 /* Implementing Newton Raphson Method */
	 printf("\nStep\t\tx0\t\tf(x0)\t\tx1\t\tf(x1)\n");
	 do
	 {
		  g0 = g(x0);
		  f0 = f(x0);
		  if(g0 == 0.0)
		  {
			   printf("Mathematical Error.");
			   exit(0);
		  }

		
		  x1 = x0 - f0/g0;

		
		  printf("%d\t\t%f\t%f\t%f\t%f\n",step,x0,f0,x1,f1);
		  x0 = x1;
		  
		  step = step+1;
		
		  if(step > N)
		  {
			   printf("Not Convergent.");
			   exit(0);
		  }
		  
		  f1 = f(x1);
		  
	 }while(fabs(f1)>e);
	
	 printf("\nRoot is: %f", x1);
	 
}


//4-Gauss Elimination-

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#define   SIZE   10
int main()
{
	 float a[SIZE][SIZE], x[SIZE], ratio;
	 int i,j,k,n;
	 /* Inputs */
	 /* 1. Reading number of unknowns */
	 printf("Enter number of unknowns: ");
	 scanf("%d", &n);
	 /* 2. Reading Augmented Matrix */
	 for(i=1;i<=n;i++)
	 {
		  for(j=1;j<=n+1;j++)
		  {
			   printf("a[%d][%d] = ",i,j);
			   scanf("%f", &a[i][j]);
		  }
	 }
	/* Applying Gauss Elimination */
	 for(i=1;i<=n-1;i++)
	 {
		  if(a[i][i] == 0.0)
		  {
			   printf("Mathematical Error!");
			   exit(0);
		  }
		  for(j=i+1;j<=n;j++)
		  {
			   ratio = a[j][i]/a[i][i];
			   
			   for(k=1;k<=n+1;k++)
			   {
			  		a[j][k] = a[j][k] - ratio*a[i][k];
			   }
		  }
	 }
	 /* Obtaining Solution by Back Subsitution */
	 x[n] = a[n][n+1]/a[n][n];
	
	 for(i=n-1;i>=1;i--)
	 {
		  x[i] = a[i][n+1];
		  for(j=i+1;j<=n;j++)
		  {
		  		x[i] = x[i] - a[i][j]*x[j];
		  }
		  x[i] = x[i]/a[i][i];
	 }
	 /* Displaying Solution */ 
	 printf("\nSolution:\n");
	 for(i=1;i<=n;i++)
	 {
	  	printf("x[%d] = %0.3f\n",i, x[i]);
	 }
	
	 return(0);
}

//5-Gauss Siedel

#include<stdio.h>
#include<math.h>

/* Arrange systems of linear
   equations to be solved in
   diagonally dominant form
   and form equation for each
   unknown and define here
*/
/* In this example we are solving
   3x + 20y - z = -18
   2x - 3y + 20z = 25
   20x + y - 2z = 17
*/
/* Arranging given system of linear
   equations in diagonally dominant
   form:
   20x + y - 2z = 17
   3x + 20y -z = -18
   2x - 3y + 20z = 25
*/
/* Equations:
   x = (17-y+2z)/20
   y = (-18-3x+z)/20
   z = (25-2x+3y)/20
*/
/* Defining function */
#define f1(x,y,z)  (17-y+2*z)/20
#define f2(x,y,z)  (-18-3*x+z)/20
#define f3(x,y,z)  (25-2*x+3*y)/20

/* Main function */
int main()
{
 float x0=0, y0=0, z0=0, x1, y1, z1, e1, e2, e3, e;
 int count=1;

 printf("Enter tolerable error:\n");
 scanf("%f", &e);

 printf("\nCount\tx\ty\tz\n");
 do
 {
  /* Calculation */
  x1 = f1(x0,y0,z0);
  y1 = f2(x1,y0,z0);
  z1 = f3(x1,y1,z0);
  printf("%d\t%0.4f\t%0.4f\t%0.4f\n",count, x1,y1,z1);

  /* Error */
  e1 = fabs(x0-x1);
  e2 = fabs(y0-y1);
  e3 = fabs(z0-z1);

  count++;

  /* Set value for next iteration */
  x0 = x1;
  y0 = y1;
  z0 = z1;

 }while(e1>e && e2>e && e3>e);

 printf("\nSolution: x=%0.3f, y=%0.3f and z = %0.3f\n",x1,y1,z1);


 return 0;
}


//6-Gauss Jordan

#include<stdio.h>
#include<math.h>
#define SIZE 10
int main()
{
		 float a[SIZE][SIZE], x[SIZE], ratio;
		 int i,j,k,n;
		 /* Inputs */
		 /* 1. Reading number of unknowns */
		 printf("Enter number of unknowns: ");
		 scanf("%d", &n);
		 /* 2. Reading Augmented Matrix */
		 printf("Enter coefficients of Augmented Matrix:\n");
		 for(i=1;i<=n;i++)
		 {
			  for(j=1;j<=n+1;j++)
			  {
				   printf("a[%d][%d] = ",i,j);
				   scanf("%f", &a[i][j]);
			  }
		 }
		 /* Applying Gauss Jordan Elimination */
		 for(i=1;i<=n;i++)
		 {
			  if(a[i][i] == 0.0)
			  {
				   printf("Mathematical Error!");
				   exit(0);
			  }
			  for(j=1;j<=n;j++)
			  {
				   if(i!=j)
				   {
					    ratio = a[j][i]/a[i][i];
					    for(k=1;k<=n+1;k++)
					    {
					     	a[j][k] = a[j][k] - ratio*a[i][k];
					    }
				   }
			  }
		 }
		 /* Obtaining Solution */
		 for(i=1;i<=n;i++)
		 {
		  	x[i] = a[i][n+1]/a[i][i];
		 }
		 /* Displaying Solution */
		 printf("\nSolution:\n");
		 for(i=1;i<=n;i++)
		 {
		  	printf("x[%d] = %0.3f\n",i, x[i]);
		 }
		 return(0);
}

//7-Lagrange Interpolation

#include<stdio.h>
void main()
{
	 float x[100], y[100], xp, yp=0, p;
	 int i,j,n;
	 /* Input Section */
	 printf("Enter number of data: ");
	 scanf("%d", &n);
	 printf("Enter data:\n");
	 for(i=1;i<=n;i++)
	 {
		  printf("x[%d] = ", i);
		  scanf("%f", &x[i]);
		  printf("y[%d] = ", i);
		  scanf("%f", &y[i]);
	 }
	 printf("Enter interpolation point: ");
	 scanf("%f", &xp);
	 /* Implementing Lagrange Interpolation */
	 for(i=1;i<=n;i++)
	 {
		  p=1;
		  for(j=1;j<=n;j++)
		  {
			   if(i!=j)
			   {
			    	p = p* (xp - x[j])/(x[i] - x[j]);
			   }
		  }
		  yp = yp + p * y[i];
	 }
	 printf("Interpolated value at %.3f is %.3f.", xp, yp);
}

//8-Newton Divided Difference

#include<stdio.h>

int main() {
    float x[10], y[10][10], sum, p;
    int i, n, j, k = 0, f;
    
    printf("\nhow many records will you enter: ");
    scanf("%d", &n);
    
    for (i = 0; i < n; i++) {
        printf("\n\nenter the value of x%d: ", i);
        scanf("%f", &x[i]);
        printf("\n\nenter the value of f(x%d): ", i);
        scanf("%f", &y[0][i]);
    }
    
    printf("\n\nEnter X for finding f(x): ");
    scanf("%f", &p);
    
    for (i = 1; i < n; i++) {
        for (j = 0; j < n - i; j++) {
            y[i][j] = (y[i - 1][j + 1] - y[i - 1][j]) / (x[j + i] - x[j]);
        }
    }
    
    printf("\n_____________________________________________________\n");
    printf("\n  x(i)\t   y(i)\t    y1(i)    y2(i)    y3(i)    y4(i)");
    printf("\n_____________________________________________________\n");
    
    for (i = 0; i < n; i++) {
        printf("\n %.3f", x[i]);
        for (j = 0; j < n - i; j++) {
            printf("   ");
            printf(" %.3f", y[j][i]);
        }
        printf("\n");
    }
    
    i = 0;
    do {
        if (x[i] < p && p < x[i + 1])
            k = 1;
        else
            i++;
    } while (k != 1);
    
    f = i;
    
    sum = y[0][0]; // Initialize sum with y[0][0]
    float term = 1.0; // Initialize the first term of the polynomial
    for (i = 1; i < n; i++) {
        term *= (p - x[i - 1]); // Update the term
        sum += y[i][0] * term; // Add the term to the sum
    }
    
    printf("\n\n f(%.2f) = %f ", p, sum);
    return 0;
}

//9- Newton Forward Interpolation-

#include <stdio.h> 
void forward(float x[3], float y[3][3], int n); 
//void BackWard(float x[4], float y[4][4], int n); 
int main() 
{ 
    /*Enter value of x in array of x and change the first column of y correspondingly*/
	int i, j; 
	int n = 3; // number of arguments 
	float x[3] = { 0, 1, 2 }; 
	float y[3][3] = { 
		{ 0, 0, 0, 0 }, 
		{ 1, 0, 0, 0 }, 
		{ 8, 0, 0, 0 }, 
		 
	}; 

	forward(x, y, n); 
	//BackWard(x, y, n); 

	return 0; 
} 
void forward(float x[3], float y[3][3], int n) 
{ 
	int i, j; 
	float a = 1.2; // interpolation point 
	float h, u, sum, p; 
	for (j = 1; j < n; j++) { 
		for (i = 0; i < n - j; i++) { 
			y[i][j] = y[i + 1][j - 1] - y[i][j - 1]; 
		} 
	} 
	printf("\n The forward difference table is:\n"); 
	for (i = 0; i < n; i++) { 
		printf("\n"); 
		for (j = 0; j < n - i; j++) { 
			printf("%f\t", y[i][j]); 
		} 
	} 

	p = 1.0; 
	sum = y[0][0]; 
	h = x[1] - x[0]; 
	u = (a - x[0]) / h; 
	for (j = 1; j < n; j++) { 
		p = p * (u - j + 1) / j; 
		sum = sum + p * y[0][j]; 
	} 
	printf("\nThe value of y at x=%0.1f is %0.3f", a, sum); 
} 

//10-Newton Backward Interpolation-

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


//11-Trapezoidal Rule

#include<stdio.h>
#include<math.h>
/* Define function here */
#define f(x) 1/(1+pow(x,2))
int main()
{
 float lower, upper, integration=0.0, stepSize, k;
 int i, subInterval;
 /* Input */
 printf("Enter lower limit of integration: ");
 scanf("%f", &lower);
 printf("Enter upper limit of integration: ");
 scanf("%f", &upper);
 printf("Enter number of sub intervals: ");
 scanf("%d", &subInterval);

 /* Calculation */
 /* Finding step size */
 stepSize = (upper - lower)/subInterval;

 /* Finding Integration Value */
 integration = f(lower) + f(upper);
 for(i=1; i<= subInterval-1; i++)
 {
  k = lower + i*stepSize;
  integration = integration + 2 * f(k);
 }
 integration = integration * stepSize/2;
 printf("\nRequired value of integration is: %.3f", integration);
 
 return 0;
}

//12-Simpson 1/3 rd rule

#include<stdio.h>
#include<math.h>
/* Define function here */
#define f(x) 1/(1+x*x)

int main()
{
 float lower, upper, integration=0.0, stepSize, k;
 int i, subInterval;
 /* Input */
 printf("Enter lower limit of integration: ");
 scanf("%f", &lower);
 printf("Enter upper limit of integration: ");
 scanf("%f", &upper);
 printf("Enter number of sub intervals: ");
 scanf("%d", &subInterval);

 /* Calculation */
 /* Finding step size */
 stepSize = (upper - lower)/subInterval;

 /* Finding Integration Value */
 integration = f(lower) + f(upper);
 for(i=1; i<= subInterval-1; i++)
 {
  k = lower + i*stepSize;
  if(i%2==0)
  {
   integration = integration + 2 * f(k);
  }
  else
  {
   integration = integration + 4 * f(k);
  }
 }
 integration = integration * stepSize/3;
 printf("\nRequired value of integration is: %.3f", integration);

 return 0;
}

//Simpson 3/8 rule

#include<stdio.h>
#include<math.h>
/* Define function here */
#define f(x) 1/(1+x*x)

int main()
{
 float lower, upper, integration=0.0, stepSize, k;
 int i, subInterval;
 /* Input */
 printf("Enter lower limit of integration: ");
 scanf("%f", &lower);
 printf("Enter upper limit of integration: ");
 scanf("%f", &upper);
 printf("Enter number of sub intervals: ");
 scanf("%d", &subInterval);

 /* Calculation */
 /* Finding step size */
 stepSize = (upper - lower)/subInterval;

 /* Finding Integration Value */
 integration = f(lower) + f(upper);
 for(i=1; i<= subInterval-1; i++)
 {
  k = lower + i*stepSize;
  if(i%3 == 0)
  {
   integration = integration + 2 * f(k);
  }
  else
  {
   integration = integration + 3 * f(k);
  }
 }
 integration = integration * stepSize*3/8;
 printf("\nRequired value of integration is: %.6f", integration);
  return 0;
}

//All the Best