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