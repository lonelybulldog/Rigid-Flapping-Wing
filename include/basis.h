#ifndef BASIS_H
#define BASIS_H

double Det3(double a0[], double a1[], double a2[]);

double Distance(double x1, double y1, double z1, double x2, double y2, double z2);

void Inverse(double a[3][3],double b[3][3]);

void Transpose(double a[3][3],double b[3][3]);

void transpose(double **a, double **b, int m, int n);

void matrix_multiply(double **A, double **B, double **C, int m, int n, int s);

void Multiply(double a[3][3],double b[3][3],double c[3][3]);

double DetA(double a[3][3]);

int gauss(double **a,double *b,int n);

int dcinv(double **a,int nn);

void first_segment(double cc, double *coef, double *previous);

void second_segment(double cc, double *coef);

void third_segment(double cc, double *coef);

#endif /*BASIS_H*/
