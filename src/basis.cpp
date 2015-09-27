#include <math.h>
#include <iostream>
#include <stdio.h>
#ifdef _OPENMP
#	include <omp.h>
#endif

using namespace std;

#define  PI 3.1415926535897932

double Det3(double a0[], double a1[], double a2[])
{
	return(a0[0]*(a2[2]*a1[1]-a2[1]*a1[2]) - a1[0]*(a2[2]*a0[1]-a2[1]*a0[2]) + a2[0]*(a1[2]*a0[1]-a1[1]*a0[2]));
}

double Distance(double x1, double y1, double z1, double x2, double y2, double z2)
{
	return sqrt(pow(x1-x2,2)+pow(y1-y2,2)+pow(z1-z2,2));
}


void Inverse(double a[3][3],double b[3][3])
{
	double detA;

	detA=a[0][0]*a[1][1]*a[2][2]-a[0][0]*a[1][2]*a[2][1]-a[0][1]*a[1][0]*a[2][2]
      +a[0][1]*a[1][2]*a[2][0]+a[0][2]*a[1][0]*a[2][1]-a[0][2]*a[1][1]*a[2][0];

    if(fabs(detA)<1e-20) 
	{
		cout<<"Ill Matrix for body inertia, can not find its inversion!"<<endl;
//		exit(0);
	}
	b[0][0]=(a[1][1]*a[2][2]-a[1][2]*a[2][1])/detA;
	b[0][1]=(a[0][2]*a[2][1]-a[0][1]*a[2][2])/detA;
	b[0][2]=(a[0][1]*a[1][2]-a[0][2]*a[1][1])/detA;
	
	b[1][0]=(a[1][2]*a[2][0]-a[1][0]*a[2][2])/detA;
	b[1][1]=(a[0][0]*a[2][2]-a[0][2]*a[2][0])/detA;
	b[1][2]=(a[0][2]*a[1][0]-a[0][0]*a[1][2])/detA;
	
	b[2][0]=(a[1][0]*a[2][1]-a[1][1]*a[2][0])/detA;
	b[2][1]=(a[0][1]*a[2][0]-a[0][0]*a[2][1])/detA;
	b[2][2]=(a[0][0]*a[1][1]-a[0][1]*a[1][0])/detA;

}

void Transpose(double a[3][3],double b[3][3])
{
	
    int i,j;
#	pragma omp parallel for private(j)
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
			b[i][j]=a[j][i];
}

void transpose(double **a, double **b, int m, int n)
{
	int i,j;
#	pragma omp parallel for private(j)
	for(i=0;i<m;i++)
	{
		for(j=0;j<n;j++)
		{
			b[j][i]=a[i][j];
		}
	}
}

void matrix_multiply(double **A, double **B, double **C, int m, int n, int s)
{
#	pragma omp parallel for
	for(int i=0;i<m;i++)
	{
		for(int j=0;j<s;j++)
		{
			C[i][j]=0;
			for(int l=0;l<n;l++)
			{
				C[i][j]=C[i][j]+A[i][l]*B[l][j];
			}
		}
	}
}


void Multiply(double a[3][3],double b[3][3],double c[3][3])
{
	int i,j;
	for(i=0;i<3;i++)
		for(j=0;j<3;j++)
		{
			c[i][j]= a[i][0]*b[0][j]+a[i][1]*b[1][j]+a[i][2]*b[2][j];
		}	
}

double DetA(double a[3][3])
{
	double detA;

	detA=a[0][0]*a[1][1]*a[2][2]-a[0][0]*a[1][2]*a[2][1]-a[0][1]*a[1][0]*a[2][2]
       +a[0][1]*a[1][2]*a[2][0]+a[0][2]*a[1][0]*a[2][1]-a[0][2]*a[1][1]*a[2][0];
	return detA;
}

int gauss(double **a,double *b,int n) 
{ 
	int *js; 
	int l,k,i,j,is; 
	double d,t; 
 	js=new int[n]; 
	l=1;
	for (k=0;k<=n-2;k++) {
		d=0.0; 
		for (i=k;i<=n-1;i++) 
		for (j=k;j<=n-1;j++) {
			t=fabs(*(*(a+i)+j));    
			if (t>d) { d=t; js[k]=j; is=i;} 
		} 
		if (d+1.0==1.0) 
			l=0; 
		else {
			if (js[k]!=k) 
				for (i=0;i<=n-1;i++) {  
				t=*(*(a+i)+k);*(*(a+i)+k)=*(*(a+i)+js[k]) ; *(*(a+i)+js[k])=t; 
				} 
			if (is!=k) { 
				for (j=k;j<=n-1;j++) {
				t=*(*(a+k)+j) ; *(*(a+k)+j)=*(*(a+is)+j); *(*(a+is)+j)=t; 
				} 
				t=b[k]; b[k]=b[is]; b[is]=t; 
			} 
		} 
		if (l==0) {
			delete[] js; printf("fail\n"); 
			return(0); 
		} 
		d=*(*(a+k)+k); 
		for (j=k+1;j<=n-1;j++) { 
			*(*(a+k)+j)=*(*(a+k)+j)/d;
		} 
		b[k]=b[k]/d; 
		for (i=k+1;i<=n-1;i++) {
			for (j=k+1;j<=n-1;j++) { 
				*(*(a+i)+j)=*(*(a+i)+j)-*(*(a+i)+k)**(*(a+k)+j); 
			} 
			b[i]=b[i]-*(*(a+i)+k)*b[k]; 
		} 
	}

	d=*(*(a+n-1)+n-1);
	if (fabs(d)+1.0==1.0) {
		delete[] js; printf("fail\n"); 
		return(0); 
	} 
	b[n-1]=b[n-1]/d; 
	for (i=n-2;i>=0;i--) {
		t=0.0; 
		for (j=i+1;j<=n-1;j++) 
			t=t+*(*(a+i)+j)*b[j]; 
		b[i]=b[i]-t; 
	} 
	js[n-1]=n-1; 
	for (k=n-1;k>=0;k--) 
		if (js[k]!=k) 
		{ t=b[k]; b[k]=b[js[k]]; b[js[k]]=t;} 
	delete(js); 
	return(1); 
} 


int dcinv(double **a,int nn)
{ 
	int n=nn;
int *is,*js,i,j,k;
double d,p;
is=new int[n];
js=new int[n];
for (k=0; k<=n-1; k++)
{ 
  d=0.0;
  for (i=k; i<=n-1; i++)
  {
   for (j=k; j<=n-1; j++)
   { 
    p=fabs(*(*(a+i)+j));
    if (p>d) 
    { 
     d=p; 
     is[k]=i; 
     js[k]=j;
    }
   }
  }
  if (d+1.0==1.0)
  {
   delete []is;
   delete []js;
   printf("err***not inv\n");
   return(0);
  }
  if (is[k]!=k)
  {
   for (j=0; j<=n-1; j++)
   { 
    p=*(*(a+k)+j);*(*(a+k)+j)=*(*(a+is[k])+j);*(*(a+is[k])+j)=p;
   }
  }
  if (js[k]!=k)
  {
   for (i=0; i<=n-1; i++)
   { 
    p=*(*(a+i)+k);*(*(a+i)+k)=*(*(a+i)+js[k]);*(*(a+i)+js[k])=p;
   }
  }
  *(*(a+k)+k)=1.0/(*(*(a+k)+k));
  for (j=0; j<=n-1; j++)
  {
   if (j!=k)
   { 
   *(*(a+k)+j)=*(*(a+k)+j)**(*(a+k)+k);
   }
  }
  for (i=0; i<=n-1; i++)
  {
   if (i!=k)
   {
    for (j=0; j<=n-1; j++)
    {
     if (j!=k)
     { 
		*(*(a+i)+j)=*(*(a+i)+j)-*(*(a+i)+k)**(*(a+k)+j);
     }
    }
   }
  }
  for (i=0; i<=n-1; i++)
  {
   if (i!=k)
   {    
		*(*(a+i)+k)=-*(*(a+i)+k)**(*(a+k)+k);
   }   
  } 
}
for (k=n-1; k>=0; k--)
{ 
  if (js[k]!=k)
  {
   for (j=0; j<=n-1; j++)
   { 
    p=*(*(a+k)+j); *(*(a+k)+j)=*(*(a+js[k])+j); *(*(a+js[k])+j)=p;
   }
  }
  if (is[k]!=k)
  {
   for (i=0; i<=n-1; i++)
   { 
    p=*(*(a+i)+k); *(*(a+i)+k)=*(*(a+i)+is[k]); *(*(a+i)+is[k])=p;
   }
  }
}
delete []is;
delete []js;
return(1);
}

/*
void first_segment(double cc, double ss, double delta, double *coef, double *previous)
{
    double **a, *b;
    a=new double* [8]; for(int i=0;i<8;i++) a[i]=new double [8];
    b=new double [8];
    double x1=0,x2=delta;
  
    a[0][0]=pow(x1,7);     a[0][1]=pow(x1,6);     a[0][2]=pow(x1,5);     a[0][3]=pow(x1,4);     a[0][4]=pow(x1,3);     a[0][5]=pow(x1,2);     a[0][6]=x1;     a[0][7]=1;
    a[1][0]=7*pow(x1,6);   a[1][1]=6*pow(x1,5);   a[1][2]=5*pow(x1,4);   a[1][3]=4*pow(x1,3);   a[1][4]=3*pow(x1,2);   a[1][5]=2*x1;          a[1][6]=1;      a[1][7]=0;
    a[2][0]=42*pow(x1,5);  a[2][1]=30*pow(x1,4);  a[2][2]=20*pow(x1,3);  a[2][3]=12*pow(x1,2);  a[2][4]=6*x1;          a[2][5]=2;             a[2][6]=0;      a[2][7]=0;
    a[3][0]=210*pow(x1,4); a[3][1]=120*pow(x1,3); a[3][2]=60*pow(x1,2);  a[3][3]=24*x1;         a[3][4]=6;             a[3][5]=0;             a[3][6]=0;      a[3][7]=0;

    a[4][0]=pow(x2,7);     a[4][1]=pow(x2,6);     a[4][2]=pow(x2,5);     a[4][3]=pow(x2,4);     a[4][4]=pow(x2,3);     a[4][5]=pow(x2,2);     a[4][6]=x2;     a[4][7]=1;
    a[5][0]=7*pow(x2,6);   a[5][1]=6*pow(x2,5);   a[5][2]=5*pow(x2,4);   a[5][3]=4*pow(x2,3);   a[5][4]=3*pow(x2,2);   a[5][5]=2*x2;          a[5][6]=1;      a[5][7]=0;
    a[6][0]=42*pow(x2,5);  a[6][1]=30*pow(x2,4);  a[6][2]=20*pow(x2,3);  a[6][3]=12*pow(x2,2);  a[6][4]=6*x2;          a[6][5]=2;             a[6][6]=0;      a[6][7]=0;
    a[7][0]=210*pow(x2,4); a[7][1]=120*pow(x2,3); a[7][2]=60*pow(x2,2);  a[7][3]=24*x2;         a[7][4]=6;             a[7][5]=0;             a[7][6]=0;      a[7][7]=0;
    
    b[0]=PI/2;
    b[1]=4*previous[0]+3*previous[1]+2*previous[2]+previous[3];
    b[2]=12*previous[0]+6*previous[1]+2*previous[2];
    b[3]=24*previous[0]+6*previous[1];
    
    b[4]=PI/2-(PI/4+cc)*sin(2*PI*x2/(1-4*ss));
    b[5]=-(2*PI/(1-4*ss))*(PI/4+cc)*cos(2*PI*x2/(1-4*ss));
    b[6]=pow(2*PI/(1-4*ss),2)*(PI/4+cc)*sin(2*PI*x2/(1-4*ss));
    b[7]=pow(2*PI/(1-4*ss),3)*(PI/4+cc)*cos(2*PI*x2/(1-4*ss));
    
    gauss(a,b,8);
    
    for(int i=0;i<8;i++)
    {
	coef[i]=b[i];
	delete[] a[i];
    }
    delete[] a;
    delete[] b;
}

void second_segment(double cc, double ss, double delta, double *coef)
{
    double **a, *b;
    a=new double* [9]; for(int i=0;i<9;i++) a[i]=new double [9];
    b=new double [9];
    double x1=0.25-ss-delta, x2=0.25-ss, x3=0.25-ss+delta;
    
    a[0][0]=pow(x1,8);     a[0][1]=pow(x1,7);     a[0][2]=pow(x1,6);     a[0][3]=pow(x1,5);     a[0][4]=pow(x1,4);     a[0][5]=pow(x1,3);     a[0][6]=pow(x1,2);  a[0][7]=x1;a[0][8]=1;
    a[1][0]=8*pow(x1,7);   a[1][1]=7*pow(x1,6);   a[1][2]=6*pow(x1,5);   a[1][3]=5*pow(x1,4);   a[1][4]=4*pow(x1,3);   a[1][5]=3*pow(x1,2);   a[1][6]=2*x1;       a[1][7]=1;a[1][8]=0;
    a[2][0]=56*pow(x1,6);  a[2][1]=42*pow(x1,5);  a[2][2]=30*pow(x1,4);  a[2][3]=20*pow(x1,3);  a[2][4]=12*pow(x1,2);  a[2][5]=6*x1;          a[2][6]=2;          a[2][7]=0;a[2][8]=0;
    a[3][0]=336*pow(x1,5); a[3][1]=210*pow(x1,4); a[3][2]=120*pow(x1,3); a[3][3]=60*pow(x1,2);  a[3][4]=24*x1;         a[3][5]=6;             a[3][6]=0;          a[3][7]=0;a[3][8]=0;
	
    a[4][0]=pow(x2,8);     a[4][1]=pow(x2,7);     a[4][2]=pow(x2,6);     a[4][3]=pow(x2,5);     a[4][4]=pow(x2,4);     a[4][5]=pow(x2,3);     a[4][6]=pow(x2,2);  a[4][7]=x2;a[4][8]=1;
	
    a[5][0]=pow(x3,8);     a[5][1]=pow(x3,7);     a[5][2]=pow(x3,6);     a[5][3]=pow(x3,5);     a[5][4]=pow(x3,4);     a[5][5]=pow(x3,3);     a[5][6]=pow(x3,2);  a[5][7]=x3;a[5][8]=1;
    a[6][0]=8*pow(x3,7);   a[6][1]=7*pow(x3,6);   a[6][2]=6*pow(x3,5);   a[6][3]=5*pow(x3,4);   a[6][4]=4*pow(x3,3);   a[6][5]=3*pow(x3,2);   a[6][6]=2*x3;       a[6][7]=1;a[6][8]=0;
    a[7][0]=56*pow(x3,6);  a[7][1]=42*pow(x3,5);  a[7][2]=30*pow(x3,4);  a[7][3]=20*pow(x3,3);  a[7][4]=12*pow(x3,2);  a[7][5]=6*x3;          a[7][6]=2;          a[7][7]=0;a[7][8]=0;
    a[8][0]=336*pow(x3,5); a[8][1]=210*pow(x3,4); a[8][2]=120*pow(x3,3); a[8][3]=60*pow(x3,2);  a[8][4]=24*x3;         a[8][5]=6;             a[8][6]=0;          a[8][7]=0;a[8][8]=0;
  
    b[0]=PI/2-(PI/4+cc)*sin(2*PI*x1/(1-4*ss));
    b[1]=-(2*PI/(1-4*ss))*(PI/4+cc)*cos(2*PI*x1/(1-4*ss));
    b[2]=pow(2*PI/(1-4*ss),2)*(PI/4+cc)*sin(2*PI*x1/(1-4*ss));
    b[3]=pow(2*PI/(1-4*ss),3)*(PI/4+cc)*cos(2*PI*x1/(1-4*ss));
    
    b[4]=PI/4-cc;
    
    b[5]=PI/2-(PI/4+cc)*sin(2*PI*(x3+2*ss)/(1+4*ss));
    b[6]=-(2*PI/(1+4*ss))*(PI/4+cc)*cos(2*PI*(x3+2*ss)/(1+4*ss));
    b[7]=pow(2*PI/(1+4*ss),2)*(PI/4+cc)*sin(2*PI*(x3+2*ss)/(1+4*ss));
    b[8]=pow(2*PI/(1+4*ss),3)*(PI/4+cc)*cos(2*PI*(x3+2*ss)/(1+4*ss));
    
    gauss(a,b,9);
    
    for(int i=0;i<9;i++)
    {
	coef[i]=b[i];
	delete[] a[i];
    }
    delete[] a;
    delete[] b;
}

void third_segment(double cc, double ss, double delta, double *coef)
{
    double **a, *b;
    a=new double* [9]; for(int i=0;i<9;i++) a[i]=new double [9];
    b=new double [9];
    double x1=0.5-delta, x2=0.5, x3=0.5+delta;
    
    a[0][0]=pow(x1,8);     a[0][1]=pow(x1,7);     a[0][2]=pow(x1,6);     a[0][3]=pow(x1,5);     a[0][4]=pow(x1,4);     a[0][5]=pow(x1,3);     a[0][6]=pow(x1,2);  a[0][7]=x1;a[0][8]=1;
    a[1][0]=8*pow(x1,7);   a[1][1]=7*pow(x1,6);   a[1][2]=6*pow(x1,5);   a[1][3]=5*pow(x1,4);   a[1][4]=4*pow(x1,3);   a[1][5]=3*pow(x1,2);   a[1][6]=2*x1;       a[1][7]=1;a[1][8]=0;
    a[2][0]=56*pow(x1,6);  a[2][1]=42*pow(x1,5);  a[2][2]=30*pow(x1,4);  a[2][3]=20*pow(x1,3);  a[2][4]=12*pow(x1,2);  a[2][5]=6*x1;          a[2][6]=2;          a[2][7]=0;a[2][8]=0;
    a[3][0]=336*pow(x1,5); a[3][1]=210*pow(x1,4); a[3][2]=120*pow(x1,3); a[3][3]=60*pow(x1,2);  a[3][4]=24*x1;         a[3][5]=6;             a[3][6]=0;          a[3][7]=0;a[3][8]=0;
	
    a[4][0]=pow(x2,8);     a[4][1]=pow(x2,7);     a[4][2]=pow(x2,6);     a[4][3]=pow(x2,5);     a[4][4]=pow(x2,4);     a[4][5]=pow(x2,3);     a[4][6]=pow(x2,2);  a[4][7]=x2;a[4][8]=1;
	
    a[5][0]=pow(x3,8);     a[5][1]=pow(x3,7);     a[5][2]=pow(x3,6);     a[5][3]=pow(x3,5);     a[5][4]=pow(x3,4);     a[5][5]=pow(x3,3);     a[5][6]=pow(x3,2);  a[5][7]=x3;a[5][8]=1;
    a[6][0]=8*pow(x3,7);   a[6][1]=7*pow(x3,6);   a[6][2]=6*pow(x3,5);   a[6][3]=5*pow(x3,4);   a[6][4]=4*pow(x3,3);   a[6][5]=3*pow(x3,2);   a[6][6]=2*x3;       a[6][7]=1;a[6][8]=0;
    a[7][0]=56*pow(x3,6);  a[7][1]=42*pow(x3,5);  a[7][2]=30*pow(x3,4);  a[7][3]=20*pow(x3,3);  a[7][4]=12*pow(x3,2);  a[7][5]=6*x3;          a[7][6]=2;          a[7][7]=0;a[7][8]=0;
    a[8][0]=336*pow(x3,5); a[8][1]=210*pow(x3,4); a[8][2]=120*pow(x3,3); a[8][3]=60*pow(x3,2);  a[8][4]=24*x3;         a[8][5]=6;             a[8][6]=0;          a[8][7]=0;a[8][8]=0;

    b[0]=PI/2-(PI/4+cc)*sin(2*PI*(x1+2*ss)/(1+4*ss));
    b[1]=-(2*PI/(1+4*ss))*(PI/4+cc)*cos(2*PI*(x1+2*ss)/(1+4*ss));
    b[2]=pow(2*PI/(1+4*ss),2)*(PI/4+cc)*sin(2*PI*(x1+2*ss)/(1+4*ss));
    b[3]=pow(2*PI/(1+4*ss),3)*(PI/4+cc)*cos(2*PI*(x1+2*ss)/(1+4*ss));
	
    b[4]=PI/2;
	
    b[5]=PI/2-(PI/4-cc)*sin(2*PI*(x3+2*ss)/(1+4*ss));
    b[6]=-(2*PI/(1+4*ss))*(PI/4-cc)*cos(2*PI*(x3+2*ss)/(1+4*ss));
    b[7]=pow(2*PI/(1+4*ss),2)*(PI/4-cc)*sin(2*PI*(x3+2*ss)/(1+4*ss));
    b[8]=pow(2*PI/(1+4*ss),3)*(PI/4-cc)*cos(2*PI*(x3+2*ss)/(1+4*ss));
    
    gauss(a,b,9);
    
    for(int i=0;i<9;i++)
    {
	coef[i]=b[i];
	delete[] a[i];
    }
    delete[] a;
    delete[] b;
}

void fourth_segment(double cc, double ss, double delta, double *coef)
{
    double **a, *b;
    a=new double* [9]; for(int i=0;i<9;i++) a[i]=new double [9];
    b=new double [9];
    double x1=0.75+ss-delta, x2=0.75+ss, x3=0.75+ss+delta;
	
    a[0][0]=pow(x1,8);     a[0][1]=pow(x1,7);     a[0][2]=pow(x1,6);     a[0][3]=pow(x1,5);     a[0][4]=pow(x1,4);     a[0][5]=pow(x1,3);     a[0][6]=pow(x1,2);  a[0][7]=x1;a[0][8]=1;
    a[1][0]=8*pow(x1,7);   a[1][1]=7*pow(x1,6);   a[1][2]=6*pow(x1,5);   a[1][3]=5*pow(x1,4);   a[1][4]=4*pow(x1,3);   a[1][5]=3*pow(x1,2);   a[1][6]=2*x1;       a[1][7]=1;a[1][8]=0;
    a[2][0]=56*pow(x1,6);  a[2][1]=42*pow(x1,5);  a[2][2]=30*pow(x1,4);  a[2][3]=20*pow(x1,3);  a[2][4]=12*pow(x1,2);  a[2][5]=6*x1;          a[2][6]=2;          a[2][7]=0;a[2][8]=0;
    a[3][0]=336*pow(x1,5); a[3][1]=210*pow(x1,4); a[3][2]=120*pow(x1,3); a[3][3]=60*pow(x1,2);  a[3][4]=24*x1;         a[3][5]=6;             a[3][6]=0;          a[3][7]=0;a[3][8]=0;
	
    a[4][0]=pow(x2,8);     a[4][1]=pow(x2,7);     a[4][2]=pow(x2,6);     a[4][3]=pow(x2,5);     a[4][4]=pow(x2,4);     a[4][5]=pow(x2,3);     a[4][6]=pow(x2,2);  a[4][7]=x2;a[4][8]=1;
	
    a[5][0]=pow(x3,8);     a[5][1]=pow(x3,7);     a[5][2]=pow(x3,6);     a[5][3]=pow(x3,5);     a[5][4]=pow(x3,4);     a[5][5]=pow(x3,3);     a[5][6]=pow(x3,2);  a[5][7]=x3;a[5][8]=1;
    a[6][0]=8*pow(x3,7);   a[6][1]=7*pow(x3,6);   a[6][2]=6*pow(x3,5);   a[6][3]=5*pow(x3,4);   a[6][4]=4*pow(x3,3);   a[6][5]=3*pow(x3,2);   a[6][6]=2*x3;       a[6][7]=1;a[6][8]=0;
    a[7][0]=56*pow(x3,6);  a[7][1]=42*pow(x3,5);  a[7][2]=30*pow(x3,4);  a[7][3]=20*pow(x3,3);  a[7][4]=12*pow(x3,2);  a[7][5]=6*x3;          a[7][6]=2;          a[7][7]=0;a[7][8]=0;
    a[8][0]=336*pow(x3,5); a[8][1]=210*pow(x3,4); a[8][2]=120*pow(x3,3); a[8][3]=60*pow(x3,2);  a[8][4]=24*x3;         a[8][5]=6;             a[8][6]=0;          a[8][7]=0;a[8][8]=0;
	
    b[0]=PI/2-(PI/4-cc)*sin(2*PI*(x1+2*ss)/(1+4*ss));
    b[1]=-(2*PI/(1+4*ss))*(PI/4-cc)*cos(2*PI*(x1+2*ss)/(1+4*ss));
    b[2]=pow(2*PI/(1+4*ss),2)*(PI/4-cc)*sin(2*PI*(x1+2*ss)/(1+4*ss));
    b[3]=pow(2*PI/(1+4*ss),3)*(PI/4-cc)*cos(2*PI*(x1+2*ss)/(1+4*ss));
	
    b[4]=3*PI/4-cc;
	
    b[5]=PI/2-(PI/4-cc)*sin(2*PI*(x3-4*ss)/(1-4*ss));
    b[6]=-(2*PI/(1-4*ss))*(PI/4-cc)*cos(2*PI*(x3-4*ss)/(1-4*ss));
    b[7]=pow(2*PI/(1-4*ss),2)*(PI/4-cc)*sin(2*PI*(x3-4*ss)/(1-4*ss));
    b[8]=pow(2*PI/(1-4*ss),3)*(PI/4-cc)*cos(2*PI*(x3-4*ss)/(1-4*ss));
	
    gauss(a,b,9);
    
    for(int i=0;i<9;i++)
    {
	coef[i]=b[i];
	delete[] a[i];
    }
    delete[] a;
    delete[] b;
}

void fifth_segment(double cc, double ss, double delta, double *coef)
{
    double **a,*b;
    a=new double* [5]; for(int i=0;i<5;i++) a[i]=new double [5];
    b=new double [5];
    double x1=1.0-delta,x2=1.0;
      
    a[0][0]=pow(x1,4);     a[0][1]=pow(x1,3);     a[0][2]=pow(x1,2);     a[0][3]=x1;            a[0][4]=1;
    a[1][0]=4*pow(x1,3);   a[1][1]=3*pow(x1,2);   a[1][2]=2*x1;          a[1][3]=1;             a[1][4]=0;
    a[2][0]=12*pow(x1,2);  a[2][1]=6*x1;          a[2][2]=2;             a[2][3]=0;             a[2][4]=0;
    a[3][0]=24*x1;         a[3][1]=6;             a[3][2]=0;             a[3][3]=0;             a[3][4]=0;
	
    a[4][0]=pow(x2,4);     a[4][1]=pow(x2,3);     a[4][2]=pow(x2,2);     a[4][3]=x2;            a[4][4]=0;
  
    
    b[0]=PI/2-(PI/4-cc)*sin(2*PI*(x1-4*ss)/(1-4*ss));
    b[1]=-(2*PI/(1-4*ss))*(PI/4-cc)*cos(2*PI*(x1-4*ss)/(1-4*ss));
    b[2]=pow(2*PI/(1-4*ss),2)*(PI/4-cc)*sin(2*PI*(x1-4*ss)/(1-4*ss));
    b[3]=pow(2*PI/(1-4*ss),3)*(PI/4-cc)*cos(2*PI*(x1-4*ss)/(1-4*ss));
    
    b[4]=PI/2.0;
    
    gauss(a,b,5);
	
    for(int i=0;i<5;i++)
    {
	coef[i]=b[i];
	delete[] a[i];
    }
    delete[] a;
    delete[] b;
}
*/  



void first_segment(double cc, double *coef, double *previous)
{
    double **a, *b;
    a=new double* [8]; for(int i=0;i<8;i++) a[i]=new double [8];
    b=new double [8];
    double x1=0,x2=0.1;
    
    a[0][0]=pow(x1,7);     a[0][1]=pow(x1,6);     a[0][2]=pow(x1,5);     a[0][3]=pow(x1,4);     a[0][4]=pow(x1,3);     a[0][5]=pow(x1,2);     a[0][6]=x1;     a[0][7]=1;
    a[1][0]=7*pow(x1,6);   a[1][1]=6*pow(x1,5);   a[1][2]=5*pow(x1,4);   a[1][3]=4*pow(x1,3);   a[1][4]=3*pow(x1,2);   a[1][5]=2*x1;          a[1][6]=1;      a[1][7]=0;
    a[2][0]=42*pow(x1,5);  a[2][1]=30*pow(x1,4);  a[2][2]=20*pow(x1,3);  a[2][3]=12*pow(x1,2);  a[2][4]=6*x1;          a[2][5]=2;             a[2][6]=0;      a[2][7]=0;
    a[3][0]=210*pow(x1,4); a[3][1]=120*pow(x1,3); a[3][2]=60*pow(x1,2);  a[3][3]=24*x1;         a[3][4]=6;             a[3][5]=0;             a[3][6]=0;      a[3][7]=0;

    a[4][0]=pow(x2,7);     a[4][1]=pow(x2,6);     a[4][2]=pow(x2,5);     a[4][3]=pow(x2,4);     a[4][4]=pow(x2,3);     a[4][5]=pow(x2,2);     a[4][6]=x2;     a[4][7]=1;
    a[5][0]=7*pow(x2,6);   a[5][1]=6*pow(x2,5);   a[5][2]=5*pow(x2,4);   a[5][3]=4*pow(x2,3);   a[5][4]=3*pow(x2,2);   a[5][5]=2*x2;          a[5][6]=1;      a[5][7]=0;
    a[6][0]=42*pow(x2,5);  a[6][1]=30*pow(x2,4);  a[6][2]=20*pow(x2,3);  a[6][3]=12*pow(x2,2);  a[6][4]=6*x2;          a[6][5]=2;             a[6][6]=0;      a[6][7]=0;
    a[7][0]=210*pow(x2,4); a[7][1]=120*pow(x2,3); a[7][2]=60*pow(x2,2);  a[7][3]=24*x2;         a[7][4]=6;             a[7][5]=0;             a[7][6]=0;      a[7][7]=0;
      

    b[0]=PI/2;
    b[1]=4*previous[0]+3*previous[1]+2*previous[2]+previous[3];
    b[2]=12*previous[0]+6*previous[1]+2*previous[2];
    b[3]=24*previous[0]+6*previous[1];
    b[4]=PI/2.0-(PI/4+cc)*sin(2*PI*x2);b[5]=-2*PI*(PI/4+cc)*cos(2*PI*x2);b[6]=4*PI*PI*(PI/4+cc)*sin(2*PI*x2);b[7]=8*PI*PI*PI*(PI/4+cc)*cos(2*PI*x2);
  
    gauss(a,b,8);
    
    for(int i=0;i<8;i++)
    {
	coef[i]=b[i];
	delete[] a[i];
    }
    delete[] a;
    delete[] b;
    
}

void second_segment(double cc, double *coef)
{
    double **a, *b;
    a=new double* [9]; for(int i=0;i<9;i++) a[i]=new double [9];
    b=new double [9];
    double x1=0.4,x2=0.5,x3=0.6;

    a[0][0]=pow(x1,8);     a[0][1]=pow(x1,7);     a[0][2]=pow(x1,6);     a[0][3]=pow(x1,5);     a[0][4]=pow(x1,4);     a[0][5]=pow(x1,3);     a[0][6]=pow(x1,2);  a[0][7]=x1;a[0][8]=1;
    a[1][0]=8*pow(x1,7);   a[1][1]=7*pow(x1,6);   a[1][2]=6*pow(x1,5);   a[1][3]=5*pow(x1,4);   a[1][4]=4*pow(x1,3);   a[1][5]=3*pow(x1,2);   a[1][6]=2*x1;       a[1][7]=1;a[1][8]=0;
    a[2][0]=56*pow(x1,6);  a[2][1]=42*pow(x1,5);  a[2][2]=30*pow(x1,4);  a[2][3]=20*pow(x1,3);  a[2][4]=12*pow(x1,2);  a[2][5]=6*x1;          a[2][6]=2;          a[2][7]=0;a[2][8]=0;
    a[3][0]=336*pow(x1,5); a[3][1]=210*pow(x1,4); a[3][2]=120*pow(x1,3); a[3][3]=60*pow(x1,2);  a[3][4]=24*x1;         a[3][5]=6;             a[3][6]=0;          a[3][7]=0;a[3][8]=0;
	
    a[4][0]=pow(x2,8);     a[4][1]=pow(x2,7);     a[4][2]=pow(x2,6);     a[4][3]=pow(x2,5);     a[4][4]=pow(x2,4);     a[4][5]=pow(x2,3);     a[4][6]=pow(x2,2);  a[4][7]=x2;a[4][8]=1;
	
    a[5][0]=pow(x3,8);     a[5][1]=pow(x3,7);     a[5][2]=pow(x3,6);     a[5][3]=pow(x3,5);     a[5][4]=pow(x3,4);     a[5][5]=pow(x3,3);     a[5][6]=pow(x3,2);  a[5][7]=x3;a[5][8]=1;
    a[6][0]=8*pow(x3,7);   a[6][1]=7*pow(x3,6);   a[6][2]=6*pow(x3,5);   a[6][3]=5*pow(x3,4);   a[6][4]=4*pow(x3,3);   a[6][5]=3*pow(x3,2);   a[6][6]=2*x3;       a[6][7]=1;a[6][8]=0;
    a[7][0]=56*pow(x3,6);  a[7][1]=42*pow(x3,5);  a[7][2]=30*pow(x3,4);  a[7][3]=20*pow(x3,3);  a[7][4]=12*pow(x3,2);  a[7][5]=6*x3;          a[7][6]=2;          a[7][7]=0;a[7][8]=0;
    a[8][0]=336*pow(x3,5); a[8][1]=210*pow(x3,4); a[8][2]=120*pow(x3,3); a[8][3]=60*pow(x3,2);  a[8][4]=24*x3;         a[8][5]=6;             a[8][6]=0;          a[8][7]=0;a[8][8]=0;
      
    b[0]=PI/2.0-(PI/4+cc)*sin(2*PI*x1);b[1]=-2*PI*(PI/4+cc)*cos(2*PI*x1);b[2]=4*PI*PI*(PI/4+cc)*sin(2*PI*x1);b[3]=8*PI*PI*PI*(PI/4+cc)*cos(2*PI*x1);
    b[4]=PI/2.0;
    b[5]=PI/2.0-(PI/4-cc)*sin(2*PI*x3);b[6]=-2*PI*(PI/4-cc)*cos(2*PI*x3);b[7]=4*PI*PI*(PI/4-cc)*sin(2*PI*x3);b[8]=8*PI*PI*PI*(PI/4-cc)*cos(2*PI*x3);
    
    gauss(a,b,9);
    
    for(int i=0;i<9;i++)
    {
	coef[i]=b[i];
	delete[] a[i];
    }
    delete[] a;
    delete[] b;
}

void third_segment(double cc, double *coef)
{
    double **a,*b;
    a=new double* [5]; for(int i=0;i<5;i++) a[i]=new double [5];
    b=new double [5];
    double x1=0.9,x2=1.0;
      
    a[0][0]=pow(x1,4);     a[0][1]=pow(x1,3);     a[0][2]=pow(x1,2);     a[0][3]=x1;            a[0][4]=1;
    a[1][0]=4*pow(x1,3);   a[1][1]=3*pow(x1,2);   a[1][2]=2*x1;          a[1][3]=1;             a[1][4]=0;
    a[2][0]=12*pow(x1,2);  a[2][1]=6*x1;          a[2][2]=2;             a[2][3]=0;             a[2][4]=0;
    a[3][0]=24*x1;         a[3][1]=6;             a[3][2]=0;             a[3][3]=0;             a[3][4]=0;
	
    a[4][0]=pow(x2,4);     a[4][1]=pow(x2,3);     a[4][2]=pow(x2,2);     a[4][3]=x2;            a[4][4]=0;
  
    b[0]=PI/2.0-(PI/4-cc)*sin(2*PI*x1);b[1]=-2*PI*(PI/4-cc)*cos(2*PI*x1);b[2]=4*PI*PI*(PI/4-cc)*sin(2*PI*x1);b[3]=8*PI*PI*PI*(PI/4-cc)*cos(2*PI*x1);
    b[4]=PI/2.0;
    
    gauss(a,b,5);
	
    for(int i=0;i<5;i++)
    {
	coef[i]=b[i];
	delete[] a[i];
    }
    delete[] a;
    delete[] b;
  
}






