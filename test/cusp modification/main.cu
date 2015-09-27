#include <cusp/array2d.h>
#include <cusp/blas/blas.h>
#include <cusp/print.h>

#include <iostream>


int main()
{
    cusp::array2d<double, cusp::device_memory, cusp::column_major> A(2,2),B(2,3),C(2,3);
    
    A(0,0) = 1; A(0,1) = 2;
    A(1,0) = 3; A(1,1) = 4;
    
    B(0,0) = 5; B(0,1) = 6; B(0,2) = 7;
    B(1,0) = 8; B(1,1) = 9; B(1,2) = 10;

    cusp::blas::gemm(A,B,C);
    
//    cusp::blas::geam(A,B,A,-1.0,100.0);
    
    cusp::print(A);
    
    cusp::print(B);
    
    cusp::print(C);
    
//    C=A;
    
//    cusp::print(C);
    
    
    return 0;
}