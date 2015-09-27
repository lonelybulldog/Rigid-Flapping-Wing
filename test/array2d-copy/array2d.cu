#include <cusp/array1d.h>
#include <cusp/array2d.h>
#include <cusp/print.h>
#include <iostream>

int main()
{
    cusp::array2d<int, cusp::device_memory, cusp::column_major> a(3,3);
    cusp::array1d<int, cusp::device_memory> b;
    
    a(0,0) = 0; a(0,1) = 1; a(0,2) = 2;
    a(1,0) = 3; a(1,1) = 4; a(1,2) = 5;
    a(2,0) = 6; a(2,1) = 7; a(2,2) = 8;
    
    b.resize(9);
    

    for(int i=0;i<3;i++)
    {
	thrust::copy( a.row(i).begin(), a.row(i).end(), b.begin()+3*i);
    } 
    

    cusp::print(a);
    
    cusp::print(b);
    
    return 0;
}
