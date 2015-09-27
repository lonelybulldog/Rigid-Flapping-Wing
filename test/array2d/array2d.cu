#include <cusp/array1d.h>
#include <cusp/array2d.h>
#include <cusp/print.h>
#include <iostream>

int main()
{
    cusp::array2d<int, cusp::device_memory> a(3,3);
    
    a(0,0) = 0; a(0,1) = 1; a(0,2) = 2;
    a(1,0) = 3; a(1,1) = 4; a(1,2) = 5;
    a(2,0) = 6; a(2,1) = 7; a(2,2) = 8;
    
    cusp::array1d<int, cusp::device_memory> b(3);
    
    b[0] = 10; b[1] = 20; b[2] = 30;
    
    int *raw;

    cusp::print(a.values);

//    raw = thrust::raw_pointer_cast( &(a.row(1)[0]) );
//    raw = thrust::raw_pointer_cast( b.data() );
    
//    std::cout<< raw[0] << " "<< raw[1] << " " << raw[2] << std::endl;
//    std::cout<< a.row(0)[0] << " "<< a.row(0)[1] << " " << a.row(0)[2] << std::endl;
    raw = new int[3];
    
    thrust::copy(a.row(0).begin(),a.row(0).end(),raw);
    std::cout<< raw[0] << " "<< raw[1] << " " << raw[2] << std::endl;
    
    

    
//    thrust::transform(a.row(0).begin(), a.row(0).end(), thrust::make_constant_iterator(-100), a.row(0).begin(), thrust::plus<double>() );
//    thrust::transform(a.row(1).begin(), a.row(1).end(), thrust::make_constant_iterator(-200), a.row(1).begin(), thrust::plus<double>() );
//    thrust::transform(a.row(2).begin(), a.row(2).end(), thrust::make_constant_iterator(-300), a.row(2).begin(), thrust::plus<double>() );
    
    cusp::print(a);
    
    
    
    return 0;
}
