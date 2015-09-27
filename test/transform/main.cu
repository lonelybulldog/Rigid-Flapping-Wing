#include <iostream>
#include <fstream>

#include <cusp/print.h>


int main()
{
    cusp::array1d<int,cusp::device_memory> POINTTYPE;

    POINTTYPE.resize(10);
    
    thrust::negate<int> op;


    thrust::transform(thrust::make_counting_iterator(1), 
                      thrust::make_counting_iterator(10),
                      POINTTYPE.begin(),
                      op);

    
    cusp::print(POINTTYPE);
    
    return 0;
}