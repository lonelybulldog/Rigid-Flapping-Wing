#include <thrust/copy.h>

#include <cusp/array1d.h>
#include <cusp/print.h>

struct is_even
{
    __host__ __device__
    bool operator()(const int x)
    {
	return (x%2) == 0;
    }
};



int main()
{
    cusp::array1d<int, cusp::device_memory> a(6);
    
    a[0]=-2; a[1]=0; a[2]=-1; a[3]=0; a[4]=1; a[5]=2;
    
    cusp::array1d<int, cusp::device_memory> b(6);
    
    b.erase( thrust::copy_if(a.begin(), a.end(), b.begin(), is_even()) , b.end() )  ;
    
    std::cout<<"Size of b is "<<b.size()<<std::endl;
    
    cusp::print(b);


}