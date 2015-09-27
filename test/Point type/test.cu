#include <thrust/transform.h>
#include <thrust/functional.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/device_vector.h>
#include <iostream>

struct functor : public thrust::unary_function<int,int>
{
    const int IPOINT, JPOINT, KPOINT;
    functor(int _IPOINT, int _JPOINT, int _KPOINT) : IPOINT(_IPOINT), JPOINT(_JPOINT), KPOINT(_KPOINT) {}
    __host__ __device__
    int operator() (int INDEX) const
    {
	const int i = INDEX / (JPOINT * KPOINT);
	const int j = (INDEX % (JPOINT * KPOINT)) / KPOINT;
	const int k = (INDEX % (JPOINT * KPOINT)) % KPOINT;
	
	if ( i==0 || i==(IPOINT-1) || j==0 || j==(JPOINT-1) || k==0 || k==(KPOINT-1) ) return 0;
	else return 1;
    }
};


int main()
{
    int IPOINT=5, JPOINT=4, KPOINT=3;
    thrust::device_vector<double> data(IPOINT*JPOINT*KPOINT);
    
    thrust::transform(thrust::make_counting_iterator(0),thrust::make_counting_iterator(IPOINT*JPOINT*KPOINT),data.begin(),functor(IPOINT,JPOINT,KPOINT));
    
    for(int i = 0; i < data.size(); i++)
        std::cout << "data[" << i << "] = " << data[i] << std::endl;
}