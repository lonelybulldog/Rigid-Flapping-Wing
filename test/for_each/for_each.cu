#include <thrust/for_each.h>
#include <thrust/device_vector.h>
#include <cusp/array1d.h>
#include <cusp/print.h>
#include <iostream>

#define NB 40

typedef struct LINKLIST_MEMBER
{
    int Meshless_Ind;
    int Nb_Points[NB];
    int Nb_Points_Old[NB];
    double Csvd[9][NB];
    double Csvd_Old[9][NB];
}LINKLIST_MEMBER;



struct functor 
{
    __host__ __device__
    void operator()(LINKLIST_MEMBER& x)
    {
	int i = x.Meshless_Ind;
	x.Meshless_Ind = -i;
    }
};
  
int main()
{

    cusp::array1d<LINKLIST_MEMBER,cusp::host_memory> LINKLIST_h(3);

    LINKLIST_h[0].Meshless_Ind = 10;
    LINKLIST_h[1].Meshless_Ind = 20;
    LINKLIST_h[2].Meshless_Ind = 30;
    
    cusp::array1d<LINKLIST_MEMBER,cusp::device_memory> LINKLIST = LINKLIST_h;

    thrust::for_each(LINKLIST.begin(), LINKLIST.end(), functor());
  
    std::cout<<" it contains: "<< std::endl;
    for(thrust::device_vector<LINKLIST_MEMBER>::iterator iter = LINKLIST.begin(); iter != LINKLIST.end(); iter++)  
    {
	std::cout << (static_cast<LINKLIST_MEMBER>(*iter)).Meshless_Ind << std::endl;
    }
   
  return 0;
}