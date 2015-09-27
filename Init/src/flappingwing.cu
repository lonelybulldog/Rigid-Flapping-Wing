#include <iostream>
#include <fstream>
#include <omp.h>

#include <cusp/print.h>

#include "../include/flow_field.h"

#include "../include/timer.h"

#include "./mesh.hpp"


bool GPU_detection();
void CPU_detection();



int main(int argc, char *argv[])
{
    std::cout << "Task Starting......" << std::endl << std::endl;
    
    GPU_detection();
    CPU_detection();
    
//    #pragma omp parallel
//    {
//	std::cout<<" Threads id: "<< omp_get_thread_num() << std::endl;
//    }
    
    
/******************************************************/
    typedef cusp::host_memory MemorySpace;
/******************************************************/
    
    timer time;

    FLOW_FIELD<MemorySpace> flow_field(2,2);
    
    flow_field.U(0,0)=0; flow_field.U(0,1)=1; flow_field.U(0,2)=2; flow_field.U(0,3)=3;
    flow_field.U(1,0)=4; flow_field.U(1,1)=5; flow_field.U(1,2)=6; flow_field.U(1,3)=7;
    
    cusp::print(flow_field.U);

    MESH_CARTESIAN<MemorySpace> cartesian;
//    cusp::print(cartesian.DELTA[1]);
    
//    std::cout<< cartesian.POINTTYPE[0].size() << " " << cartesian.POINTTYPE[1].size() << " " << cartesian.POINTTYPE[2].size() << std::endl;
    
    std::cout << "Elapsed time: " << time.milliseconds_elapsed() << "ms" << std::endl;

}
