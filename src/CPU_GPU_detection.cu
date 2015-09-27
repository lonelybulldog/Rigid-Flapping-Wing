#include <stdio.h>
#include <iostream>
#include <omp.h>

bool GPU_detection()
{
    int deviceCount = 0;
    cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

    if (error_id != cudaSuccess)
    {
        printf("cudaGetDeviceCount returned %d\n-> %s\n", (int)error_id, cudaGetErrorString(error_id));
        printf("Result = FAIL\n");
        exit(EXIT_FAILURE);
    }

    if (deviceCount == 0)
    {
        printf("There are no available device(s) that support CUDA\n");
	return false;
    }
    else
    {
        printf("Detected %d CUDA Capable device(s)\n", deviceCount);
    }
    
    int dev;
    for (dev = 0; dev < deviceCount; ++dev)
    {
        cudaSetDevice(dev);
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, dev);

        printf("Device %d: \"%s\"  |  ", dev, deviceProp.name);
    }
    std::cout << std::endl << std::endl;
    
    //manually specfiy the GPU device.
    dev = 0;
    cudaSetDevice(dev);
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, dev);
    std::cout << "Set CUDA device(s) to " << deviceProp.name << "." << std::endl << std::endl;

    return true;
}

void CPU_detection()
{
    int No_CPU = omp_get_num_procs();
    
    omp_set_num_threads(No_CPU);
    
    std::cout << "Set OpenMP thread(s) to " << No_CPU << "." << std::endl << std::endl;
}