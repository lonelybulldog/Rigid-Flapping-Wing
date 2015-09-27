#include <thrust/system/omp/vector.h>
#include <thrust/system/tbb/vector.h>
#include <thrust/iterator/retag.h>
#include <cstdio>

struct omp_hello
{
  void operator()(int x)
  {
    printf("Hello, world from OpenMP!\n");
  }
};

struct tbb_hello
{
  void operator()(int x)
  {
    printf("Hello, world from TBB!\n");
  }
};

int main()
{
  thrust::omp::vector<int> omp_vec(2, 7);
  thrust::tbb::vector<int> tbb_vec(2, 13);

  thrust::for_each(thrust::reinterpret_tag<thrust::tbb::tag>(omp_vec.begin()),
                   thrust::reinterpret_tag<thrust::tbb::tag>(omp_vec.end()), 
                   tbb_hello());

  thrust::for_each(thrust::reinterpret_tag<thrust::omp::tag>(tbb_vec.begin()),
                   thrust::reinterpret_tag<thrust::omp::tag>(tbb_vec.end()),
                   omp_hello());
}
