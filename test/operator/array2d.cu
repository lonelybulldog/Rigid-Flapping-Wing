#include <cusp/array1d.h>
#include <cusp/print.h>
#include <thrust/sequence.h>
#include <iostream>


class MESH_CARTESIAN
{
public:
    int POINT_CARTESIAN;
    cusp::array1d<double, cusp::device_memory> a;
    
    MESH_CARTESIAN(int s)
    {
	POINT_CARTESIAN = s;
	a.resize(POINT_CARTESIAN);
	thrust::sequence(a.begin(),a.end(),1,1);
    }
    friend class MESH;
};

class MESH_LESS
{
public:
    int POINT_MESHLESS;
    cusp::array1d<double, cusp::device_memory> b;
    
    MESH_LESS(int s)
    {
	POINT_MESHLESS = s;
	b.resize(POINT_MESHLESS);
	thrust::sequence(b.begin(),b.end(),10,10);
    }
    friend class MESH;
};




class MESH
{
public:
    int POINT_CARTESIAN, POINT_MESHLESS, POINT_ALL;
   
    class X
    {
    public:
//	cusp::array1d<double, cusp::device_memory> a;
//	cusp::array1d<double, cusp::device_memory> b;
	double *a, *b;
	int POINT_CARTESIAN;
	
	X()
	{
	   std::cout << "X is initilized." << std::endl;
	}
	
	~X()
	{
	    std::cout << "X is deleted." << std::endl;
	}
	
	__host__ __device__
	double operator[](int i)
	{
	    if(i<POINT_CARTESIAN) return a[i];
	    else return b[i-POINT_CARTESIAN];
	}
    };
    
    X _X;
    
    MESH(MESH_CARTESIAN& cartesian, MESH_LESS& meshless)
    {
	std::cout << "MESH is initilized." << std::endl;
	POINT_CARTESIAN = cartesian.POINT_CARTESIAN;
	POINT_MESHLESS = meshless.POINT_MESHLESS;
	POINT_ALL = POINT_CARTESIAN + POINT_MESHLESS;
	
	_X.POINT_CARTESIAN = POINT_CARTESIAN;
//	_X.a = cartesian.a;
//	_X.b = meshless.b;
	
	_X.a = thrust::raw_pointer_cast( cartesian.a.data() );
	_X.b = thrust::raw_pointer_cast( meshless.b.data() );
	
//	_X.a = new double [POINT_CARTESIAN];
//	_X.b = new double [POINT_MESHLESS];
	
//	thrust::copy(cartesian.a.begin(), cartesian.a.end(), _X.a);
//	thrust::copy(meshless.b.begin(),meshless.b.end(), _X.b);
	
	
	std::cout << "MESH initilization is done." << std::endl;
    }
    
    ~MESH()
    {
	delete[] _X.a;
	delete[] _X.b;
	std::cout << "MESH is deleted." << std::endl;
    }
    
};









int main()
{
//    MESH mesh(5,5);
//    mesh.-X(5);
  
    MESH_CARTESIAN cartesian(5);
    MESH_LESS      meshless(5);
    MESH           mesh(cartesian,meshless);
  
//    cusp::print(mesh._X.a);
//    cusp::print(mesh._X.b);
    
    std::cout<< "mark" << std::endl;
    
    cusp::array1d<double,cusp::device_memory> result(10);
    
//    std::cout << mesh._X[0] << std::endl;
//    for(int i=0; i<10;i++) std::cout<< mesh._X[i] << std::endl;
    thrust::copy(&mesh._X.a[0],&mesh._X.a[0]+5,result.begin());
    thrust::copy(&mesh._X.b[0],&mesh._X.b[0]+5,result.begin()+5);
    
    cusp::print(result);
    
    std::cout<< "end" << std::endl;
    
    thrust::sequence(cartesian.a.begin(),cartesian.a.end(),100,100);
    
    
    thrust::copy(&mesh._X.a[0],&mesh._X.a[0]+5,result.begin());
    thrust::copy(&mesh._X.b[0],&mesh._X.b[0]+5,result.begin()+5);
    for(int i=0; i<10;i++) std::cout<< mesh._X[i] << std::endl;
//    cusp::print(result);
    
//    std::cout<<result[0]<<" "<<result[1]<<" "<<std::endl;
    
    return 0;
}
