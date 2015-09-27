


template<typename MemorySpace> class MESH_CARTESIAN
{
public:
    int IPOINT, JPOINT, KPOINT, POINT_CARTESIAN;
    
    cusp::array1d<double,MemorySpace> XYZ[3];
    cusp::array1d<double,MemorySpace> DELTA[3];
    double VELOCITY[3],ACCELERATION[3];

    cusp::array1d<int,MemorySpace> IINDEX;
    cusp::array1d<int,MemorySpace> JINDEX;
    cusp::array1d<int,MemorySpace> POINTTYPE[3];
    
    cusp::array1d<int,MemorySpace> COUNTTYPE;  /*thurst::exclusive_scan*/
    
    MESH_CARTESIAN();
    ~MESH_CARTESIAN();
    
    
private:
    void Div(double start, int segment, double *segmentlength, double *ratio, int *segmentpoint, int POINT, double *temp);
};

template<typename MemorySpace> class MESH_LESS
{
public:
    int POINT_MESSLESS;
    
    cusp::array1d<double,MemorySpace> POSITION[3];
    cusp::array1d<double,MemorySpace> VELOCITY[3];
    cusp::array1d<double,MemorySpace> ACCELERATION[3];

    MESH_LESS(int Point_Messless);
    ~MESH_LESS();
    
};
