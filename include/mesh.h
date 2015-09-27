#ifndef MESH_H
#define MESH_H

#include <stdio.h>
#include <cusp/array1d.h>

#include "frame_structure.h"

class LINKLIST;
class FLOW_FIELD;

class MESH_CARTESIAN
{
public:
    int IPOINT, JPOINT, KPOINT, POINT_CARTESIAN;
    double XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, MESHSIZE, SAFEDISTANCE;
    
    cusp::array1d<double,cusp::device_memory> XYZ[3];
    cusp::array1d<double,cusp::device_memory> DELTA[3];
    cusp::array1d<double,cusp::device_memory> VELOCITY;
    cusp::array1d<double,cusp::device_memory> ACCELERATION;

    cusp::array1d<int,cusp::device_memory> IINDEX;
    cusp::array1d<int,cusp::device_memory> JINDEX;
    cusp::array1d<int,cusp::device_memory> POINTTYPE[3]; /* 0 : previous step, 1 : current step, 2 : temporary step */
    
    cusp::array1d<int,cusp::device_memory> COUNTTYPE;  /*thrust::exclusive_scan*/
    
    MESH_CARTESIAN();
    ~MESH_CARTESIAN();
    
    friend class LINKLIST;
    friend class MESH;
    
private:
    void Div(double start, int segment, double *segmentlength, double *ratio, int *segmentpoint, int POINT, double *temp);
};

class MESH_LESS
{
public:
    int POINT_MESSLESS, INNER_POINT, OUTER_POINT, *INNER_POINT_OFFSET, *OUTER_POINT_OFFSET, *OFFSET;
    
    cusp::array1d<double,cusp::device_memory> POSITION[3];
    cusp::array1d<double,cusp::device_memory> VELOCITY[3];
    cusp::array1d<double,cusp::device_memory> ACCELERATION[3];
    cusp::array2d<double,cusp::device_memory,cusp::column_major> OUTER_NORMAL_VECTOR;
    cusp::array1d<int,cusp::device_memory> POINTTYPE[3];

    MESH_LESS(REF_FRAME& global);
    ~MESH_LESS();
    
    cusp::array1d<int,cusp::device_memory> INNER_NODE_INDEX;
    cusp::array1d<int,cusp::device_memory> OUTER_NODE_INDEX;
    cusp::array1d<double,cusp::device_memory> SURFACE_ELE_AREA;
    
    friend class MESH;
};

class MESH
{
public:   
    MESH_CARTESIAN   *cartesian;
    MESH_LESS        *meshless;
    int IPOINT, JPOINT, KPOINT, POINT_CARTESIAN, POINT_MESSLESS, POINT_ALL;
    
    class UNIFIED_POSITION
    {
    public:
        int IPOINT, JPOINT, KPOINT, POINT_CARTESIAN, POINT_MESSLESS, POINT_ALL;
        
        double *XYZ_raw[3], *POSITION_MESHLESS_raw[3];
        
        int ijk[3];
        
	__host__ __device__
        double operator()( int s, int INDEX )
        {
            if( 0 <= INDEX && INDEX < POINT_CARTESIAN )
            {
                ijk[0] = INDEX / (JPOINT * KPOINT);
                ijk[1] = (INDEX % (JPOINT * KPOINT)) / KPOINT;
                ijk[2] = (INDEX % (JPOINT * KPOINT)) % KPOINT;
                
                return XYZ_raw[s][ijk[s]];
            }
            else if ( POINT_CARTESIAN <= INDEX && INDEX < POINT_ALL )
            {
                return POSITION_MESHLESS_raw[s][INDEX-POINT_CARTESIAN];
            }
            else
            {
		printf("Index beyond the scope in mesh operator.\n");
//                std::cout << "Index beyond the scope in mesh operator." << std::endl;
                return 0;
            }
        }
    } XYZ;
    
    class UNIFIED_VELOCITY
    {
    public:
        int IPOINT, JPOINT, KPOINT, POINT_CARTESIAN, POINT_MESSLESS, POINT_ALL;
        
        double *VELOCITY_CARTESIAN_raw, *VELOCITY_MESHLESS_raw[3];
        
	__host__ __device__
        double operator()( int s, int INDEX )
        {
            if( 0 <= INDEX && INDEX < POINT_CARTESIAN )
            {
                return VELOCITY_CARTESIAN_raw[s];
            }
            else if ( POINT_CARTESIAN <= INDEX && INDEX < POINT_ALL )
            {
                return VELOCITY_MESHLESS_raw[s][INDEX-POINT_CARTESIAN];
            }
            else
            {
		printf("Index beyond the scope in mesh operator.\n");
//                std::cout << "Index beyond the scope in mesh operator." << std::endl;
                return 0;
            }
        }
    } UVW;
    
    class UNIFIED_ACCELERATION
    {
    public:
        int IPOINT, JPOINT, KPOINT, POINT_CARTESIAN, POINT_MESSLESS, POINT_ALL;
        
        double *ACCELERATION_CARTESIAN_raw, *ACCELERATION_MESHLESS_raw[3];
        
	__host__ __device__
        double operator()( int s, int INDEX)
        {
            if( 0 <= INDEX && INDEX < POINT_CARTESIAN )
            {
                return ACCELERATION_CARTESIAN_raw[s];
            }
            else if ( POINT_CARTESIAN <= INDEX && INDEX < POINT_ALL )
            {
                return ACCELERATION_MESHLESS_raw[s][INDEX-POINT_CARTESIAN];
            }
            else
            {
		printf("Index beyond the scope in mesh operator.\n");
//                std::cout << "Index beyond the scope in mesh operator." << std::endl;
                return 0;
            }

        }
    } ACC;
    
    MESH( MESH_CARTESIAN& Cartesian, MESH_LESS& Meshless );
    ~MESH();
    void UPDATE_MESHLESS(REF_FRAME& global);
    void SEARCH_TYPE(bool Initialization, LINKLIST& linklist, FLOW_FIELD& flow_field);
    void UPDATE_POINTTYPE();
    void UPDATE_POINTTYPE_IMPLICIT();
    void TEST_FUNCTION();
    
    friend class FLOW_FIELD;
    friend class SOLVER;
    
};

#endif /*MESH_H*/