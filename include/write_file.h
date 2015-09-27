#ifndef WRITE_FILE_H
#define WRITE_FILE_H

#include "mesh.h"
#include "frame_structure.h"
#include "flow_field.h"
#include <cusp/array1d.h>
#include <cusp/array2d.h>

#define LOCAL       1
#define CLUSTER     2


class WRITE_FILE
{
public:
    int FREQ_CARTESIAN, FREQ_MESHLESS, FREQ_REALSURFACE;
    
    int N;
    
    REF_FRAME     *global;
    MESH          *mesh;
    FLOW_FIELD    *flow_field;
    SOLVER        *solver;
    FLAPPATTERN   *flappattern;
    
    WRITE_FILE(REF_FRAME& global, MESH& mesh, FLOW_FIELD& flow_field, SOLVER& solver, FLAPPATTERN& flappattern);
    ~WRITE_FILE();
    
    void   OUTPUT_RESULT(int it);

private:
    cusp::array1d<int,cusp::host_memory> IINDEX;
    cusp::array1d<int,cusp::host_memory> JINDEX;
    cusp::array1d<double,cusp::host_memory> XYZ[3];
  
    cusp::array1d<double,cusp::host_memory> X,Y,Z;
    
    cusp::array2d<int,cusp::host_memory,cusp::column_major> *TETRAHEDRON;
    cusp::array2d<int,cusp::host_memory,cusp::column_major> *TRIANGLE;
    
    void   OUTPUT_PLT(int it);
    void   OUTPUT_PLT_CARTESIAN(int it);
    void   OUTPUT_PLT_MESHLESS(int it);
    void   OUTPUT_PLT_REALSURFACE(int it);
    
    void   OUTPUT_DAT(int it);
    void   OUTPUT_DAT_CARTESIAN(int it);
    void   OUTPUT_DAT_MESHLESS(int it);
    void   OUTPUT_DAT_REALSURFACE(int it);
};

#endif /*WRITE_FILE_H*/