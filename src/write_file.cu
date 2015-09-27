#include <iostream>
#include <fstream>
#include <string.h>
#include <omp.h>

#if MACHINE_SYSTEM == LOCAL
#include "/usr/local/tec360/include/MASTER.h" /*--Tecplot Head--*/
#include "/usr/local/tec360/include/TECIO.h"  /*--Tecplot Head--*/
#endif

#include "../include/write_file.h"

#define  N_VARS 7


WRITE_FILE :: WRITE_FILE(REF_FRAME& Global, MESH& Mesh, FLOW_FIELD& Flow_field, SOLVER& Solver, FLAPPATTERN& Flappattern)
{
    std::cout << "Write_file is being initialized......";
    
    global      = &Global;
    mesh        = &Mesh;
    flow_field  = &Flow_field;
    solver      = &Solver;
    flappattern = &Flappattern;
    
    FREQ_CARTESIAN   = 50;
    FREQ_MESHLESS    = 50;
    FREQ_REALSURFACE = 50;
    
    IINDEX = mesh->cartesian->IINDEX;
    JINDEX = mesh->cartesian->JINDEX;
    for(int i=0;i<3;i++) XYZ[i] = mesh->cartesian->XYZ[i];
    
    X.resize(mesh->POINT_CARTESIAN);
    Y.resize(mesh->POINT_CARTESIAN);
    Z.resize(mesh->POINT_CARTESIAN);
    
    int ijk;
    #pragma omp parallel for private(j,k,ijk)
    for(int i = 0; i < mesh->IPOINT; i++)
    {
	for(int j = 0; j < mesh->JPOINT; j++)
	{
	    for(int k = 0; k < mesh->KPOINT; k++)
	    {
		ijk=IINDEX[i]+JINDEX[j]+k;
		X[ijk]=XYZ[0][i];
		Y[ijk]=XYZ[1][j];
		Z[ijk]=XYZ[2][k];
	    }
	}
    }
    
    N = global->obj_number;
    TETRAHEDRON = new cusp::array2d<int,cusp::host_memory,cusp::column_major> [N];
    TRIANGLE    = new cusp::array2d<int,cusp::host_memory,cusp::column_major> [N];
    
    for(int i=0;i<N;i++)
    {
	TETRAHEDRON[i]=global->rigid_body[i]->TETRAHEDRON;
	TRIANGLE[i]=global->rigid_body[i]->TRIANGLE;
    }
    
    std::cout << "Initialization done." << std::endl;
}

WRITE_FILE :: ~WRITE_FILE()
{
    delete [] TETRAHEDRON;
    delete [] TRIANGLE;
    std::cout << "Write_file is being deleted." << std::endl;
}

void WRITE_FILE :: OUTPUT_RESULT(int it)
{
    if( it%FREQ_CARTESIAN == 0 || it%FREQ_CARTESIAN == 0 || it%FREQ_REALSURFACE == 0 )
    {
#if   MACHINE_SYSTEM == CLUSTER
	OUTPUT_DAT(it);
#elif MACHINE_SYSTEM == LOCAL
	OUTPUT_PLT(it);
#else
	std::cout << "Undefined MACHINE_SYSTEM. Output failure." << std::endl;
#endif
    }
}

void WRITE_FILE :: OUTPUT_DAT(int it)
{
    if( it%FREQ_CARTESIAN == 0 )   OUTPUT_DAT_CARTESIAN(it);
    if( it%FREQ_MESHLESS == 0 )    OUTPUT_DAT_MESHLESS(it);
    if( it%FREQ_REALSURFACE == 0 ) OUTPUT_DAT_REALSURFACE(it);
}

void WRITE_FILE :: OUTPUT_DAT_CARTESIAN(int iter)
{
    char *filename = new char [30];	
    char Char_Buffer[30];
    strcpy(filename, "../result/Cartesian-");
    sprintf(Char_Buffer, "%d", int(iter/FREQ_CARTESIAN));
    strcat(filename, Char_Buffer);
    strcat(filename, ".dat");

    std::ofstream out;
    out.open(filename);
    out<<"TITLE= \"Cartesian Point\""<<std::endl;
    out<<"VARIABLES = \"X\", \"Y\", \"Z\", \"U\", \"V\", \"W\", \"P\", \"POINTTYPE\""<<std::endl;
    out<<"Zone T="<<"\""<<4<<"\""<<", "<<"SOLUTIONTIME="<<global->Time<<", "
		  <<", I="<<mesh->IPOINT<<", J="<<mesh->JPOINT<<", K="<<mesh->KPOINT<<", F=POINT"<<std::endl;
    int ijk;
    cusp::array1d<double,cusp::host_memory> TEMPU(mesh->POINT_CARTESIAN), TEMPV(mesh->POINT_CARTESIAN), TEMPW(mesh->POINT_CARTESIAN), TEMPP(mesh->POINT_CARTESIAN);
    thrust::copy(flow_field->U.row(0).begin(), flow_field->U.row(0).begin() + mesh->POINT_CARTESIAN, TEMPU.begin());
    thrust::copy(flow_field->V.row(0).begin(), flow_field->V.row(0).begin() + mesh->POINT_CARTESIAN, TEMPV.begin());
    thrust::copy(flow_field->W.row(0).begin(), flow_field->W.row(0).begin() + mesh->POINT_CARTESIAN, TEMPW.begin());
    thrust::copy(flow_field->P.row(0).begin(), flow_field->P.row(0).begin() + mesh->POINT_CARTESIAN, TEMPP.begin());
    cusp::array1d<int,cusp::host_memory> TEMPT = mesh->cartesian->POINTTYPE[1];
    
    for(int k=0;k<mesh->KPOINT;k++) {
	for(int j=0;j<mesh->JPOINT;j++) {
	    for(int i=0;i<mesh->IPOINT;i++) {
		ijk=IINDEX[i]+JINDEX[j]+k;
//		if(i==108&&j==85&&k==56) std::cout<<mesh->cartesian->POINTTYPE[1][ijk]<<std::endl;
		out<<XYZ[0][i]<<" "<<XYZ[1][j]<<" "<<XYZ[2][k]<<" "
		<<TEMPU[ijk]<<" "<<TEMPV[ijk]<<" "<<TEMPW[ijk]<<" "<<TEMPP[ijk]<<" "<<TEMPT[ijk]<<std::endl;
//		<<TEMPU[ijk]<<" "<<TEMPV[ijk]<<" "<<TEMPW[ijk]<<" "<<TEMPP[ijk]<<" "<<mesh->cartesian->POINTTYPE[1][ijk]<<std::endl;
	    }
	}
    }
    out.close();
    delete filename;
}

void WRITE_FILE :: OUTPUT_DAT_MESHLESS(int iter)
{
    char *filename = new char [30];
    char Char_Buffer[30];
    strcpy(filename, "../result/Meshless-");
    sprintf(Char_Buffer, "%d", int(iter/FREQ_MESHLESS));
    strcat(filename, Char_Buffer);
    strcat(filename, ".dat");

    int NumNodes;

    std::ofstream out;
    out.open(filename);
    out<<"TITLE= \"Meshless Point\""<<std::endl;
    out<<"VARIABLES = \"X\", \"Y\", \"Z\", \"U\", \"V\", \"W\", \"P\", \"POINTTYPE\", \"ONVx\", \"ONVy\", \"ONVz\""<<std::endl;
    for(int i = 0; i < global->obj_number; i++)
    {
	out<<"ZONE T="<<"\""<<i<<"\""<<", "<<"SOLUTIONTIME="<<global->Time<<", "
		<<"N="<<global->rigid_body[i]->POINT_NUMBER<<", "
		<<"E="<<global->rigid_body[i]->TETRAHEDRON_NUMBER<<", "
		<<"DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON"<<std::endl;

	NumNodes = global->rigid_body[i]->POINT_NUMBER;
	
	cusp::array1d<double,cusp::host_memory> TEMPX(NumNodes),TEMPY(NumNodes),TEMPZ(NumNodes),TEMPU(NumNodes),TEMPV(NumNodes),TEMPW(NumNodes),TEMPP(NumNodes);
	cusp::array1d<double,cusp::host_memory> TEMPONVx(NumNodes), TEMPONVy(NumNodes), TEMPONVz(NumNodes);
	cusp::array1d<int,cusp::host_memory> TEMPT(NumNodes);
	
	TEMPX = global->rigid_body[i]->XYZ.row(0);
	TEMPY = global->rigid_body[i]->XYZ.row(1);
	TEMPZ = global->rigid_body[i]->XYZ.row(2);
	thrust::copy(flow_field->U.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i], 
		     flow_field->U.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i] + NumNodes,
		     TEMPU.begin());
	thrust::copy(flow_field->V.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i], 
		     flow_field->V.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i] + NumNodes,
		     TEMPV.begin());
	thrust::copy(flow_field->W.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i], 
		     flow_field->W.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i] + NumNodes,
		     TEMPW.begin());
	thrust::copy(flow_field->P.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i], 
		     flow_field->P.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i] + NumNodes,
		     TEMPP.begin());
	thrust::copy(mesh->meshless->POINTTYPE[1].begin() + mesh->meshless->OFFSET[i], 
		     mesh->meshless->POINTTYPE[1].begin() + mesh->meshless->OFFSET[i] + NumNodes, 
		     TEMPT.begin());
	thrust::copy(mesh->meshless->OUTER_NORMAL_VECTOR.row(0).begin() + mesh->meshless->OFFSET[i],
		     mesh->meshless->OUTER_NORMAL_VECTOR.row(0).begin() + mesh->meshless->OFFSET[i] + NumNodes, 
		     TEMPONVx.begin());
	thrust::copy(mesh->meshless->OUTER_NORMAL_VECTOR.row(1).begin() + mesh->meshless->OFFSET[i],
		     mesh->meshless->OUTER_NORMAL_VECTOR.row(1).begin() + mesh->meshless->OFFSET[i] + NumNodes, 
		     TEMPONVy.begin());
	thrust::copy(mesh->meshless->OUTER_NORMAL_VECTOR.row(2).begin() + mesh->meshless->OFFSET[i],
		     mesh->meshless->OUTER_NORMAL_VECTOR.row(2).begin() + mesh->meshless->OFFSET[i] + NumNodes, 
		     TEMPONVz.begin());
	
	for(int j=0;j<global->rigid_body[i]->POINT_NUMBER;j++) 
	{
	    out<<TEMPX[j]<<" "<<TEMPY[j]<<" "<<TEMPZ[j]<<" "<<TEMPU[j]<<" "<<TEMPV[j]<<" "<<TEMPW[j]<<" "<<TEMPP[j]<<" "<<TEMPT[j]<<" ";
//	    out<<mesh->meshless->POSITION[0][j+mesh->meshless->OFFSET[i]]<<" ";
//	    out<<mesh->meshless->POSITION[1][j+mesh->meshless->OFFSET[i]]<<" ";
//	    out<<mesh->meshless->POSITION[2][j+mesh->meshless->OFFSET[i]]<<" ";
//	    out<<TEMPU[j]<<" "<<TEMPV[j]<<" "<<TEMPW[j]<<" "<<TEMPP[j]<<" "<<TEMPT[j]<<" ";
	    out<<TEMPONVx[j]<<" "<<TEMPONVy[j]<<" "<<TEMPONVz[j]<<std::endl;  
	}

	for(int j=0;j<global->rigid_body[i]->TETRAHEDRON_NUMBER;j++)
	{
	    out<<TETRAHEDRON[i](0,j)<<" "<<TETRAHEDRON[i](1,j)<<" "<<TETRAHEDRON[i](2,j)<<" "<<TETRAHEDRON[i](3,j)<<std::endl;
	}
    }
    out.close();
    delete filename;
}

void WRITE_FILE :: OUTPUT_DAT_REALSURFACE(int iter)
{
    char *filename = new char [30];
    char Char_Buffer[30];
    strcpy(filename, "../result/RealSurface-");
    sprintf(Char_Buffer, "%d", int(iter/FREQ_REALSURFACE));
    strcat(filename, Char_Buffer);
    strcat(filename, ".dat");
    
    int NumNodes;
    
    std::ofstream out; 
    out.open(filename);
    out<<"TITLE= \"Surface Point\""<<std::endl;
    out<<"VARIABLES = \"X\", \"Y\", \"Z\", \"U\", \"V\", \"W\", \"P\", \"POINTTYPE\""<<std::endl;
    for(int i = 0; i < global->obj_number; i++) 
    {
	out<<"ZONE T="<<"\""<<i<<"\""<<", "<<"SOLUTIONTIME="<<global->Time<<", "
		<<"N="<<global->rigid_body[i]->POINT_NUMBER<<", "
		<<"E="<<global->rigid_body[i]->TRIANGLE_NUMBER<<", "
		<<"DATAPACKING=POINT, ZONETYPE=FETRIANGLE"<<std::endl;
	
	NumNodes = global->rigid_body[i]->POINT_NUMBER;
	
	cusp::array1d<double,cusp::host_memory> TEMPX(NumNodes),TEMPY(NumNodes),TEMPZ(NumNodes),TEMPU(NumNodes),TEMPV(NumNodes),TEMPW(NumNodes),TEMPP(NumNodes);
	cusp::array1d<int,cusp::host_memory> TEMPT(NumNodes);
	
	TEMPX = global->rigid_body[i]->XYZ.row(0);
	TEMPY = global->rigid_body[i]->XYZ.row(1);
	TEMPZ = global->rigid_body[i]->XYZ.row(2);
	thrust::copy(flow_field->U.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i], 
		     flow_field->U.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i] + NumNodes,
		     TEMPU.begin());
	thrust::copy(flow_field->V.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i], 
		     flow_field->V.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i] + NumNodes,
		     TEMPV.begin());
	thrust::copy(flow_field->W.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i], 
		     flow_field->W.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i] + NumNodes,
		     TEMPW.begin());
	thrust::copy(flow_field->P.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i], 
		     flow_field->P.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i] + NumNodes,
		     TEMPP.begin());
	thrust::copy(mesh->meshless->POINTTYPE[1].begin() + mesh->meshless->OFFSET[i], 
		     mesh->meshless->POINTTYPE[1].begin() + mesh->meshless->OFFSET[i] + NumNodes, 
		     TEMPT.begin());
	
	for(int j=0;j<global->rigid_body[i]->POINT_NUMBER;j++) 
	{
	    out<<TEMPX[j]-flappattern->bodycentre[0]+solver->position[0]<<" "
	       <<TEMPY[j]-flappattern->bodycentre[1]+solver->position[1]<<" "
	       <<TEMPZ[j]-flappattern->bodycentre[2]+solver->position[2]<<" "
	       <<TEMPU[j]<<" "<<TEMPV[j]<<" "<<TEMPW[j]<<" "<<TEMPP[j]<<" "<<TEMPT[j]<<std::endl;   
	}

	for(int j=0;j<global->rigid_body[i]->TRIANGLE_NUMBER;j++)
	{
	    out<<TRIANGLE[i](0,j)<<" "<<TRIANGLE[i](1,j)<<" "<<TRIANGLE[i](2,j)<<std::endl;
	}
    }
    out.close();
    delete filename;
}


#if MACHINE_SYSTEM == LOCAL
void WRITE_FILE :: OUTPUT_PLT(int it)
{
    if( it%FREQ_CARTESIAN == 0 )   OUTPUT_PLT_CARTESIAN(it);
    if( it%FREQ_MESHLESS == 0 )    OUTPUT_PLT_MESHLESS(it);
    if( it%FREQ_REALSURFACE == 0 ) OUTPUT_PLT_REALSURFACE(it);
}

void WRITE_FILE :: OUTPUT_PLT_CARTESIAN(int iter)
{
    char *filename = new char [64];
    char Char_Buffer[33];
    strcpy(filename, "../result/Cartesian-");
    sprintf(Char_Buffer, "%d", int(iter/FREQ_CARTESIAN));
    strcat(filename, Char_Buffer);
    strcat(filename, ".plt");

    INTEGER4 FileType = 0;
    INTEGER4 Debug = 1;
    INTEGER4 VIsDouble = 1;
//    INTEGER4 I = 0; /* Used to track return codes */

    char TITLE[30] = "Cartesian Point";
    char VARIABLES[30] = "X Y Z U V W P";
    char SCRATCH[30] = ".";
    /* Open the file and write the tecplot datafile header information */
    TECINI112(TITLE, /* TITLE */
	      VARIABLES, /* VARIABLES */
	      filename, /* FILENAME */
	      SCRATCH, /* Scratch Directory */
	      &FileType,
	      &Debug,
	      &VIsDouble);
	
    /* Ordered Zone Parameters */
    INTEGER4 ZoneType = 0; 
    INTEGER4 IMax = mesh->IPOINT;
    INTEGER4 JMax = mesh->JPOINT;
    INTEGER4 KMax = mesh->KPOINT;
    INTEGER4 ICellMax = 0;
    INTEGER4 JCellMax = 0;
    INTEGER4 KCellMax = 0;
    double SolTime = global->Time;
    INTEGER4 StrandID = 100;
    INTEGER4 ParentZn = 0;
    INTEGER4 IsBlock = 1;  /* 1 - Block */
    INTEGER4 NFConns = 0;
    INTEGER4 FNMode = 0;
    INTEGER4 TotalNumFaceNodes = 1;
    INTEGER4 TotalNumBndryFaces = 1;
    INTEGER4 TotalNumBndryConnections = 1;
    INTEGER4 *PassiveVarArray = NULL;
    INTEGER4 *ValueLocArray = NULL;
    INTEGER4 *VarShareArray = NULL;
    INTEGER4 ShrConn = 0;
	
    /* Ordered Zone */
    TECZNE112(TITLE,  /* ZoneTitle */
	      &ZoneType,
	      &IMax,
	      &JMax,
	      &KMax,
	      &ICellMax,
	      &JCellMax,
	      &KCellMax,
	      &SolTime, /* SolutionTime */
	      &StrandID,
	      &ParentZn,
	      &IsBlock,
	      &NFConns,
	      &FNMode,
	      &TotalNumFaceNodes,
	      &TotalNumBndryFaces,
	      &TotalNumBndryConnections,
	      PassiveVarArray,
	      ValueLocArray,
	      VarShareArray,
	      &ShrConn);
	
    INTEGER4 DIsDouble = 1;
    INTEGER4 NumPts = mesh->POINT_CARTESIAN;

    TECDAT112(&NumPts, thrust::raw_pointer_cast( X.data() ), &DIsDouble);
    TECDAT112(&NumPts, thrust::raw_pointer_cast( Y.data() ), &DIsDouble);
    TECDAT112(&NumPts, thrust::raw_pointer_cast( Z.data() ), &DIsDouble);

    cusp::array1d<double,cusp::host_memory> TEMPU(mesh->POINT_CARTESIAN), TEMPV(mesh->POINT_CARTESIAN), TEMPW(mesh->POINT_CARTESIAN), TEMPP(mesh->POINT_CARTESIAN);
    thrust::copy(flow_field->U.row(0).begin(), flow_field->U.row(0).begin() + mesh->POINT_CARTESIAN, TEMPU.begin());
    thrust::copy(flow_field->V.row(0).begin(), flow_field->V.row(0).begin() + mesh->POINT_CARTESIAN, TEMPV.begin());
    thrust::copy(flow_field->W.row(0).begin(), flow_field->W.row(0).begin() + mesh->POINT_CARTESIAN, TEMPW.begin());
    thrust::copy(flow_field->P.row(0).begin(), flow_field->P.row(0).begin() + mesh->POINT_CARTESIAN, TEMPP.begin());
    
    TECDAT112(&NumPts, thrust::raw_pointer_cast( TEMPU.data() ), &DIsDouble);
    TECDAT112(&NumPts, thrust::raw_pointer_cast( TEMPV.data() ), &DIsDouble);
    TECDAT112(&NumPts, thrust::raw_pointer_cast( TEMPW.data() ), &DIsDouble);
    TECDAT112(&NumPts, thrust::raw_pointer_cast( TEMPP.data() ), &DIsDouble);

    TECEND112();
    delete filename;
}

void WRITE_FILE :: OUTPUT_PLT_MESHLESS(int iter)
{
    char *filename = new char [64];	
    char Char_Buffer[33];
    strcpy(filename, "../result/Meshless-");
    sprintf(Char_Buffer, "%d", int(iter/FREQ_MESHLESS));
    strcat(filename, Char_Buffer);
    strcat(filename, ".plt");

    INTEGER4 FileType = 0;
    INTEGER4 Debug = 1;
    INTEGER4 VIsDouble = 1;
//    INTEGER4 I = 0; /* Used to track return codes */
    
     char TITLE[30] = "Meshless Point";
     char VARIABLES[30] = "X Y Z U V W P";
     char SCRATCH[30] = ".";
     
     TECINI112(	TITLE, /* TITLE */
		VARIABLES, /* VARIABLES */
		filename, /* FILENAME */
		SCRATCH, /* Scratch Directory */
		&FileType,
		&Debug,
		&VIsDouble);
    /* Zone Parameters */
    INTEGER4 ZoneType;
    INTEGER4 NumNodes, NumElems, NumFaces = 0;
    INTEGER4 ICellMax = 0, JCellMax = 0, KCellMax = 0;
    double SolTime = global->Time;
    INTEGER4 StrandID = 0;
    INTEGER4 ParentZn = 0;
    INTEGER4 IsBlock = 1;  /* 1 - Block */
    INTEGER4 NFConns = 0;
    INTEGER4 FNMode = 0;
    INTEGER4 TotalNumFaceNodes = 1;
    INTEGER4 TotalNumBndryFaces = 1;
    INTEGER4 TotalNumBndryConnections = 1;
//    INTEGER4 PassiveVarArray[N_VARS];
//    INTEGER4 ValueLocArray[N_VARS];
    INTEGER4 VarShareArray[N_VARS];
    INTEGER4 ShrConn = 0;
	
    /* Data Parameters*/
    INTEGER4 DIsDouble = 1;
	
    /*Meshless Points*/
    ZoneType = 4;
    for(int i = 0; i < global->obj_number; i++) 
    {
	NumNodes = global->rigid_body[i]->POINT_NUMBER;
	NumElems = global->rigid_body[i]->TETRAHEDRON_NUMBER;
	StrandID = i+1;
		
	sprintf(Char_Buffer, "%d", StrandID);
	TECZNE112(	Char_Buffer,  /* ZoneTitle */
			&ZoneType,
			&NumNodes,
			&NumElems,
			&NumFaces,
			&ICellMax,
			&JCellMax,
			&KCellMax,
			&SolTime, /* SolutionTime */
			&StrandID,
			&ParentZn,
			&IsBlock,
			&NFConns,
			&FNMode,
			&TotalNumFaceNodes,
			&TotalNumBndryFaces,
			&TotalNumBndryConnections,
			NULL, // PassiveVarArray
			NULL, // ValueLocArray
			NULL, // VarShareArray
			&ShrConn);
	
	cusp::array1d<double,cusp::host_memory> TEMP;
	TEMP.resize(NumNodes);
	
	TEMP = global->rigid_body[i]->XYZ.row(0);
	TECDAT112(&NumNodes, thrust::raw_pointer_cast( TEMP.data() ), &DIsDouble);
	
	TEMP = global->rigid_body[i]->XYZ.row(1);
	TECDAT112(&NumNodes, thrust::raw_pointer_cast( TEMP.data() ), &DIsDouble);
	
	TEMP = global->rigid_body[i]->XYZ.row(2);
	TECDAT112(&NumNodes, thrust::raw_pointer_cast( TEMP.data() ), &DIsDouble);
	
	thrust::copy(flow_field->U.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i], 
		     flow_field->U.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i] + NumNodes,
		     TEMP.begin());
	TECDAT112(&NumNodes, thrust::raw_pointer_cast( TEMP.data() ), &DIsDouble);
	
	thrust::copy(flow_field->V.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i], 
		     flow_field->V.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i] + NumNodes,
		     TEMP.begin());
	TECDAT112(&NumNodes, thrust::raw_pointer_cast( TEMP.data() ), &DIsDouble);
	
	thrust::copy(flow_field->W.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i], 
		     flow_field->W.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i] + NumNodes,
		     TEMP.begin());
	TECDAT112(&NumNodes, thrust::raw_pointer_cast( TEMP.data() ), &DIsDouble);

	thrust::copy(flow_field->P.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i], 
		     flow_field->P.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i] + NumNodes,
		     TEMP.begin());
	TECDAT112(&NumNodes, thrust::raw_pointer_cast( TEMP.data() ), &DIsDouble);

	TECNOD112( thrust::raw_pointer_cast( TETRAHEDRON[i].values.data() ) );
    }
    
    /* Inner Surface */
    ZoneType = 2;
    for(int i=0;i<global->obj_number;i++) 
    {
	NumNodes = global->rigid_body[i]->POINT_NUMBER;
	NumElems = global->rigid_body[i]->TRIANGLE_NUMBER;
	StrandID = global->obj_number+i+1;
	for(int j=0;j<N_VARS;j++) VarShareArray[j] = i+1;
		
	sprintf(Char_Buffer, "%d", StrandID);
	TECZNE112(	Char_Buffer,  // ZoneTitle 
			&ZoneType,
			&NumNodes,
			&NumElems,
			&NumFaces,
			&ICellMax,
			&JCellMax,
			&KCellMax,
			&SolTime, // SolutionTime 
			&StrandID,
			&ParentZn,
			&IsBlock,
			&NFConns,
			&FNMode,
			&TotalNumFaceNodes,
			&TotalNumBndryFaces,
			&TotalNumBndryConnections,
			NULL, // PassiveVarArray
			NULL, // ValueLocArray
			VarShareArray,
			&ShrConn);
		
//	cusp::array2d<int,cusp::host_memory,cusp::column_major> TEMP2d(global->rigid_body[i]->TRIANGLE);
//	TECNOD112( thrust::raw_pointer_cast( TEMP2d.values.data() ) );
	TECNOD112( thrust::raw_pointer_cast( TRIANGLE[i].values.data() ) );
    }
	
    TECEND112();
    delete filename;
}

void WRITE_FILE :: OUTPUT_PLT_REALSURFACE(int iter)
{
    char *filename = new char [64];	
    char Char_Buffer[33];
    strcpy(filename, "../result/RealSurface-");
    sprintf(Char_Buffer, "%d", int(iter/FREQ_REALSURFACE));
    strcat(filename, Char_Buffer);
    strcat(filename, ".plt");
    
    INTEGER4 FileType = 0;
    INTEGER4 Debug = 1;
    INTEGER4 VIsDouble = 1;
//    INTEGER4 I = 0; /* Used to track return codes */

    char TITLE[30] = "Real Surface Points";
    char VARIABLES[30] = "X Y Z U V W P";
    char SCRATCH[30] = ".";

    TECINI112(TITLE, /* TITLE */
	      VARIABLES, /* VARIABLES */
	      filename, /* FILENAME */
	      SCRATCH, /* Scratch Directory */
	      &FileType,
	      &Debug,
	      &VIsDouble);
	
    /* Zone Parameters */
    INTEGER4 ZoneType;
    INTEGER4 NumNodes, NumElems, NumFaces = 0;
    INTEGER4 ICellMax = 0, JCellMax = 0, KCellMax = 0;
    double SolTime = global->Time;
    INTEGER4 StrandID = 0;
    INTEGER4 ParentZn = 0;
    INTEGER4 IsBlock = 1;  /* 1 - Block */
    INTEGER4 NFConns = 0;
    INTEGER4 FNMode = 0;
    INTEGER4 TotalNumFaceNodes = 1;
    INTEGER4 TotalNumBndryFaces = 1;
    INTEGER4 TotalNumBndryConnections = 1;
//    INTEGER4 PassiveVarArray[N_VARS];
//    INTEGER4 ValueLocArray[N_VARS];
//    INTEGER4 VarShareArray[N_VARS];
    INTEGER4 ShrConn = 0;
	
    /* Data Parameters*/
    INTEGER4 DIsDouble = 1;
	
    /* Inner Surface */
    ZoneType = 2;
    for(int i = 0; i < global->obj_number; i++)
    {
	NumNodes = global->rigid_body[i]->POINT_NUMBER;
	NumElems = global->rigid_body[i]->TRIANGLE_NUMBER;
	StrandID = global->obj_number+i+1;
		
	sprintf(Char_Buffer, "%d", StrandID);
	TECZNE112(Char_Buffer,  // ZoneTitle 
		  &ZoneType,
		  &NumNodes,
		  &NumElems,
		  &NumFaces,
		  &ICellMax,
		  &JCellMax,
		  &KCellMax,
		  &SolTime, // SolutionTime 
		  &StrandID,
		  &ParentZn,
		  &IsBlock,
		  &NFConns,
		  &FNMode,
		  &TotalNumFaceNodes,
		  &TotalNumBndryFaces,
		  &TotalNumBndryConnections,
		  NULL, // PassiveVarArray
		  NULL, // ValueLocArray
		  NULL, // VarShareArray
		  &ShrConn);
		
	cusp::array1d<double,cusp::host_memory> TEMP;
	TEMP.resize(NumNodes);
	
	TEMP = global->rigid_body[i]->XYZ.row(0);
	TECDAT112(&NumNodes, thrust::raw_pointer_cast( TEMP.data() ), &DIsDouble);
	
	TEMP = global->rigid_body[i]->XYZ.row(1);
	TECDAT112(&NumNodes, thrust::raw_pointer_cast( TEMP.data() ), &DIsDouble);
	
	TEMP = global->rigid_body[i]->XYZ.row(2);
	TECDAT112(&NumNodes, thrust::raw_pointer_cast( TEMP.data() ), &DIsDouble);
	
	thrust::copy(flow_field->U.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i], 
		     flow_field->U.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i] + NumNodes,
		     TEMP.begin());
	TECDAT112(&NumNodes, thrust::raw_pointer_cast( TEMP.data() ), &DIsDouble);
	
	thrust::copy(flow_field->V.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i], 
		     flow_field->V.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i] + NumNodes,
		     TEMP.begin());
	TECDAT112(&NumNodes, thrust::raw_pointer_cast( TEMP.data() ), &DIsDouble);
	
	thrust::copy(flow_field->W.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i], 
		     flow_field->W.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i] + NumNodes,
		     TEMP.begin());
	TECDAT112(&NumNodes, thrust::raw_pointer_cast( TEMP.data() ), &DIsDouble);

	thrust::copy(flow_field->P.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i], 
		     flow_field->P.row(0).begin() + mesh->POINT_CARTESIAN + mesh->meshless->OFFSET[i] + NumNodes,
		     TEMP.begin());
	TECDAT112(&NumNodes, thrust::raw_pointer_cast( TEMP.data() ), &DIsDouble);
	
	TECNOD112( thrust::raw_pointer_cast( TRIANGLE[i].values.data() ) );
    }
    TECEND112();
    delete filename;
}


#endif
