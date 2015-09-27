#include "../include/mesh.h"
#include <thrust/adjacent_difference.h>
#include <thrust/sequence.h>
#include <thrust/reduce.h>
#include <iostream>
#include <fstream>

struct type_functor : public thrust::unary_function<int,int>
{
    const int IPOINT, JPOINT, KPOINT;
    
    type_functor(int _IPOINT, int _JPOINT, int _KPOINT) : IPOINT(_IPOINT), JPOINT(_JPOINT), KPOINT(_KPOINT) {}
    
    __device__
    int operator() (int INDEX) const
    {
	const int i = INDEX / (JPOINT * KPOINT);
	const int j = (INDEX % (JPOINT * KPOINT)) / KPOINT;
	const int k = (INDEX % (JPOINT * KPOINT)) % KPOINT;
	
	if( (i == 0) || (j == 0) || (k == 0) || (i == (IPOINT-1)) || (j == (JPOINT-1)) || (k == (KPOINT-1)) )
	{
	    if( ((i==0||i==IPOINT-1)&&(j==0)&&(0<=k)&&(k<=KPOINT-1))||
		((i==0||i==IPOINT-1)&&(j==JPOINT-1)&&(0<=k)&&(k<=KPOINT-1))||
		((i==0||i==IPOINT-1)&&(1<=j)&&(j<=JPOINT-2)&&(k==0))||
		((i==0||i==IPOINT-1)&&(1<=j)&&(j<=JPOINT-2)&&(k==KPOINT-1))||
		((1<=i)&&(i<=IPOINT-2)&&(j==0)&&(k==0))||
		((1<=i)&&(i<=IPOINT-2)&&(j==0)&&(k==KPOINT-1))||
		((1<=i)&&(i<=IPOINT-2)&&(j==JPOINT-1)&&(k==0))||
		((1<=i)&&(i<=IPOINT-2)&&(j==JPOINT-1)&&(k==KPOINT-1)) )
	    {
		return 6;
	    }
	    else
	    {
		return 2;
	    }
	}
	else
	{
	    return 1;
	}
    }
};


MESH_CARTESIAN :: MESH_CARTESIAN()
{
    std::cout << "Cartesian mesh is being initialized......";
	
    double Xstart, Xend, Ystart, Yend, Zstart, Zend;
    int Xsegment, Ysegment, Zsegment, XPoint = 1, YPoint = 1, ZPoint = 1;
        
    std::ifstream read;
    read.open("../flapping/cart.dat");
    if(!read.is_open()) std::cout << "Cart.dat is unable to open." << std::endl;
    read >> IPOINT >> JPOINT >> KPOINT;
    POINT_CARTESIAN = IPOINT * JPOINT * KPOINT;
    
    read >> Xstart >> Xsegment;
    double Xsegmentlength[Xsegment], Xratio[Xsegment];
    int Xsegmentpoint[Xsegment];
    Xend = Xstart;
    for(int i = 0; i < Xsegment; i++)
    {
	read >> Xsegmentlength[i] >> Xratio[i] >> Xsegmentpoint[i];
	Xend = Xend + Xsegmentlength[i];
	XPoint = XPoint + Xsegmentpoint[i];
    }
    if(XPoint != IPOINT) std::cout << "Setting of Cartesian mesh is wrong: X_direction." << std::endl;

    read >> Ystart >> Ysegment;
    double Ysegmentlength[Ysegment], Yratio[Ysegment];
    int Ysegmentpoint[Ysegment];
    Yend = Ystart;
    for(int i = 0; i < Ysegment; i++)
    {
	read >> Ysegmentlength[i] >> Yratio[i] >> Ysegmentpoint[i];
	Yend = Yend + Ysegmentlength[i];
	YPoint = YPoint + Ysegmentpoint[i];
    }
    if(YPoint != JPOINT) std::cout << "Setting of Cartesian mesh is wrong: Y_direction" << std::endl;
	
    read >> Zstart >> Zsegment;
    double Zsegmentlength[Zsegment], Zratio[Zsegment];
    int Zsegmentpoint[Zsegment];
    Zend = Zstart;
    for(int i = 0; i < Zsegment; i++)
    {
	read >> Zsegmentlength[i] >> Zratio[i] >> Zsegmentpoint[i];
	Zend = Zend + Zsegmentlength[i];
	ZPoint = ZPoint + Zsegmentpoint[i];
    }
    if(ZPoint != KPOINT) std::cout << "Setting of Cartesian mesh is wrong: Z_direction" << std::endl;
        
    read.close();
	
    XYZ[0].resize(IPOINT);
    XYZ[1].resize(JPOINT);
    XYZ[2].resize(KPOINT);
    
    DELTA[0].resize(IPOINT);
    DELTA[1].resize(JPOINT);
    DELTA[2].resize(KPOINT);
    
    VELOCITY.resize(3);
    ACCELERATION.resize(3);
    
    IINDEX.resize(IPOINT);
    JINDEX.resize(JPOINT);
    
    for(int i = 0; i < 3; i++) POINTTYPE[i].resize(POINT_CARTESIAN);
    
    COUNTTYPE.resize(POINT_CARTESIAN);  /*Need Initialization*/
    
    double tempX[IPOINT], tempY[JPOINT], tempZ[KPOINT];
    Div(Xstart, Xsegment, Xsegmentlength, Xratio, Xsegmentpoint, IPOINT, tempX);
    Div(Ystart, Ysegment, Ysegmentlength, Yratio, Ysegmentpoint, JPOINT, tempY);
    Div(Zstart, Zsegment, Zsegmentlength, Zratio, Zsegmentpoint, KPOINT, tempZ);
    thrust::copy(tempX, tempX+IPOINT, XYZ[0].begin());
    thrust::copy(tempY, tempY+JPOINT, XYZ[1].begin());
    thrust::copy(tempZ, tempZ+KPOINT, XYZ[2].begin());
    thrust::adjacent_difference(XYZ[0].begin(), XYZ[0].end(), DELTA[0].begin()); DELTA[0][0]=0;
    thrust::adjacent_difference(XYZ[1].begin(), XYZ[1].end(), DELTA[1].begin()); DELTA[1][0]=0;
    thrust::adjacent_difference(XYZ[2].begin(), XYZ[2].end(), DELTA[2].begin()); DELTA[2][0]=0;
    
    for(int i = 0; i < 3; i++)
    {
        VELOCITY[i]     = 0;
        ACCELERATION[i] = 0;
    }

    XMIN = XYZ[0][0]; XMAX = XYZ[0][IPOINT-1];
    YMIN = XYZ[1][0]; YMAX = XYZ[1][JPOINT-1];
    ZMIN = XYZ[2][0]; ZMAX = XYZ[2][KPOINT-1];
    MESHSIZE = pow( ( (XMAX-XMIN)/(IPOINT-1) * (YMAX-YMIN)/(JPOINT-1) * (ZMAX-ZMIN)/(KPOINT-1) ), 1.0/3.0 );
    SAFEDISTANCE = 0.1;

    thrust::sequence(IINDEX.begin(), IINDEX.end(), 0, KPOINT*JPOINT);
    thrust::sequence(JINDEX.begin(), JINDEX.end(), 0, KPOINT);
    
    thrust::transform(thrust::make_counting_iterator(0), 
                      thrust::make_counting_iterator(POINT_CARTESIAN),
                      POINTTYPE[0].begin(),
                      type_functor(IPOINT,JPOINT,KPOINT));
    thrust::copy(POINTTYPE[0].begin(), POINTTYPE[0].end(), POINTTYPE[1].begin());
    thrust::copy(POINTTYPE[0].begin(), POINTTYPE[0].end(), POINTTYPE[2].begin());
    
    
    
    
    
    
    std::cout << "Initialization done." << std::endl;
}


MESH_CARTESIAN :: ~MESH_CARTESIAN()
{
    std::cout << "Cartesian mesh is being deleted." << std::endl;
}
    

void MESH_CARTESIAN :: Div(double start, int segment, double *segmentlength, double *ratio, int *segmentpoint, int POINT, double *temp)
{
    double firstsegment;
    int offset = 1; 
    
    temp[0] = start;
    
    for(int i = 0; i < segment; i++)
    {
	if ( fabs(ratio[i] - 1 ) < 1e-4 )
	{
	    firstsegment = segmentlength[i] / segmentpoint[i];
	    
	    for (int j = 0; j < segmentpoint[i]; j++)
	    {
		start = start + firstsegment;
		
		temp[offset + j] = start;
	    }
	    
	    offset = offset + segmentpoint[i];
	}
	else
	{
	    firstsegment = segmentlength[i] * (1 - ratio[i]) / (1 - pow(ratio[i] , segmentpoint[i]));
	    
	    for (int j = 0; j < segmentpoint[i]; j++)
	    {
		start = start + firstsegment * pow(ratio[i], j);
		
		temp[offset + j] = start;
	    }
	    
	    offset = offset + segmentpoint[i];
	}
    }
}

struct is_inner
{
    __device__
    bool operator()(const int x)
    {
	return ( x == 1 );
    }
};

struct is_outer
{
    __device__
    bool operator()(const int x)
    {
	return ( x == 2 );
    }
};

MESH_LESS :: MESH_LESS( REF_FRAME& global )
{
    std::cout << "Meshless is being initialized......";

    POINT_MESSLESS = 0, INNER_POINT = 0, OUTER_POINT = 0;
    OFFSET = new int [global.obj_number];
    INNER_POINT_OFFSET = new int [global.obj_number];
    OUTER_POINT_OFFSET = new int [global.obj_number];
   
    for ( int i = 0; i < global.obj_number; i++ )
    {
	OFFSET[i] = POINT_MESSLESS;
	INNER_POINT_OFFSET[i] = INNER_POINT;
	OUTER_POINT_OFFSET[i] = OUTER_POINT;
	POINT_MESSLESS = POINT_MESSLESS + global.rigid_body[i]->POINT_NUMBER;
	INNER_POINT = INNER_POINT + global.rigid_body[i]->INNER_POINT_NUMBER;
	OUTER_POINT = OUTER_POINT + global.rigid_body[i]->OUTER_POINT_NUMBER;
	
    }
	
    for( int i = 0; i < 3; i++ )
    {
	POSITION[i].resize(POINT_MESSLESS);
	VELOCITY[i].resize(POINT_MESSLESS);
	ACCELERATION[i].resize(POINT_MESSLESS);
	INNER_NODE_INDEX.resize(INNER_POINT);
	OUTER_NODE_INDEX.resize(OUTER_POINT);
	SURFACE_ELE_AREA.resize(INNER_POINT);
    }
	
    std::cout<<"OUTER_POINT is "<<OUTER_POINT<<std::endl;
    for( int i = 0; i < global.obj_number; i++ )
    {
	for( int j = 0; j < 3; j++ )
	{
	    thrust::copy(global.rigid_body[i]->XYZ.row(j).begin(), global.rigid_body[i]->XYZ.row(j).end(), POSITION[j].begin() + OFFSET[i]);
	    thrust::copy(global.rigid_body[i]->UVW.row(j).begin(), global.rigid_body[i]->UVW.row(j).end(), VELOCITY[j].begin() + OFFSET[i]);
	    thrust::copy(global.rigid_body[i]->ACC.row(j).begin(), global.rigid_body[i]->ACC.row(j).end(), ACCELERATION[j].begin() + OFFSET[i]);
	}
	
	thrust::copy_if(thrust::make_counting_iterator(OFFSET[i]), 
			thrust::make_counting_iterator(OFFSET[i] + global.rigid_body[i]->POINT_NUMBER),
			global.rigid_body[i]->INNERMARK.begin(),
			INNER_NODE_INDEX.begin() + INNER_POINT_OFFSET[i],
			is_inner());
	
	thrust::copy_if(thrust::make_counting_iterator(OFFSET[i]),
			thrust::make_counting_iterator(OFFSET[i] + global.rigid_body[i]->POINT_NUMBER),
			global.rigid_body[i]->OUTERMARK.begin(),
			OUTER_NODE_INDEX.begin() + OUTER_POINT_OFFSET[i],
			is_outer());
	
	thrust::copy_if(global.rigid_body[i]->AREA.begin(),
			global.rigid_body[i]->AREA.end(),
			global.rigid_body[i]->INNERMARK.begin(),
			SURFACE_ELE_AREA.begin() + INNER_POINT_OFFSET[i],
			is_inner());
    }

    std::cout << "Initialization done." << std::endl;
}


MESH_LESS :: ~MESH_LESS()
{
    delete[] OFFSET;
    std::cout << "Meshless is being deleted." << std::endl;
}


MESH :: MESH( MESH_CARTESIAN& Cartesian, MESH_LESS& Meshless )
{
    std::cout << "Mesh is being initialized......";
    
    cartesian = &Cartesian;
    meshless  = &Meshless;

    IPOINT = cartesian->IPOINT;
    JPOINT = cartesian->JPOINT;
    KPOINT = cartesian->KPOINT;
    POINT_CARTESIAN = cartesian->POINT_CARTESIAN;
    POINT_MESSLESS = meshless->POINT_MESSLESS;
    POINT_ALL = POINT_CARTESIAN + POINT_MESSLESS;
    
    XYZ.IPOINT          = IPOINT;
    XYZ.JPOINT          = JPOINT;
    XYZ.KPOINT          = KPOINT;
    XYZ.POINT_CARTESIAN = POINT_CARTESIAN;
    XYZ.POINT_MESSLESS  = POINT_MESSLESS;
    XYZ.POINT_ALL       = POINT_CARTESIAN + POINT_MESSLESS;
    
    UVW.IPOINT          = IPOINT;
    UVW.JPOINT          = JPOINT;
    UVW.KPOINT          = KPOINT;
    UVW.POINT_CARTESIAN = POINT_CARTESIAN;
    UVW.POINT_MESSLESS  = POINT_MESSLESS;
    UVW.POINT_ALL       = POINT_CARTESIAN + POINT_MESSLESS;
    
    ACC.IPOINT          = IPOINT;
    ACC.JPOINT          = JPOINT;
    ACC.KPOINT          = KPOINT;
    ACC.POINT_CARTESIAN = POINT_CARTESIAN;
    ACC.POINT_MESSLESS  = POINT_MESSLESS;
    ACC.POINT_ALL       = POINT_CARTESIAN + POINT_MESSLESS;
        
    for( int s = 0; s < 3; s++ )
    {
        XYZ.XYZ_raw[s]                   = thrust::raw_pointer_cast( cartesian->XYZ[s].data() );
        XYZ.POSITION_MESHLESS_raw[s]     = thrust::raw_pointer_cast( meshless->POSITION[s].data() );
        
        UVW.VELOCITY_MESHLESS_raw[s]     = thrust::raw_pointer_cast( meshless->VELOCITY[s].data() );
        
        ACC.ACCELERATION_MESHLESS_raw[s] = thrust::raw_pointer_cast( meshless->ACCELERATION[s].data() );
    }
    UVW.VELOCITY_CARTESIAN_raw     = thrust::raw_pointer_cast( cartesian->VELOCITY.data() );
    ACC.ACCELERATION_CARTESIAN_raw = thrust::raw_pointer_cast( cartesian->ACCELERATION.data() );
    
    std::cout << "Initialization done." << std::endl;
}

MESH :: ~MESH()
{
    std::cout << "Mesh is being deleted." << std::endl;
}

void MESH :: UPDATE_MESHLESS( REF_FRAME& global )
{
    for( int i = 0; i < global.obj_number; i++ )
    {
	for( int j = 0; j < 3; j++ )
	{
	    thrust::copy(global.rigid_body[i]->XYZ.row(j).begin(), global.rigid_body[i]->XYZ.row(j).end(), meshless->POSITION[j].begin() + meshless->OFFSET[i]);
	    thrust::copy(global.rigid_body[i]->UVW.row(j).begin(), global.rigid_body[i]->UVW.row(j).end(), meshless->VELOCITY[j].begin() + meshless->OFFSET[i]);
	    thrust::copy(global.rigid_body[i]->ACC.row(j).begin(), global.rigid_body[i]->ACC.row(j).end(), meshless->ACCELERATION[j].begin() + meshless->OFFSET[i]);
	}
    }
}

struct key_functor_SEARCH_TYPE0 : public thrust :: unary_function<int, int>
{
    const int ONI_SIZE;
    
    key_functor_SEARCH_TYPE0(int _ONI_SIZE) : ONI_SIZE(_ONI_SIZE) {}
    
    __device__
    int operator()(int x) { return x/ONI_SIZE; }
};

typedef thrust::tuple<int, int, int, int, double> LOCATE_POINT;
/*
typedef struct LOCATE_POINT
{
    int    I;
    int    J;
    int    K;
    int    IM;
    double DISTANCE;
}LOCATE_POINT;
*/

struct functor_LOCATE_POINT : public thrust::unary_function<int, LOCATE_POINT>
{
    int XLENGTH, YLENGTH, ZLENGTH, ISTART, JSTART, KSTART, MLENGTH;
    MESH_CARTESIAN   *cartesian;
    MESH_LESS        *meshless;
    int              *OUTER_NODE_INDEX_raw;
    double           *cartesian_XYZ_raw[3];
    double           *meshless_XYZ_raw[3];
    
    
    functor_LOCATE_POINT(int _XLENGTH, int _YLENGTH, int _ZLENGTH, 
			 int _ISTART, int _JSTART, int _KSTART, int _MLENGTH, 
			 MESH_CARTESIAN *Cartesian, MESH_LESS *Meshless)
    {
	XLENGTH = _XLENGTH; YLENGTH = _YLENGTH; ZLENGTH = _ZLENGTH;
	ISTART = _ISTART; JSTART = _JSTART; KSTART = _KSTART; MLENGTH = _MLENGTH;
	cartesian = Cartesian; meshless = Meshless;
	
	OUTER_NODE_INDEX_raw = thrust::raw_pointer_cast( meshless->OUTER_NODE_INDEX.data() );
	
	for(int s = 0; s < 3; s++)
	{
	    cartesian_XYZ_raw[s] = thrust::raw_pointer_cast( cartesian->XYZ[s].data() );
	    meshless_XYZ_raw[s]  = thrust::raw_pointer_cast( meshless->POSITION[s].data() );
	}
    }
    
    __device__
    LOCATE_POINT operator() (int x) const
    {
	LOCATE_POINT locate_point;
	
	const int im = x % MLENGTH;
	const int i  = (x / MLENGTH) / (YLENGTH * ZLENGTH) + ISTART;
	const int j  = ( (x / MLENGTH) % (YLENGTH * ZLENGTH) ) / ZLENGTH + JSTART;
	const int k  = ( (x / MLENGTH) % (YLENGTH * ZLENGTH) ) % ZLENGTH + KSTART;
	
	double x1,x2,y1,y2,z1,z2;
//	int ONI = meshless->OUTER_NODE_INDEX[im];
	int ONI = OUTER_NODE_INDEX_raw[im];
	
//	x1 = cartesian->XYZ[0][i]; y1 = cartesian->XYZ[1][j]; z1 = cartesian->XYZ[2][k];
//	x2 = meshless->POSITION[0][ONI]; y2 = meshless->POSITION[1][ONI]; z2 = meshless->POSITION[2][ONI];
	
	x1 = cartesian_XYZ_raw[0][i]; y1 = cartesian_XYZ_raw[1][j]; z1 = cartesian_XYZ_raw[2][k];
	x2 = meshless_XYZ_raw[0][ONI]; y2 = meshless_XYZ_raw[1][ONI]; z2 = meshless_XYZ_raw[2][ONI];
	
	const double d = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );
	
	thrust::get<0>(locate_point) = i;
	thrust::get<1>(locate_point) = j;
	thrust::get<2>(locate_point) = k;
	thrust::get<3>(locate_point) = im;
	thrust::get<4>(locate_point) = d;
	
	return locate_point;
/*
	locate_point.I = i;
	locate_point.J = j;
	locate_point.K = k;
	locate_point.IM = im;
	locate_point.DISTANCE = d;
	return locate_point;
*/
    }
};

struct nearest_point : public thrust::binary_function<LOCATE_POINT, LOCATE_POINT, LOCATE_POINT>
{
    __device__
    LOCATE_POINT operator()(LOCATE_POINT x, LOCATE_POINT y)
    {
	if( thrust::get<4>(x) < thrust::get<4>(y) ) return x;
//	if( x.DISTANCE < y.DISTANCE ) return x;
	else return y;
    }
};


void MESH :: SEARCH_TYPE0()
{
    double Xmax,Ymax,Zmax,Xmin,Ymin,Zmin;
    int Istart, Iend, Jstart,Jend, Kstart, Kend;
    std::cout<<"Reach here 0"<<std::endl;
    thrust::pair<thrust::device_vector<double>::iterator, thrust::device_vector<double>::iterator> minmax;
    
    minmax = thrust::minmax_element( meshless->POSITION[0].begin(), meshless->POSITION[0].end() );
    Xmin = *minmax.first; Xmax = *minmax.second;
    
    minmax = thrust::minmax_element( meshless->POSITION[1].begin(), meshless->POSITION[1].end() );
    Ymin = *minmax.first; Ymax = *minmax.second;
    
    minmax = thrust::minmax_element( meshless->POSITION[2].begin(), meshless->POSITION[2].end() );
    Zmin = *minmax.first; Zmax = *minmax.second;
    std::cout<<"Reach here 0.1"<<std::endl;
    if( Xmin <= ( cartesian->XMIN + cartesian->SAFEDISTANCE ) || Xmax >= ( cartesian->XMAX - cartesian->SAFEDISTANCE ) ||
        Ymin <= ( cartesian->YMIN + cartesian->SAFEDISTANCE ) || Ymax >= ( cartesian->YMAX - cartesian->SAFEDISTANCE ) ||
        Zmin <= ( cartesian->ZMIN + cartesian->SAFEDISTANCE ) || Zmax >= ( cartesian->ZMAX - cartesian->SAFEDISTANCE ) )
    {
	std::cout<<"Meshless points are out of range! Need bigger box! Code terminate!"<<std::endl;
    }
    std::cout<<"Reach here 0.2"<<std::endl;
    Kstart = int( (Zmin - cartesian->ZMIN)/cartesian->MESHSIZE ) - 1;
    Kend   = int( (Zmax - cartesian->ZMIN)/cartesian->MESHSIZE ) + 2;
    Jstart = int( (Ymin - cartesian->YMIN)/cartesian->MESHSIZE ) - 1;
    Jend   = int( (Ymax - cartesian->YMIN)/cartesian->MESHSIZE ) + 2;
    Istart = int( (Xmin - cartesian->XMIN)/cartesian->MESHSIZE ) - 1;
    Iend   = int( (Xmax - cartesian->XMIN)/cartesian->MESHSIZE ) + 2;
    std::cout<<"Reach here 0.24"<<std::endl;
    thrust::host_vector<LOCATE_POINT> ijk_im_d_host( (Iend-Istart+1)*(Jend-Jstart+1)*(Kend-Kstart+1)*(meshless->OUTER_POINT) );
    //(Iend-Istart+1)*(Jend-Jstart+1)*(Kend-Kstart+1)*meshless->OUTER_POINT
//    cusp::array1d<LOCATE_POINT,cusp::host_memory> ijk_im_d_host;
    std::cout<<"Reach here 0.25"<<std::endl;
//    ijk_im_d_host.resize((Iend-Istart+1)*(Jend-Jstart+1)*(Kend-Kstart+1)*meshless->OUTER_POINT);
    std::cout<<"Reach here 0.26"<<std::endl;
    cusp::array1d<LOCATE_POINT,cusp::device_memory> ijk_im_d_device = ijk_im_d_host;
    
    std::cout<<"Reach here 0.3"<<std::endl;
    thrust::device_vector<LOCATE_POINT> ijk_im_d_selected((Iend-Istart+1)*(Jend-Jstart+1)*(Kend-Kstart+1));
    
    std::cout<<"Reach here 1"<<std::endl;
    
    thrust::transform( thrust::make_counting_iterator(0), 
		       thrust::make_counting_iterator( (Iend-Istart+1)*(Jend-Jstart+1)*(Kend-Kstart+1)*meshless->OUTER_POINT ),
		       ijk_im_d_device.begin(),
		       functor_LOCATE_POINT(Iend-Istart+1, Jend-Jstart+1, Kend-Kstart+1, Istart, Jstart, Kstart, meshless->OUTER_POINT, cartesian, meshless) );
    
    std::cout<<"Reach here 2"<<std::endl;
    
    thrust::equal_to<int> binary_pred;
    thrust::device_vector<int> key_value( (Iend-Istart+1)*(Jend-Jstart+1)*(Kend-Kstart+1) );
    
    thrust::reduce_by_key( thrust::make_transform_iterator( thrust::make_counting_iterator(0), key_functor_SEARCH_TYPE0( meshless->OUTER_POINT ) ), 
			   thrust::make_transform_iterator( thrust::make_counting_iterator( (Iend-Istart+1)*(Jend-Jstart+1)*(Kend-Kstart+1)*meshless->OUTER_POINT ), key_functor_SEARCH_TYPE0( meshless->OUTER_POINT ) ),
			   ijk_im_d_device.begin(),
//			   thrust::make_discard_iterator(),
			   key_value.begin(),
			   ijk_im_d_selected.begin(),
			   binary_pred,
			   nearest_point() );
    
    
    std::cout<<Istart<<" "<<Iend<<" "<<Jstart<<" "<<Jend<<" "<<Kstart<<" "<<Kend<<std::endl;
    
    
}

