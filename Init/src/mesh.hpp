#include "../include/mesh.h"

struct type_functor : public thrust::unary_function<int,int>
{
    const int IPOINT, JPOINT, KPOINT;
    
    type_functor(int _IPOINT, int _JPOINT, int _KPOINT) : IPOINT(_IPOINT), JPOINT(_JPOINT), KPOINT(_KPOINT) {}
    
    __host__ __device__
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

template<typename MemorySpace>
MESH_CARTESIAN<MemorySpace> :: MESH_CARTESIAN()
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

    IINDEX.resize(IPOINT);
    JINDEX.resize(JPOINT);
    
    for(int i = 0; i < 3; i++) POINTTYPE[i].resize(POINT_CARTESIAN);
    
    COUNTTYPE.resize(POINT_CARTESIAN);
    
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

template<typename MemorySpace>
MESH_CARTESIAN<MemorySpace> :: ~MESH_CARTESIAN()
{
    std::cout << "Cartesian mesh is being deleted." << std::endl;
}
    
template<typename MemorySpace>
void MESH_CARTESIAN<MemorySpace> :: Div(double start, int segment, double *segmentlength, double *ratio, int *segmentpoint, int POINT, double *temp)
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









template<typename MemorySpace>
MESH_LESS<MemorySpace> :: MESH_LESS(int Point_Messless)
{
	std::cout << "Meshless is being initialized......";
	
	POINT_MESSLESS = Point_Messless;
	
	for(int i = 0; i < 3; i++)
	{
	    POSITION[i].resize(POINT_MESSLESS);
	    VELOCITY[i].resize(POINT_MESSLESS);
	    ACCELERATION[i].resize(POINT_MESSLESS);
	}
	
	std::cout << "Initialization done." << std::endl;
}

template<typename MemorySpace>
MESH_LESS<MemorySpace> :: ~MESH_LESS()
{
    std::cout << "Meshless is being deleted." << std::endl;
}