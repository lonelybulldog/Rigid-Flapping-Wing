

template<typename MemorySpace> class FLOW_FIELD
{
public:
    int    POINT_ALL, DERIVATIVE;
    
    cusp::array2d<double,MemorySpace,cusp::row_major> U;
    cusp::array2d<double,MemorySpace,cusp::row_major> V;
    cusp::array2d<double,MemorySpace,cusp::row_major> W;
    cusp::array2d<double,MemorySpace,cusp::row_major> P;
    
    cusp::array2d<double,MemorySpace,cusp::row_major> U_old;
    cusp::array2d<double,MemorySpace,cusp::row_major> V_old;
    cusp::array2d<double,MemorySpace,cusp::row_major> W_old;
    cusp::array2d<double,MemorySpace,cusp::row_major> P_old;
        
    cusp::array2d<double,MemorySpace,cusp::row_major> U_temp;
    cusp::array2d<double,MemorySpace,cusp::row_major> V_temp;
    cusp::array2d<double,MemorySpace,cusp::row_major> W_temp;
    cusp::array2d<double,MemorySpace,cusp::row_major> P_temp;
        
    FLOW_FIELD(int Point_Car, int Point_Messless)
    {
        std::cout << "Flow field is being initialized......";
        
        POINT_ALL = Point_Car + Point_Messless;
        DERIVATIVE = 2;
	
        U.resize(DERIVATIVE,POINT_ALL);
	V.resize(DERIVATIVE,POINT_ALL);
	W.resize(DERIVATIVE,POINT_ALL);
	P.resize(DERIVATIVE,POINT_ALL);
	
	U_old.resize(DERIVATIVE,POINT_ALL);
	V_old.resize(DERIVATIVE,POINT_ALL);
	W_old.resize(DERIVATIVE,POINT_ALL);
	P_old.resize(DERIVATIVE,POINT_ALL);
	
	U_temp.resize(DERIVATIVE,POINT_ALL);
	V_temp.resize(DERIVATIVE,POINT_ALL);
	W_temp.resize(DERIVATIVE,POINT_ALL);
	P_temp.resize(DERIVATIVE,POINT_ALL);

	std::cout << "Initialization done." << std::endl;
    }
    
    ~FLOW_FIELD()
    {
        std::cout << "Flow field is being deleted." << std::endl;
    }

};

