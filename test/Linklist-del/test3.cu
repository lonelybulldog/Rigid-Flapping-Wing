#include <cusp/array1d.h>
#include <cusp/print.h>
#include <thrust/sequence.h>
#include <thrust/functional.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/remove.h>
#include <iostream>

#define NB 40

//class LINKLIST;  /*try deleted*/

class MESH
{
public:
    //thrust::device_vector<int> stencil;
    cusp::array1d<int,cusp::device_memory> stencil;
    MESH()
    {
        std::cout << "Mesh is created...." << std::endl;
	stencil.resize(10);
        thrust::sequence(stencil.begin(),stencil.end(),1,1);
    }
    ~MESH()
    {
        std::cout << "Mesh is deleted." << std::endl;
    }
    friend class LINKLIST;
};

class FLOW_FIELD
{
public:
    cusp::array2d<int, cusp::device_memory> U;
    FLOW_FIELD()
    {
	std::cout << "Flow field is created...." << std::endl;
	U.resize(2,10);
	U(0,0)=10; U(0,1)=10; U(0,2)=10; U(0,3)=10; U(0,4)=10; U(0,5)=10; U(0,6)=10; U(0,7)=10; U(0,8)=10; U(0,9)=10;
	U(1,0)=20; U(1,1)=20; U(1,2)=20; U(1,3)=20; U(1,4)=20; U(1,5)=20; U(1,6)=20; U(1,7)=20; U(1,8)=20; U(1,9)=20;
//	thrust::sequence(U.begin(),U.end(),10,10);
    }
    ~FLOW_FIELD()
    {
	std::cout << "Flow field is deleted." << std::endl;
    }
    friend class LINKLIST;
};


typedef struct LINKLIST_MEMBER
{
    int Meshless_Ind;
    int Nb_Points[NB];
    int Nb_Points_Old[NB];
    double Csvd[9][NB];
    double Csvd_Old[9][NB];
}LINKLIST_MEMBER;


struct is_del
{
//    thrust::device_vector<int> stencil;
//    thrust::omp::vector<int> stencil;
    int *raw;
    int *U_raw;
  
    is_del(MESH& mesh, FLOW_FIELD& flow_field)
    {
	raw=thrust::raw_pointer_cast(mesh.stencil.data());
	
//	U_raw=thrust::raw_pointer_cast(flow_field.U.data());
	U_raw=thrust::raw_pointer_cast(  &flow_field.U.row(1)[0]  );
//	stencil.resize(10);
//	thrust::sequence(stencil.begin(),stencil.end(),0,1);
    }
  
    __host__ __device__
    bool operator()(const LINKLIST_MEMBER& x)
    {
//	return (x.Meshless_Ind % 2) == 0;
	bool s=(raw[x.Meshless_Ind]) % 2 == 0;
	if(s) U_raw[x.Meshless_Ind]=0;
	
	return s;
    }
};

  

class LINKLIST
{
public:
//    thrust::device_vector<LINKLIST_MEMBER> LINKLISTPOINT;
    cusp::array1d<LINKLIST_MEMBER,cusp::device_memory> LINKLISTPOINT;
    LINKLIST()
    {
        std::cout << "Linklist is created...." << std::endl;
    }
    ~LINKLIST()
    {
        std::cout << "Linklist is deleted." << std::endl;
    }
    void Insert_point(int i)
    {
        LINKLIST_MEMBER temp;
        temp.Meshless_Ind = i;
        LINKLISTPOINT.push_back(temp);
    }
    
    void Delete_point(MESH& mesh, FLOW_FIELD& flow_field)
    {
//	thrust::remove_if(LINKLISTPOINT.begin(), LINKLISTPOINT.end(), is_del());
//	is_del ISDEL;
        LINKLISTPOINT.erase(thrust::remove_if(LINKLISTPOINT.begin(),
                                              LINKLISTPOINT.end(),
					      is_del(mesh,flow_field)), LINKLISTPOINT.end());
    }

};





int main()
{
      int i=3, j=6, k=7;
      
      MESH mesh;
      FLOW_FIELD flow_field;
      LINKLIST linklist;
      linklist.Insert_point(i);
      linklist.Insert_point(j);
      linklist.Insert_point(k);
      
//      linklist.Delete_point(mesh,flow_field);
      
      std::cout<<linklist.LINKLISTPOINT.size()<<std::endl;
      std::cout<<linklist.LINKLISTPOINT[2].Meshless_Ind<<std::endl;
      std::cout<<" it contains: "<< std::endl;
//      for(thrust::device_vector<LINKLIST_MEMBER>::iterator iter = linklist.LINKLISTPOINT.begin(); iter != linklist.LINKLISTPOINT.end(); iter++)  {
//                std::cout << (static_cast<LINKLIST_MEMBER>(*iter)).Meshless_Ind << std::endl;
//        }

//      cusp::print(flow_field.U);
      
}