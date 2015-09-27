#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sequence.h>
#include <thrust/functional.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/remove.h>
#include <iostream>
#include <thrust/system/omp/vector.h>

#define NB 40

//class LINKLIST;
//struct is_del;

class MESH
{
public:
    thrust::device_vector<int> stencil;
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
//    friend struct is_del;
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
  
    is_del(MESH& mesh)
    {
	raw=thrust::raw_pointer_cast(mesh.stencil.data());
//	stencil.resize(10);
//	thrust::sequence(stencil.begin(),stencil.end(),0,1);
    }
  
    __device__
    bool operator()(const LINKLIST_MEMBER& x)
    {
//	return (x.Meshless_Ind % 2) == 0;
//	return (raw[x.Meshless_Ind]) % 2 == 0;
        return raw[x.Meshless_Ind] == 4;
    }
};

  

class LINKLIST
{
public:
    thrust::device_vector<LINKLIST_MEMBER> LINKLISTPOINT;
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
    
    void Delete_point(MESH& mesh)
    {
//	thrust::remove_if(LINKLISTPOINT.begin(), LINKLISTPOINT.end(), is_del());
//	is_del ISDEL;
        LINKLISTPOINT.erase(thrust::remove_if(LINKLISTPOINT.begin(),
                                              LINKLISTPOINT.end(),
					      is_del(mesh)), LINKLISTPOINT.end());
    }

};





int main()
{
      int i=3, j=6, k=7;
      
      MESH mesh;
      LINKLIST linklist;
      linklist.Insert_point(i);
      linklist.Insert_point(j);
      linklist.Insert_point(k);
      
      linklist.Delete_point(mesh);
      
      std::cout<<linklist.LINKLISTPOINT.size()<<std::endl;
      std::cout<<" it contains: "<< std::endl;
      for(thrust::device_vector<LINKLIST_MEMBER>::iterator iter = linklist.LINKLISTPOINT.begin(); iter != linklist.LINKLISTPOINT.end(); iter++)  {
                std::cout << (static_cast<LINKLIST_MEMBER>(*iter)).Meshless_Ind << std::endl;
        }

}