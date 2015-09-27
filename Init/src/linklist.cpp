
typedef struct LINKLIST_MEMBER;
{
    int Meshless_Ind;
    int Nb_Points[NB];
    int Nb_Points_Old[NB];
    double Csvd[9][NB];
    double Csvd_Old[9][NB];
}LINKLIST_MEMBER;


template<typename MemorySpace> class LINKLIST
{
public:
    cusp::array1d<LINKLIST_MEMBER,MemorySpace> LINKLISTPOINT;
    
    LINKLIST();
    ~LINKLIST();
    Insert_point(int i);
    
};

template<typename MemorySpace>
LINKLIST<MemorySpace> :: LINKLIST()
{
    std::cout << "Linklist is being initialized......";
    
    std::cout << "Initialization done." << std::endl;
}

template<typename MemorySpace>
LINKLIST<MemorySpace> :: ~LINKLIST()
{
    std::cout << "Linklist is being deleted." << std::endl;
}

template<typename MemorySpace>
void LINKLIST<MemorySpace> :: Insert_point(int i)
{
    /*maybe add filter later*/
    LINKLIST_MEMBER temp;
    temp.Meshless_Ind = i;
    LINKLISTPOINT.push_back(temp);
}

struct is_del
{
    __host__ __device__
    bool operator()(const Meshless_Member& x)
    {
        int i=x.Meshless_Ind;
        return ((GType2_DEV[i]==0||GType2_DEV[i]==1)&&GTypeT_DEV[i]==3);
    }
};

template<typename MemorySpace>
void LINKLIST<MemorySpace> :: Delete_point()
{
    LINKLISTPOINT.erase(thrust::remove_if(LINKLISTPOINT.begin(),
                                          LINKLISTPOINT.end(),
                                          is_del()), LINKLISTPOINT.end());
    /*or do remove_if first, then for_each make UVWP to be 0, then delete unnecessary element.*/
}






#include "../common/svd.h"

vector<int> listindex;
int FILTER=0;

extern vector<Meshless_Member> Meshless_Member_List;
extern short  *GType2, *GTypeT;
extern double *U, *V, *W, *P, *U_Old, *V_Old, *W_Old, *P_Old, *X, *Y, *Z;

void Copy_New2Old_C2M()
{
    int i,io,j; 
    #pragma omp parallel for private(io,j)
    for(i=0;i<Meshless_Member_List.size();i++)
    {
	for(io=0;io<NB;io++)
	{
	    Meshless_Member_List[i].Nb_Points_Old[io]=Meshless_Member_List[i].Nb_Points[io];
	    for(j=0;j<6;j++) Meshless_Member_List[i].Csvd_Old[j][io]=Meshless_Member_List[i].Csvd[j][io];
	}
    }
}

void Insert_CarPoint(int i)
{
    Meshless_Member temp;
    if(Meshless_Member_List.size()==0)
    {
	temp.Meshless_Ind = i;
	Meshless_Member_List.push_back(temp);
    }
    else
    {
	if(FILTER==1)
	{
	    for(int s=0;s<Meshless_Member_List.size();s++)
	    {
		if(Meshless_Member_List[s].Meshless_Ind==i) 
		{
			cout<<"Link list point repeat."<<endl;
			return;
		}
	    }
	    temp.Meshless_Ind=i;
	    Meshless_Member_List.push_back(temp);
	}
	else
	{
	    temp.Meshless_Ind=i;
	    Meshless_Member_List.push_back(temp);
	}
    }
}

void Del_OldFSI()
{ 
    int i;
    vector<Meshless_Member>::iterator itr=Meshless_Member_List.begin();
    while(itr != Meshless_Member_List.end())
    {
	i=itr->Meshless_Ind;
	if((GType2[i]==0||GType2[i]==1)&&GTypeT[i]==3)
	{
	    if(GType2[i]==0)
	    {
		U[i]=0;
		V[i]=0;
		W[i]=0;
		P[i]=0;
	    }
	    itr=Meshless_Member_List.erase(itr);
	}
	else
	{
	    itr++;
	}
    }
}

void Fresh_Node()
{
    Meshless_Member *temp;
    int i, j, ijk, ijk_p, io, s;
    double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
    double tmp1_u, tmp2_u, tmp3_u, tmp4_u, tmp5_u, tmp6_u, tmp7_u, tmp8_u, tmp9_u;
    double tmp1_v, tmp2_v, tmp3_v, tmp4_v, tmp5_v, tmp6_v, tmp7_v, tmp8_v, tmp9_v;
    double tmp1_w, tmp2_w, tmp3_w, tmp4_w, tmp5_w, tmp6_w, tmp7_w, tmp8_w, tmp9_w;
    double tmp1_p, tmp2_p, tmp3_p, tmp4_p, tmp5_p, tmp6_p, tmp7_p, tmp8_p, tmp9_p;
    double drp[2], dx[2], dy[2], dz[2], tmpu[2], tmpv[2], tmpw[2], tmpp[2]; 
    double dr;
    
    for(s=0;s<Meshless_Member_List.size();s++)
    {
	temp=&Meshless_Member_List[s];
	ijk=temp->Meshless_Ind;
	if(GType2[ijk]==4)
	{
	    for(i=0;i<2;i++)
	    {
		ijk_p=temp->Nb_Points[i];
		dx[i]=X[ijk_p]-X[ijk];
		dy[i]=Y[ijk_p]-Y[ijk];
		dz[i]=Z[ijk_p]-Z[ijk];
		drp[i] = sqrt(dx[i]*dx[i]+dy[i]*dy[i]+dz[i]*dz[i]);
		tmp1_u=0;tmp2_u=0;tmp3_u=0;tmp4_u=0;tmp5_u=0;tmp6_u=0;tmp7_u=0;tmp8_u=0;tmp9_u=0;
		tmp1_v=0;tmp2_v=0;tmp3_v=0;tmp4_v=0;tmp5_v=0;tmp6_v=0;tmp7_v=0;tmp8_v=0;tmp9_v=0;
		tmp1_w=0;tmp2_w=0;tmp3_w=0;tmp4_w=0;tmp5_w=0;tmp6_w=0;tmp7_w=0;tmp8_w=0;tmp9_w=0;
		tmp1_p=0;tmp2_p=0;tmp3_p=0;tmp4_p=0;tmp5_p=0;tmp6_p=0;tmp7_p=0;tmp8_p=0;tmp9_p=0;
		tmp1=0;tmp2=0;tmp3=0;tmp4=0;tmp5=0;tmp6=0;tmp7=0;tmp8=0;tmp9=0;
		for(io=0;io<NB;io++)
		{
		    j=temp->Nb_Points[io];
			       
		    tmp1_u=tmp1_u+temp->Csvd[0][io]*U[j];
		    tmp2_u=tmp2_u+temp->Csvd[1][io]*U[j];
		    tmp3_u=tmp3_u+temp->Csvd[2][io]*U[j];
		    tmp4_u=tmp4_u+temp->Csvd[3][io]*U[j];
		    tmp5_u=tmp5_u+temp->Csvd[4][io]*U[j];
		    tmp6_u=tmp6_u+temp->Csvd[5][io]*U[j];
		    tmp7_u=tmp7_u+temp->Csvd[6][io]*U[j];
		    tmp8_u=tmp8_u+temp->Csvd[7][io]*U[j];
		    tmp9_u=tmp9_u+temp->Csvd[8][io]*U[j];
		
		    tmp1_v=tmp1_v+temp->Csvd[0][io]*V[j];
		    tmp2_v=tmp2_v+temp->Csvd[1][io]*V[j];
		    tmp3_v=tmp3_v+temp->Csvd[2][io]*V[j];
		    tmp4_v=tmp4_v+temp->Csvd[3][io]*V[j];
		    tmp5_v=tmp5_v+temp->Csvd[4][io]*V[j];
		    tmp6_v=tmp6_v+temp->Csvd[5][io]*V[j];
		    tmp7_v=tmp7_v+temp->Csvd[6][io]*V[j];
		    tmp8_v=tmp8_v+temp->Csvd[7][io]*V[j];
		    tmp9_v=tmp9_v+temp->Csvd[8][io]*V[j];
		
		    tmp1_w=tmp1_w+temp->Csvd[0][io]*W[j];
		    tmp2_w=tmp2_w+temp->Csvd[1][io]*W[j];
		    tmp3_w=tmp3_w+temp->Csvd[2][io]*W[j];
		    tmp4_w=tmp4_w+temp->Csvd[3][io]*W[j];
		    tmp5_w=tmp5_w+temp->Csvd[4][io]*W[j];
		    tmp6_w=tmp6_w+temp->Csvd[5][io]*W[j];
		    tmp7_w=tmp7_w+temp->Csvd[6][io]*W[j];
		    tmp8_w=tmp8_w+temp->Csvd[7][io]*W[j];
		    tmp9_w=tmp9_w+temp->Csvd[8][io]*W[j];
		
		    tmp1_p=tmp1_p+temp->Csvd[0][io]*P[j];
		    tmp2_p=tmp2_p+temp->Csvd[1][io]*P[j];
		    tmp3_p=tmp3_p+temp->Csvd[2][io]*P[j];
		    tmp4_p=tmp4_p+temp->Csvd[3][io]*P[j];
		    tmp5_p=tmp5_p+temp->Csvd[4][io]*P[j];
		    tmp6_p=tmp6_p+temp->Csvd[5][io]*P[j];
		    tmp7_p=tmp7_p+temp->Csvd[6][io]*P[j];
		    tmp8_p=tmp8_p+temp->Csvd[7][io]*P[j];
		    tmp9_p=tmp9_p+temp->Csvd[8][io]*P[j];
					
		    tmp1+=temp->Csvd[0][io];
		    tmp2+=temp->Csvd[1][io];
		    tmp3+=temp->Csvd[2][io];
		    tmp4+=temp->Csvd[3][io];
		    tmp5+=temp->Csvd[4][io];
		    tmp6+=temp->Csvd[5][io];
		    tmp7+=temp->Csvd[6][io];
		    tmp8+=temp->Csvd[7][io];
		    tmp9+=temp->Csvd[8][io];
		}
		
		
		dr=(1-dx[i]*tmp1-dy[i]*tmp2-dz[i]*tmp3-0.5*dx[i]*dx[i]*tmp4-0.5*dy[i]*dy[i]*tmp5-0.5*dz[i]*dz[i]*tmp6-dx[i]*dy[i]*tmp7-dx[i]*dz[i]*tmp8-dy[i]*dz[i]*tmp9); 
		tmpu[i]=(U[ijk_p]-dx[i]*tmp1_u-dy[i]*tmp2_u-dz[i]*tmp3_u-0.5*dx[i]*dx[i]*tmp4_u
			-0.5*dy[i]*dy[i]*tmp5_u-0.5*dz[i]*dz[i]*tmp6_u-dx[i]*dy[i]*tmp7_u
			-dx[i]*dz[i]*tmp8_u-dy[i]*dz[i]*tmp9_u)/dr;
		tmpv[i]=(V[ijk_p]-dx[i]*tmp1_v-dy[i]*tmp2_v-dz[i]*tmp3_v-0.5*dx[i]*dx[i]*tmp4_v
			-0.5*dy[i]*dy[i]*tmp5_v-0.5*dz[i]*dz[i]*tmp6_v-dx[i]*dy[i]*tmp7_v
			-dx[i]*dz[i]*tmp8_v-dy[i]*dz[i]*tmp9_v)/dr;
		tmpw[i]=(W[ijk_p]-dx[i]*tmp1_w-dy[i]*tmp2_w-dz[i]*tmp3_w-0.5*dx[i]*dx[i]*tmp4_w
			-0.5*dy[i]*dy[i]*tmp5_w-0.5*dz[i]*dz[i]*tmp6_w-dx[i]*dy[i]*tmp7_w
			-dx[i]*dz[i]*tmp8_w-dy[i]*dz[i]*tmp9_w)/dr;
		tmpp[i]=(P[ijk_p]-dx[i]*tmp1_p-dy[i]*tmp2_p-dz[i]*tmp3_p-0.5*dx[i]*dx[i]*tmp4_p
			-0.5*dy[i]*dy[i]*tmp5_p-0.5*dz[i]*dz[i]*tmp6_p-dx[i]*dy[i]*tmp7_p
			-dx[i]*dz[i]*tmp8_p-dy[i]*dz[i]*tmp9_p)/dr;
			

	    }
			
	    U[ijk]=(drp[1]*tmpu[0]+drp[0]*tmpu[1])/(drp[0]+drp[1]);
	    V[ijk]=(drp[1]*tmpv[0]+drp[0]*tmpv[1])/(drp[0]+drp[1]);
	    W[ijk]=(drp[1]*tmpw[0]+drp[0]*tmpw[1])/(drp[0]+drp[1]);
	    P[ijk]=(drp[1]*tmpp[0]+drp[0]*tmpp[1])/(drp[0]+drp[1]);
	    U_Old[ijk]=U[ijk];
	    V_Old[ijk]=V[ijk];
	    W_Old[ijk]=W[ijk];
	    P_Old[ijk]=P[ijk];
	    //cout<<"fresh Point "<<ijk<<" :"<<U[ijk]<<" "<<V[ijk]<<" "<<W[ijk]<<" "<<P[ijk]<<endl;
	}
    }
}

void Update_UT()
{
    Meshless_Member *temp;
    int i,j,ijk,ijk_p,io,s;
    double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9;
    double tmp1_u,tmp2_u,tmp3_u,tmp4_u,tmp5_u,tmp6_u,tmp7_u,tmp8_u,tmp9_u;
    double tmp1_v,tmp2_v,tmp3_v,tmp4_v,tmp5_v,tmp6_v,tmp7_v,tmp8_v,tmp9_v;
    double tmp1_w,tmp2_w,tmp3_w,tmp4_w,tmp5_w,tmp6_w,tmp7_w,tmp8_w,tmp9_w;
    double tmp1_p,tmp2_p,tmp3_p,tmp4_p,tmp5_p,tmp6_p,tmp7_p,tmp8_p,tmp9_p;
    double drp[2],dx[2],dy[2],dz[2],tmpu[2],tmpv[2],tmpw[2],tmpp[2];
    double dr;
  
    for(s=0;s<Meshless_Member_List.size();s++)
    {
	temp=&Meshless_Member_List[s];
	ijk=temp->Meshless_Ind;
	if((GType2[ijk]==1&&GTypeT[ijk]==0)||(GType2[ijk]==3&&GTypeT[ijk]==0)||(GType2[ijk]==4&&GTypeT[ijk]==3))
	{
	    for(i=0;i<2;i++)
	    {
		ijk_p=temp->Nb_Points[i];
		dx[i]=X[ijk_p]-X[ijk];
		dy[i]=Y[ijk_p]-Y[ijk];
		dz[i]=Z[ijk_p]-Z[ijk];
		drp[i]=sqrt(dx[i]*dx[i]+dy[i]*dy[i]+dz[i]*dz[i]);
		tmp1_u=0;tmp2_u=0;tmp3_u=0;tmp4_u=0;tmp5_u=0;tmp6_u=0;tmp7_u=0;tmp8_u=0;tmp9_u=0;
		tmp1_v=0;tmp2_v=0;tmp3_v=0;tmp4_v=0;tmp5_v=0;tmp6_v=0;tmp7_v=0;tmp8_v=0;tmp9_v=0;
		tmp1_w=0;tmp2_w=0;tmp3_w=0;tmp4_w=0;tmp5_w=0;tmp6_w=0;tmp7_w=0;tmp8_w=0;tmp9_w=0;
		tmp1_p=0;tmp2_p=0;tmp3_p=0;tmp4_p=0;tmp5_p=0;tmp6_p=0;tmp7_p=0;tmp8_p=0;tmp9_p=0;
		tmp1=0;tmp2=0;tmp3=0;tmp4=0;tmp5=0;tmp6=0;tmp7=0;tmp8=0;tmp9=0;
		for(io=0;io<NB;io++)
		{
		  j=temp->Nb_Points[io];
		  
		  tmp1_u=tmp1_u+temp->Csvd[0][io]*U[j];
		  tmp2_u=tmp2_u+temp->Csvd[1][io]*U[j];
		  tmp3_u=tmp3_u+temp->Csvd[2][io]*U[j];
		  tmp4_u=tmp4_u+temp->Csvd[3][io]*U[j];
		  tmp5_u=tmp5_u+temp->Csvd[4][io]*U[j];
		  tmp6_u=tmp6_u+temp->Csvd[5][io]*U[j];
		  tmp7_u=tmp7_u+temp->Csvd[6][io]*U[j];
		  tmp8_u=tmp8_u+temp->Csvd[7][io]*U[j];
		  tmp9_u=tmp9_u+temp->Csvd[8][io]*U[j];
		  
		  tmp1_v=tmp1_v+temp->Csvd[0][io]*V[j];
		  tmp2_v=tmp2_v+temp->Csvd[1][io]*V[j];
		  tmp3_v=tmp3_v+temp->Csvd[2][io]*V[j];
		  tmp4_v=tmp4_v+temp->Csvd[3][io]*V[j];
		  tmp5_v=tmp5_v+temp->Csvd[4][io]*V[j];
		  tmp6_v=tmp6_v+temp->Csvd[5][io]*V[j];
		  tmp7_v=tmp7_v+temp->Csvd[6][io]*V[j];
		  tmp8_v=tmp8_v+temp->Csvd[7][io]*V[j];
		  tmp9_v=tmp9_v+temp->Csvd[8][io]*V[j];
		  
		  tmp1_w=tmp1_w+temp->Csvd[0][io]*W[j];
		  tmp2_w=tmp2_w+temp->Csvd[1][io]*W[j];
		  tmp3_w=tmp3_w+temp->Csvd[2][io]*W[j];
		  tmp4_w=tmp4_w+temp->Csvd[3][io]*W[j];
		  tmp5_w=tmp5_w+temp->Csvd[4][io]*W[j];
		  tmp6_w=tmp6_w+temp->Csvd[5][io]*W[j];
		  tmp7_w=tmp7_w+temp->Csvd[6][io]*W[j];
		  tmp8_w=tmp8_w+temp->Csvd[7][io]*W[j];
		  tmp9_w=tmp9_w+temp->Csvd[8][io]*W[j];
		  
		  tmp1_p=tmp1_p+temp->Csvd[0][io]*P[j];
		  tmp2_p=tmp2_p+temp->Csvd[1][io]*P[j];
		  tmp3_p=tmp3_p+temp->Csvd[2][io]*P[j];
		  tmp4_p=tmp4_p+temp->Csvd[3][io]*P[j];
		  tmp5_p=tmp5_p+temp->Csvd[4][io]*P[j];
		  tmp6_p=tmp6_p+temp->Csvd[5][io]*P[j];
		  tmp7_p=tmp7_p+temp->Csvd[6][io]*P[j];
		  tmp8_p=tmp8_p+temp->Csvd[7][io]*P[j];
		  tmp9_p=tmp9_p+temp->Csvd[8][io]*P[j];
		  
		  tmp1=tmp1+temp->Csvd[0][io];
		  tmp2=tmp2+temp->Csvd[1][io];
		  tmp3=tmp3+temp->Csvd[2][io];
		  tmp4=tmp4+temp->Csvd[3][io];
		  tmp5=tmp5+temp->Csvd[4][io];
		  tmp6=tmp6+temp->Csvd[5][io];
		  tmp7=tmp7+temp->Csvd[6][io];
		  tmp8=tmp8+temp->Csvd[7][io];
		  tmp9=tmp9+temp->Csvd[8][io];
		  
		}
		dr=(1-dx[i]*tmp1-dy[i]*tmp2-dz[i]*tmp3-0.5*dx[i]*dx[i]*tmp4
		-0.5*dy[i]*dy[i]*tmp5-0.5*dz[i]*dz[i]*tmp6-dx[i]*dy[i]*tmp7
		 -dx[i]*dz[i]*tmp8-dy[i]*dz[i]*tmp9);



		tmpu[i]=(U[ijk_p]-dx[i]*tmp1_u-dy[i]*tmp2_u-dz[i]*tmp3_u
			  -0.5*dx[i]*dx[i]*tmp4_u-0.5*dy[i]*dy[i]*tmp5_u
			  -0.5*dz[i]*dz[i]*tmp6_u-dx[i]*dy[i]*tmp7_u
			  -dx[i]*dz[i]*tmp8_u-dy[i]*dz[i]*tmp9_u)/dr;

			  
			  
		tmpv[i]=(V[ijk_p]-dx[i]*tmp1_v-dy[i]*tmp2_v-dz[i]*tmp3_v
			  -0.5*dx[i]*dx[i]*tmp4_v-0.5*dy[i]*dy[i]*tmp5_v
			  -0.5*dz[i]*dz[i]*tmp6_v-dx[i]*dy[i]*tmp7_v
			  -dx[i]*dz[i]*tmp8_v-dy[i]*dz[i]*tmp9_v)/dr;

			  
			  
		tmpw[i]=(W[ijk_p]-dx[i]*tmp1_w-dy[i]*tmp2_w-dz[i]*tmp3_w
			  -0.5*dx[i]*dx[i]*tmp4_w-0.5*dy[i]*dy[i]*tmp5_w
			  -0.5*dz[i]*dz[i]*tmp6_w-dx[i]*dy[i]*tmp7_w
			  -dx[i]*dz[i]*tmp8_w-dy[i]*dz[i]*tmp9_w)/dr;
			  

			  
		tmpp[i]=(P[ijk_p]-dx[i]*tmp1_p-dy[i]*tmp2_p-dz[i]*tmp3_p
			  -0.5*dx[i]*dx[i]*tmp4_p-0.5*dy[i]*dy[i]*tmp5_p
			  -0.5*dz[i]*dz[i]*tmp6_p-dx[i]*dy[i]*tmp7_p
			  -dx[i]*dz[i]*tmp8_p-dy[i]*dz[i]*tmp9_p)/dr;

		}
		//if(LPlot)cout<<"fresh point"<<ijk<<endl;
		U[ijk]=(drp[1]*tmpu[0]+drp[0]*tmpu[1])/(drp[0]+drp[1]);
		V[ijk]=(drp[1]*tmpv[0]+drp[0]*tmpv[1])/(drp[0]+drp[1]);
		W[ijk]=(drp[1]*tmpw[0]+drp[0]*tmpw[1])/(drp[0]+drp[1]);
		P[ijk]=(drp[1]*tmpp[0]+drp[0]*tmpp[1])/(drp[0]+drp[1]);
		U_Old[ijk]=U[ijk];
		V_Old[ijk]=V[ijk];
		W_Old[ijk]=W[ijk];
	    P_Old[ijk]=P[ijk];
	}
    }
}








