#include "../include/linklist.h"
#include <thrust/remove.h>
#include <thrust/for_each.h>
#include <iostream>
#include <stdio.h>

LINKLIST :: LINKLIST()
{
    std::cout << "Linklist is being initialized......";
    
    std::cout << "Initialization done." << std::endl;
}

LINKLIST :: ~LINKLIST()
{
    std::cout << "Linklist is being deleted." << std::endl;
}

void LINKLIST :: Insert_point(int i)
{
    /*maybe add filter later*/
    LINKLIST_MEMBER temp;
    temp.Meshless_Ind = i;
    LINKLISTPOINT.push_back(temp);
}

struct is_del
{
    int *TYPE1_raw, *TYPE2_raw, *TYPET_raw;
    double *U_raw, *V_raw, *W_raw, *P_raw;
    
    is_del(MESH_CARTESIAN *cartesian, FLOW_FIELD& flow_field)
    {
        TYPE1_raw = thrust::raw_pointer_cast( cartesian->POINTTYPE[0].data() );
        TYPE2_raw = thrust::raw_pointer_cast( cartesian->POINTTYPE[1].data() );
        TYPET_raw = thrust::raw_pointer_cast( cartesian->POINTTYPE[2].data() );
	
        U_raw = thrust::raw_pointer_cast( &flow_field.U.row(0)[0] );
        V_raw = thrust::raw_pointer_cast( &flow_field.V.row(0)[0] );
        W_raw = thrust::raw_pointer_cast( &flow_field.W.row(0)[0] );
        P_raw = thrust::raw_pointer_cast( &flow_field.P.row(0)[0] );
    }
    
    __device__
    bool operator()(const LINKLIST_MEMBER& x)
    {
        int i = x.Meshless_Ind;
	
	bool s = ( TYPE2_raw[i] == 0 || TYPE2_raw[i] == 1 ) && (TYPET_raw[i] == 3);
	if (s)
	{
	    printf("Delete point %d activated.\n", i);
	    if (TYPE2_raw[i] == 0)
	    {
		U_raw[i] = 0.0; V_raw[i] = 0.0; W_raw[i] = 0.0; P_raw[i] = 0.0;
	    }
	}
	
        return s;
    }
};
/*
void LINKLIST :: Delete_point(MESH_CARTESIAN& cartesian, FLOW_FIELD& flow_field)
{
    int i;
    for(thrust::device_vector<LINKLIST_MEMBER>::iterator iter = LINKLISTPOINT.begin(); iter != LINKLISTPOINT.end(); iter++)  
    {
	i = (static_cast<LINKLIST_MEMBER>(*iter)).Meshless_Ind;
	
	if( cartesian.POINTTYPE[2][i]==3 && (cartesian.POINTTYPE[1][i]==0||cartesian.POINTTYPE[1][i]==1) ) std::cout<<"Point "<<i<<" is needed to delete. IN LINKLIST"<<std::endl;
    }
  
  
    std::cout<<"Before Deleting point, the length is "<<LINKLISTPOINT.size()<<" ";
    LINKLISTPOINT.erase( thrust::remove_if( LINKLISTPOINT.begin(),
                                            LINKLISTPOINT.end(),
                                            is_del( cartesian, flow_field ) ), LINKLISTPOINT.end() );
    std::cout<<"After deleting point, the length is "<<LINKLISTPOINT.size()<<std::endl;
}
*/


void LINKLIST :: Delete_point(MESH& mesh, FLOW_FIELD& flow_field)
{
    std::cout<<"Before Deleting point, the length is "<<LINKLISTPOINT.size()<<" ";
    LINKLISTPOINT.erase( thrust::remove_if( LINKLISTPOINT.begin(),
                                            LINKLISTPOINT.end(),
                                            is_del( mesh.cartesian, flow_field ) ), LINKLISTPOINT.end() );
    std::cout<<"After deleting point, the length is "<<LINKLISTPOINT.size()<<std::endl;
}

struct is_fresh
{
    double *U_raw, *V_raw, *W_raw, *P_raw, *U_old_raw, *V_old_raw, *W_old_raw, *P_old_raw;
    int    *TYPE2_raw;
    MESH::UNIFIED_POSITION XYZ_raw;
    
    is_fresh(MESH& mesh, MESH_CARTESIAN& cartesian, FLOW_FIELD& flow_field)
    {        
        XYZ_raw = mesh.XYZ;
        
        U_raw = thrust::raw_pointer_cast( &flow_field.U.row(0)[0] );
        V_raw = thrust::raw_pointer_cast( &flow_field.V.row(0)[0] );
        W_raw = thrust::raw_pointer_cast( &flow_field.W.row(0)[0] );
        P_raw = thrust::raw_pointer_cast( &flow_field.P.row(0)[0] );
        
        U_old_raw = thrust::raw_pointer_cast( &flow_field.U_old.row(0)[0] );
        V_old_raw = thrust::raw_pointer_cast( &flow_field.V_old.row(0)[0] );
        W_old_raw = thrust::raw_pointer_cast( &flow_field.W_old.row(0)[0] );
        P_old_raw = thrust::raw_pointer_cast( &flow_field.P_old.row(0)[0] );
	
        TYPE2_raw = thrust::raw_pointer_cast( cartesian.POINTTYPE[1].data() );
    }
    
    __device__
    void operator()( LINKLIST_MEMBER& s )
    {
        int ijk = s.Meshless_Ind, ijk_p, i, io, j;
        
        if( TYPE2_raw[ijk] == 4 )
        {
	    double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
	    double tmp1_u, tmp2_u, tmp3_u, tmp4_u, tmp5_u, tmp6_u, tmp7_u, tmp8_u, tmp9_u;
	    double tmp1_v, tmp2_v, tmp3_v, tmp4_v, tmp5_v, tmp6_v, tmp7_v, tmp8_v, tmp9_v;
	    double tmp1_w, tmp2_w, tmp3_w, tmp4_w, tmp5_w, tmp6_w, tmp7_w, tmp8_w, tmp9_w;
	    double tmp1_p, tmp2_p, tmp3_p, tmp4_p, tmp5_p, tmp6_p, tmp7_p, tmp8_p, tmp9_p;
	    double drp[2], dx[2], dy[2], dz[2], tmpu[2], tmpv[2], tmpw[2], tmpp[2]; 
	    double dr;
	  
	    for(i = 0; i < 2; i++)
	    {
		ijk_p = s.Nb_Points[i];

		dx[i] = XYZ_raw(0, ijk_p) - XYZ_raw(0, ijk);
		dy[i] = XYZ_raw(1, ijk_p) - XYZ_raw(1, ijk);
		dz[i] = XYZ_raw(2, ijk_p) - XYZ_raw(2, ijk);
		drp[i] = sqrt(dx[i]*dx[i]+dy[i]*dy[i]+dz[i]*dz[i]);

		tmp1_u=0;tmp2_u=0;tmp3_u=0;tmp4_u=0;tmp5_u=0;tmp6_u=0;tmp7_u=0;tmp8_u=0;tmp9_u=0;
		tmp1_v=0;tmp2_v=0;tmp3_v=0;tmp4_v=0;tmp5_v=0;tmp6_v=0;tmp7_v=0;tmp8_v=0;tmp9_v=0;
		tmp1_w=0;tmp2_w=0;tmp3_w=0;tmp4_w=0;tmp5_w=0;tmp6_w=0;tmp7_w=0;tmp8_w=0;tmp9_w=0;
		tmp1_p=0;tmp2_p=0;tmp3_p=0;tmp4_p=0;tmp5_p=0;tmp6_p=0;tmp7_p=0;tmp8_p=0;tmp9_p=0;
		tmp1=0;tmp2=0;tmp3=0;tmp4=0;tmp5=0;tmp6=0;tmp7=0;tmp8=0;tmp9=0;
		
		for(io = 0; io < NB; io++)
		{
		    j = s.Nb_Points[io];
		    
		    tmp1_u = tmp1_u + s.Csvd[0][io]*U_raw[j];
		    tmp2_u = tmp2_u + s.Csvd[1][io]*U_raw[j];
		    tmp3_u = tmp3_u + s.Csvd[2][io]*U_raw[j];
		    tmp4_u = tmp4_u + s.Csvd[3][io]*U_raw[j];
		    tmp5_u = tmp5_u + s.Csvd[4][io]*U_raw[j];
		    tmp6_u = tmp6_u + s.Csvd[5][io]*U_raw[j];
		    tmp7_u = tmp7_u + s.Csvd[6][io]*U_raw[j];
		    tmp8_u = tmp8_u + s.Csvd[7][io]*U_raw[j];
		    tmp9_u = tmp9_u + s.Csvd[8][io]*U_raw[j];
		    
		    tmp1_v = tmp1_v + s.Csvd[0][io]*V_raw[j];
		    tmp2_v = tmp2_v + s.Csvd[1][io]*V_raw[j];
		    tmp3_v = tmp3_v + s.Csvd[2][io]*V_raw[j];
		    tmp4_v = tmp4_v + s.Csvd[3][io]*V_raw[j];
		    tmp5_v = tmp5_v + s.Csvd[4][io]*V_raw[j];
		    tmp6_v = tmp6_v + s.Csvd[5][io]*V_raw[j];
		    tmp7_v = tmp7_v + s.Csvd[6][io]*V_raw[j];
		    tmp8_v = tmp8_v + s.Csvd[7][io]*V_raw[j];
		    tmp9_v = tmp9_v + s.Csvd[8][io]*V_raw[j];
		
		    tmp1_w = tmp1_w + s.Csvd[0][io]*W_raw[j];
		    tmp2_w = tmp2_w + s.Csvd[1][io]*W_raw[j];
		    tmp3_w = tmp3_w + s.Csvd[2][io]*W_raw[j];
		    tmp4_w = tmp4_w + s.Csvd[3][io]*W_raw[j];
		    tmp5_w = tmp5_w + s.Csvd[4][io]*W_raw[j];
		    tmp6_w = tmp6_w + s.Csvd[5][io]*W_raw[j];
		    tmp7_w = tmp7_w + s.Csvd[6][io]*W_raw[j];
		    tmp8_w = tmp8_w + s.Csvd[7][io]*W_raw[j];
		    tmp9_w = tmp9_w + s.Csvd[8][io]*W_raw[j];
		
		    tmp1_p = tmp1_p + s.Csvd[0][io]*P_raw[j];
		    tmp2_p = tmp2_p + s.Csvd[1][io]*P_raw[j];
		    tmp3_p = tmp3_p + s.Csvd[2][io]*P_raw[j];
		    tmp4_p = tmp4_p + s.Csvd[3][io]*P_raw[j];
		    tmp5_p = tmp5_p + s.Csvd[4][io]*P_raw[j];
		    tmp6_p = tmp6_p + s.Csvd[5][io]*P_raw[j];
		    tmp7_p = tmp7_p + s.Csvd[6][io]*P_raw[j];
		    tmp8_p = tmp8_p + s.Csvd[7][io]*P_raw[j];
		    tmp9_p = tmp9_p + s.Csvd[8][io]*P_raw[j];
					
		    tmp1+=s.Csvd[0][io];
		    tmp2+=s.Csvd[1][io];
		    tmp3+=s.Csvd[2][io];
		    tmp4+=s.Csvd[3][io];
		    tmp5+=s.Csvd[4][io];
		    tmp6+=s.Csvd[5][io];
		    tmp7+=s.Csvd[6][io];
		    tmp8+=s.Csvd[7][io];
		    tmp9+=s.Csvd[8][io];
		}
		
		
		dr=(1-dx[i]*tmp1-dy[i]*tmp2-dz[i]*tmp3-0.5*dx[i]*dx[i]*tmp4-0.5*dy[i]*dy[i]*tmp5-0.5*dz[i]*dz[i]*tmp6-dx[i]*dy[i]*tmp7-dx[i]*dz[i]*tmp8-dy[i]*dz[i]*tmp9); 
		tmpu[i]=(U_raw[ijk_p]-dx[i]*tmp1_u-dy[i]*tmp2_u-dz[i]*tmp3_u-0.5*dx[i]*dx[i]*tmp4_u
			-0.5*dy[i]*dy[i]*tmp5_u-0.5*dz[i]*dz[i]*tmp6_u-dx[i]*dy[i]*tmp7_u
			-dx[i]*dz[i]*tmp8_u-dy[i]*dz[i]*tmp9_u)/dr;
		tmpv[i]=(V_raw[ijk_p]-dx[i]*tmp1_v-dy[i]*tmp2_v-dz[i]*tmp3_v-0.5*dx[i]*dx[i]*tmp4_v
			-0.5*dy[i]*dy[i]*tmp5_v-0.5*dz[i]*dz[i]*tmp6_v-dx[i]*dy[i]*tmp7_v
			-dx[i]*dz[i]*tmp8_v-dy[i]*dz[i]*tmp9_v)/dr;
		tmpw[i]=(W_raw[ijk_p]-dx[i]*tmp1_w-dy[i]*tmp2_w-dz[i]*tmp3_w-0.5*dx[i]*dx[i]*tmp4_w
			-0.5*dy[i]*dy[i]*tmp5_w-0.5*dz[i]*dz[i]*tmp6_w-dx[i]*dy[i]*tmp7_w
			-dx[i]*dz[i]*tmp8_w-dy[i]*dz[i]*tmp9_w)/dr;
		tmpp[i]=(P_raw[ijk_p]-dx[i]*tmp1_p-dy[i]*tmp2_p-dz[i]*tmp3_p-0.5*dx[i]*dx[i]*tmp4_p
			-0.5*dy[i]*dy[i]*tmp5_p-0.5*dz[i]*dz[i]*tmp6_p-dx[i]*dy[i]*tmp7_p
			-dx[i]*dz[i]*tmp8_p-dy[i]*dz[i]*tmp9_p)/dr;
			

	    }
			
	    U_raw[ijk] = (drp[1]*tmpu[0]+drp[0]*tmpu[1])/(drp[0]+drp[1]);
	    V_raw[ijk] = (drp[1]*tmpv[0]+drp[0]*tmpv[1])/(drp[0]+drp[1]);
	    W_raw[ijk] = (drp[1]*tmpw[0]+drp[0]*tmpw[1])/(drp[0]+drp[1]);
	    P_raw[ijk] = (drp[1]*tmpp[0]+drp[0]*tmpp[1])/(drp[0]+drp[1]);
	    U_old_raw[ijk] = U_raw[ijk];
	    V_old_raw[ijk] = V_raw[ijk];
	    W_old_raw[ijk] = W_raw[ijk];
	    P_old_raw[ijk] = P_raw[ijk];
        }
    }
};

void LINKLIST :: Fresh_point(MESH& mesh, MESH_CARTESIAN& cartesian, FLOW_FIELD& flow_field)
{
    thrust::for_each( LINKLISTPOINT.begin(), LINKLISTPOINT.end(), is_fresh( mesh, cartesian, flow_field ) );
}

struct is_update
{
    double *U_raw, *V_raw, *W_raw, *P_raw, *U_old_raw, *V_old_raw, *W_old_raw, *P_old_raw;
    int    *TYPE2_raw, *TYPET_raw;
    MESH::UNIFIED_POSITION XYZ_raw;
    
    is_update(MESH& mesh, MESH_CARTESIAN& cartesian, FLOW_FIELD& flow_field)
    {
        XYZ_raw = mesh.XYZ;
        
        U_raw = thrust::raw_pointer_cast( &flow_field.U.row(0)[0] );
        V_raw = thrust::raw_pointer_cast( &flow_field.V.row(0)[0] );
        W_raw = thrust::raw_pointer_cast( &flow_field.W.row(0)[0] );
        P_raw = thrust::raw_pointer_cast( &flow_field.P.row(0)[0] );
        
        U_old_raw = thrust::raw_pointer_cast( &flow_field.U_old.row(0)[0] );
        V_old_raw = thrust::raw_pointer_cast( &flow_field.V_old.row(0)[0] );
        W_old_raw = thrust::raw_pointer_cast( &flow_field.W_old.row(0)[0] );
        P_old_raw = thrust::raw_pointer_cast( &flow_field.P_old.row(0)[0] );
        
        TYPE2_raw = thrust::raw_pointer_cast( cartesian.POINTTYPE[1].data() );
        TYPET_raw = thrust::raw_pointer_cast( cartesian.POINTTYPE[2].data() );
    }
    
    __device__
    void operator()( LINKLIST_MEMBER& s )
    {
        int ijk = s.Meshless_Ind, ijk_p, i, io, j;
        
        if( (TYPE2_raw[ijk] == 1 && TYPET_raw[ijk] == 0) ||
            (TYPE2_raw[ijk] == 3 && TYPET_raw[ijk] == 0) ||
            (TYPE2_raw[ijk] == 4 && TYPET_raw[ijk] == 3) )
        {
            double tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
            double tmp1_u, tmp2_u, tmp3_u, tmp4_u, tmp5_u, tmp6_u, tmp7_u, tmp8_u, tmp9_u;
            double tmp1_v, tmp2_v, tmp3_v, tmp4_v, tmp5_v, tmp6_v, tmp7_v, tmp8_v, tmp9_v;
            double tmp1_w, tmp2_w, tmp3_w, tmp4_w, tmp5_w, tmp6_w, tmp7_w, tmp8_w, tmp9_w;
            double tmp1_p, tmp2_p, tmp3_p, tmp4_p, tmp5_p, tmp6_p, tmp7_p, tmp8_p, tmp9_p;
            double drp[2], dx[2], dy[2], dz[2], tmpu[2], tmpv[2], tmpw[2], tmpp[2];
            double dr;
            
            for(i = 0; i < 2; i++)
            {
                ijk_p = s.Nb_Points[i];
                
                dx[i] = XYZ_raw(0, ijk_p) - XYZ_raw(0, ijk);
                dy[i] = XYZ_raw(1, ijk_p) - XYZ_raw(1, ijk);
                dz[i] = XYZ_raw(2, ijk_p) - XYZ_raw(2, ijk);
                drp[i] = sqrt(dx[i]*dx[i]+dy[i]*dy[i]+dz[i]*dz[i]);

                tmp1_u=0;tmp2_u=0;tmp3_u=0;tmp4_u=0;tmp5_u=0;tmp6_u=0;tmp7_u=0;tmp8_u=0;tmp9_u=0;
                tmp1_v=0;tmp2_v=0;tmp3_v=0;tmp4_v=0;tmp5_v=0;tmp6_v=0;tmp7_v=0;tmp8_v=0;tmp9_v=0;
                tmp1_w=0;tmp2_w=0;tmp3_w=0;tmp4_w=0;tmp5_w=0;tmp6_w=0;tmp7_w=0;tmp8_w=0;tmp9_w=0;
                tmp1_p=0;tmp2_p=0;tmp3_p=0;tmp4_p=0;tmp5_p=0;tmp6_p=0;tmp7_p=0;tmp8_p=0;tmp9_p=0;
                tmp1=0;tmp2=0;tmp3=0;tmp4=0;tmp5=0;tmp6=0;tmp7=0;tmp8=0;tmp9=0;
                
                for(io = 0; io < NB; io++)
                {
                    j = s.Nb_Points[io];
                    
                    tmp1_u = tmp1_u + s.Csvd[0][io]*U_raw[j];
                    tmp2_u = tmp2_u + s.Csvd[1][io]*U_raw[j];
                    tmp3_u = tmp3_u + s.Csvd[2][io]*U_raw[j];
                    tmp4_u = tmp4_u + s.Csvd[3][io]*U_raw[j];
                    tmp5_u = tmp5_u + s.Csvd[4][io]*U_raw[j];
                    tmp6_u = tmp6_u + s.Csvd[5][io]*U_raw[j];
                    tmp7_u = tmp7_u + s.Csvd[6][io]*U_raw[j];
                    tmp8_u = tmp8_u + s.Csvd[7][io]*U_raw[j];
                    tmp9_u = tmp9_u + s.Csvd[8][io]*U_raw[j];
                    
                    tmp1_v = tmp1_v + s.Csvd[0][io]*V_raw[j];
                    tmp2_v = tmp2_v + s.Csvd[1][io]*V_raw[j];
                    tmp3_v = tmp3_v + s.Csvd[2][io]*V_raw[j];
                    tmp4_v = tmp4_v + s.Csvd[3][io]*V_raw[j];
                    tmp5_v = tmp5_v + s.Csvd[4][io]*V_raw[j];
                    tmp6_v = tmp6_v + s.Csvd[5][io]*V_raw[j];
                    tmp7_v = tmp7_v + s.Csvd[6][io]*V_raw[j];
                    tmp8_v = tmp8_v + s.Csvd[7][io]*V_raw[j];
                    tmp9_v = tmp9_v + s.Csvd[8][io]*V_raw[j];
                    
                    tmp1_w = tmp1_w + s.Csvd[0][io]*W_raw[j];
                    tmp2_w = tmp2_w + s.Csvd[1][io]*W_raw[j];
                    tmp3_w = tmp3_w + s.Csvd[2][io]*W_raw[j];
                    tmp4_w = tmp4_w + s.Csvd[3][io]*W_raw[j];
                    tmp5_w = tmp5_w + s.Csvd[4][io]*W_raw[j];
                    tmp6_w = tmp6_w + s.Csvd[5][io]*W_raw[j];
                    tmp7_w = tmp7_w + s.Csvd[6][io]*W_raw[j];
                    tmp8_w = tmp8_w + s.Csvd[7][io]*W_raw[j];
                    tmp9_w = tmp9_w + s.Csvd[8][io]*W_raw[j];
                    
                    tmp1_p = tmp1_p + s.Csvd[0][io]*P_raw[j];
                    tmp2_p = tmp2_p + s.Csvd[1][io]*P_raw[j];
                    tmp3_p = tmp3_p + s.Csvd[2][io]*P_raw[j];
                    tmp4_p = tmp4_p + s.Csvd[3][io]*P_raw[j];
                    tmp5_p = tmp5_p + s.Csvd[4][io]*P_raw[j];
                    tmp6_p = tmp6_p + s.Csvd[5][io]*P_raw[j];
                    tmp7_p = tmp7_p + s.Csvd[6][io]*P_raw[j];
                    tmp8_p = tmp8_p + s.Csvd[7][io]*P_raw[j];
                    tmp9_p = tmp9_p + s.Csvd[8][io]*P_raw[j];
                    
                    tmp1+=s.Csvd[0][io];
                    tmp2+=s.Csvd[1][io];
                    tmp3+=s.Csvd[2][io];
                    tmp4+=s.Csvd[3][io];
                    tmp5+=s.Csvd[4][io];
                    tmp6+=s.Csvd[5][io];
                    tmp7+=s.Csvd[6][io];
                    tmp8+=s.Csvd[7][io];
                    tmp9+=s.Csvd[8][io];
                }
                
                
                dr=(1-dx[i]*tmp1-dy[i]*tmp2-dz[i]*tmp3-0.5*dx[i]*dx[i]*tmp4-0.5*dy[i]*dy[i]*tmp5-0.5*dz[i]*dz[i]*tmp6-dx[i]*dy[i]*tmp7-dx[i]*dz[i]*tmp8-dy[i]*dz[i]*tmp9);
                tmpu[i]=(U_raw[ijk_p]-dx[i]*tmp1_u-dy[i]*tmp2_u-dz[i]*tmp3_u-0.5*dx[i]*dx[i]*tmp4_u
                         -0.5*dy[i]*dy[i]*tmp5_u-0.5*dz[i]*dz[i]*tmp6_u-dx[i]*dy[i]*tmp7_u
                         -dx[i]*dz[i]*tmp8_u-dy[i]*dz[i]*tmp9_u)/dr;
                tmpv[i]=(V_raw[ijk_p]-dx[i]*tmp1_v-dy[i]*tmp2_v-dz[i]*tmp3_v-0.5*dx[i]*dx[i]*tmp4_v
                         -0.5*dy[i]*dy[i]*tmp5_v-0.5*dz[i]*dz[i]*tmp6_v-dx[i]*dy[i]*tmp7_v
                         -dx[i]*dz[i]*tmp8_v-dy[i]*dz[i]*tmp9_v)/dr;
                tmpw[i]=(W_raw[ijk_p]-dx[i]*tmp1_w-dy[i]*tmp2_w-dz[i]*tmp3_w-0.5*dx[i]*dx[i]*tmp4_w
                         -0.5*dy[i]*dy[i]*tmp5_w-0.5*dz[i]*dz[i]*tmp6_w-dx[i]*dy[i]*tmp7_w
                         -dx[i]*dz[i]*tmp8_w-dy[i]*dz[i]*tmp9_w)/dr;
                tmpp[i]=(P_raw[ijk_p]-dx[i]*tmp1_p-dy[i]*tmp2_p-dz[i]*tmp3_p-0.5*dx[i]*dx[i]*tmp4_p
                         -0.5*dy[i]*dy[i]*tmp5_p-0.5*dz[i]*dz[i]*tmp6_p-dx[i]*dy[i]*tmp7_p
                         -dx[i]*dz[i]*tmp8_p-dy[i]*dz[i]*tmp9_p)/dr;
                
            }
            
            U_raw[ijk] = (drp[1]*tmpu[0]+drp[0]*tmpu[1])/(drp[0]+drp[1]);
            V_raw[ijk] = (drp[1]*tmpv[0]+drp[0]*tmpv[1])/(drp[0]+drp[1]);
            W_raw[ijk] = (drp[1]*tmpw[0]+drp[0]*tmpw[1])/(drp[0]+drp[1]);
            P_raw[ijk] = (drp[1]*tmpp[0]+drp[0]*tmpp[1])/(drp[0]+drp[1]);
            U_old_raw[ijk] = U_raw[ijk];
            V_old_raw[ijk] = V_raw[ijk];
            W_old_raw[ijk] = W_raw[ijk];
            P_old_raw[ijk] = P_raw[ijk];
        }
    }
};

void LINKLIST :: Update_point(MESH& mesh, MESH_CARTESIAN& cartesian, FLOW_FIELD& flow_field)
{
    thrust::for_each( LINKLISTPOINT.begin(), LINKLISTPOINT.end(), is_update( mesh, cartesian, flow_field ) );
}