#include "../include/frame_structure.h"

#include <iostream>
#include <fstream>
#include <math.h>

#include <cusp/blas/blas.h>


double REF_FRAME :: compute_angular_velocity1(double Time)
{
    return -compute_theta(Time,1)*sin( compute_phi(Time,0) ) + compute_psi(Time,1)*cos( compute_theta(Time,0) )*cos( compute_phi(Time,0) ) ;
}

double REF_FRAME :: compute_angular_velocity2(double Time)
{
    return compute_theta(Time,1)*cos( compute_phi(Time,0) ) + compute_psi(Time,1)*cos( compute_theta(Time,0) )*sin( compute_phi(Time,0) ) ;
}

double REF_FRAME :: compute_angular_velocity3(double Time)
{
    return compute_phi(Time,1) - compute_psi(Time,1)*sin( compute_theta(Time,0) ) ;
}

double REF_FRAME :: compute_angular_acceleration1(double Time)
{
    return (  - compute_theta(Time,2) * sin( compute_phi(Time,0) )
    		  - compute_theta(Time,1) * compute_phi(Time,1) * cos( compute_phi(Time,0) )
              + compute_psi(Time,2) * cos( compute_theta(Time,0) ) * cos( compute_phi(Time,0) )
              - compute_psi(Time,1) * compute_theta(Time,1) * sin( compute_theta(Time,0) )*cos( compute_phi(Time,0) )
              - compute_psi(Time,1) * compute_phi(Time,1) * cos( compute_theta(Time,0) )*sin( compute_phi(Time,0))  );
}

double REF_FRAME :: compute_angular_acceleration2(double Time)
{
    return (  compute_theta(Time,2) * cos( compute_phi(Time,0) )
              - compute_theta(Time,1) * compute_phi(Time,1) * sin( compute_phi(Time,0) )
              + compute_psi(Time,2) * cos(compute_theta(Time,0)) * sin( compute_phi(Time,0) )
              - compute_psi(Time,1) * compute_theta(Time,1) * sin( compute_theta(Time,0) )*sin( compute_phi(Time,0) )
              + compute_psi(Time,1) * compute_phi(Time,1) * cos( compute_theta(Time,0) )*cos( compute_phi(Time,0) )  );
}

double REF_FRAME :: compute_angular_acceleration3(double Time)
{
    return (  compute_phi(Time,2)
              - compute_psi(Time,1) * compute_theta(Time,1) * cos( compute_theta(Time,0) )
              - compute_psi(Time,2) * sin( compute_theta(Time,0) )  );
}

void REF_FRAME :: update_all()
{
    position[0]=compute_x(Time,0); //std::cout<<"x "<<position[0]<<" ";
    position[1]=compute_y(Time,0); //std::cout<<"y "<<position[1]<<" ";
    position[2]=compute_z(Time,0); //std::cout<<"z "<<position[2]<<std::endl;
    velocity[0]=compute_x(Time,1);
    velocity[1]=compute_y(Time,1);
    velocity[2]=compute_z(Time,1);
    acceleration[0]=compute_x(Time,2);
    acceleration[1]=compute_y(Time,2);
    acceleration[2]=compute_z(Time,2);
    
    angle[0]=compute_phi(Time,0); //std::cout<<"phi = "<<angle[0]<<" ";
    angle[1]=compute_theta(Time,0); //std::cout<<"theta = "<<angle[1]<<" ";
    angle[2]=compute_psi(Time,0); //std::cout<<"psi = "<<angle[2]<<std::endl;
    angular_velocity[0]=compute_angular_velocity1(Time);
    angular_velocity[1]=compute_angular_velocity2(Time);
    angular_velocity[2]=compute_angular_velocity3(Time);
    compute_orientation();
    angular_acceleration[0]=compute_angular_acceleration1(Time);
    angular_acceleration[1]=compute_angular_acceleration2(Time);
    angular_acceleration[2]=compute_angular_acceleration3(Time);
}


void REF_FRAME :: copy_new2old()
{
    for(int i=0;i<3;i++)
    {
	position_old[i]             =    position[i];
	velocity_old[i]             =    velocity[i];
	acceleration_old[i]         =    acceleration[i];
	angle_old[i]                =    angle[i];
	angular_velocity_old[i]     =    angular_velocity[i];
	angular_acceleration_old[i] =    angular_acceleration[i];
	force_old[i]                =    force[i];
	torque_old[i]               =    torque[i];
	momentum_old[i]             =    momentum[i];
	angular_momentum_old[i]     =    angular_momentum[i];
	for(int j=0;j<3;j++)
	{
	    orientation_old[i][j]   =    orientation[i][j];
	    inertia_old[i][j]       =    inertia[i][j];
	}
    }
}

void REF_FRAME :: compute_orientation_tmp()
{
	double Phi,Theta,Psi;
	Phi=angle_tmp[0];Theta=angle_tmp[1];Psi=angle_tmp[2];
	orientation_tmp[0][0] = cos(Phi)*cos(Theta);
	orientation_tmp[0][1] = cos(Phi)*sin(Theta)*sin(Psi) - sin(Phi)*cos(Psi);
	orientation_tmp[0][2] = cos(Phi)*sin(Theta)*cos(Psi) + sin(Phi)*sin(Psi);
	orientation_tmp[1][0] = sin(Phi)*cos(Theta);
	orientation_tmp[1][1] = sin(Phi)*sin(Theta)*sin(Psi) + cos(Phi)*cos(Psi);
	orientation_tmp[1][2] = sin(Phi)*sin(Theta)*cos(Psi) - cos(Phi)*sin(Psi);
	orientation_tmp[2][0] = -sin(Theta);
	orientation_tmp[2][1] = cos(Theta)*sin(Psi);
	orientation_tmp[2][2] = cos(Theta)*cos(Psi);
}

void REF_FRAME :: compute_orientation()
{
	double Phi,Theta,Psi;
	Phi=angle[0];Theta=angle[1];Psi=angle[2];
	orientation[0][0] = cos(Phi)*cos(Theta);
	orientation[0][1] = cos(Phi)*sin(Theta)*sin(Psi) - sin(Phi)*cos(Psi);
	orientation[0][2] = cos(Phi)*sin(Theta)*cos(Psi) + sin(Phi)*sin(Psi);
	orientation[1][0] = sin(Phi)*cos(Theta);
	orientation[1][1] = sin(Phi)*sin(Theta)*sin(Psi) + cos(Phi)*cos(Psi);
	orientation[1][2] = sin(Phi)*sin(Theta)*cos(Psi) - cos(Phi)*sin(Psi);
	orientation[2][0] = -sin(Theta);
	orientation[2][1] = cos(Theta)*sin(Psi);
	orientation[2][2] = cos(Theta)*cos(Psi);
}

void REF_FRAME :: memory_allocate_for_rigid_body(RIGID_OBJ *temp)
{    
    temp->XYZ.resize(3,temp->POINT_NUMBER);
    temp->UVW.resize(3,temp->POINT_NUMBER);
    temp->ACC.resize(3,temp->POINT_NUMBER);
    temp->OUTER_NORMAL_VECTOR.resize(3,temp->POINT_NUMBER);
    
    temp->TRIANGLE.resize(3,temp->TRIANGLE_NUMBER);
    temp->TETRAHEDRON.resize(4,temp->TETRAHEDRON_NUMBER);

    temp->AREA.resize(temp->POINT_NUMBER);
    temp->INNERMARK.resize(temp->POINT_NUMBER);
    temp->OUTERMARK.resize(temp->POINT_NUMBER);
}

RIGID_OBJ* REF_FRAME :: get_rigid_body_from_file(const char *filename)
{
    RIGID_OBJ *temp = new RIGID_OBJ;
    std::ifstream in;
    in.open(filename);
    in>>temp->IDENTITY;
    in>>temp->LEVEL;
    in>>temp->POINTMODE;
    in>>temp->POINT_NUMBER>>temp->TETRAHEDRON_NUMBER>>temp->TRIANGLE_NUMBER;
    in>>temp->INNER_POINT_NUMBER>>temp->OUTER_POINT_NUMBER;

    memory_allocate_for_rigid_body(temp);

    double *readXYZ[3], *readONV[3], *readAREA;
    int    *readINNERMARK, *readOUTERMARK, *readTETRAHEDRON[4], *readTRIANGLE[3];
    
    for(int i = 0; i < 3; i++) readXYZ[i] = new double [temp->POINT_NUMBER];
    for(int i = 0; i < 3; i++) readONV[i] = new double [temp->POINT_NUMBER];
    readAREA = new double [temp->POINT_NUMBER];
    readINNERMARK = new int [temp->POINT_NUMBER];
    readOUTERMARK = new int [temp->POINT_NUMBER];
    for(int i = 0; i < 4; i++) readTETRAHEDRON[i] = new int [temp->TETRAHEDRON_NUMBER];
    for(int i = 0; i < 3; i++) readTRIANGLE[i] = new int [temp->TRIANGLE_NUMBER];
    
    for (int i = 0; i < temp->POINT_NUMBER; i++)
    {
	in >> readXYZ[0][i] >> readXYZ[1][i] >> readXYZ[2][i];
	in >> readONV[0][i] >> readONV[1][i] >> readONV[2][i];
	in >> readINNERMARK[i] >> readOUTERMARK[i];
	in >> readAREA[i];
	
//	std::cout<<"Read file, onv is "<<readONV[0][i]<<" "<<readONV[1][i]<<" "<<readONV[2][i]<<" "<<sqrt(readONV[0][i]*readONV[0][i]+readONV[1][i]*readONV[1][i]+readONV[2][i]*readONV[2][i])<<std::endl;
    }
    for (int j = 0; j < temp->TETRAHEDRON_NUMBER; j++)
    {
	in >> readTETRAHEDRON[0][j] >> readTETRAHEDRON[1][j] >> readTETRAHEDRON[2][j] >> readTETRAHEDRON[3][j];
    }
    for (int j = 0; j < temp->TRIANGLE_NUMBER; j++)
    {
	in >> readTRIANGLE[0][j] >> readTRIANGLE[1][j] >> readTRIANGLE[2][j];
    }
    in.close();
    
    for(int i = 0; i < 3; i++) thrust::copy(readXYZ[i], readXYZ[i] + temp->POINT_NUMBER, temp->XYZ.row(i).begin());
    for(int i = 0; i < 3; i++) thrust::fill(temp->UVW.row(i).begin(), temp->UVW.row(i).end(), 0.0);
    for(int i = 0; i < 3; i++) thrust::fill(temp->ACC.row(i).begin(), temp->ACC.row(i).end(), 0.0);
    for(int i = 0; i < 3; i++) thrust::copy(readONV[i], readONV[i] + temp->POINT_NUMBER, temp->OUTER_NORMAL_VECTOR.row(i).begin());
    thrust::copy(readINNERMARK, readINNERMARK + temp->POINT_NUMBER, temp->INNERMARK.begin());
    thrust::copy(readOUTERMARK, readOUTERMARK + temp->POINT_NUMBER, temp->OUTERMARK.begin());
    thrust::copy(readAREA, readAREA + temp->POINT_NUMBER, temp->AREA.begin());
    for(int i = 0; i < 4; i++) thrust::copy(readTETRAHEDRON[i], readTETRAHEDRON[i] + temp->TETRAHEDRON_NUMBER, temp->TETRAHEDRON.row(i).begin());
    for(int i = 0; i < 3; i++) thrust::copy(readTRIANGLE[i], readTRIANGLE[i] + temp->TRIANGLE_NUMBER, temp->TRIANGLE.row(i).begin());
    
    for(int i = 0; i < 3; i++) delete[] readXYZ[i];
    for(int i = 0; i < 3; i++) delete[] readONV[i];
    delete[] readAREA;
    delete[] readINNERMARK;
    delete[] readOUTERMARK;
    for(int i = 0; i < 4; i++) delete[] readTETRAHEDRON[i];
    for(int i = 0; i < 3; i++) delete[] readTRIANGLE[i];
    
    std::cout << "Point is: " << temp->POINT_NUMBER << std::endl;

    return temp;
}


RIGID_OBJ* REF_FRAME :: get_rigid_body_from_subframe(int no_fra,int no_obj)
{
    if(subframe_number==0){std::cout<<"There is no subframe!"<<std::endl;}
    if(no_fra<0||no_fra>=subframe_number){std::cout<<"The range of no_fra is wrong!"<<std::endl;}
    if(no_obj<0||no_obj>=sub[no_fra]->obj_number){std::cout<<"The range of no_obj is wrong!"<<std::endl;}
    RIGID_OBJ *temp=new RIGID_OBJ;
    temp->IDENTITY=sub[no_fra]->rigid_body[no_obj]->IDENTITY;
    temp->LEVEL=sub[no_fra]->rigid_body[no_obj]->LEVEL-1;
    temp->POINTMODE=sub[no_fra]->rigid_body[no_obj]->POINTMODE;
    temp->POINT_NUMBER=sub[no_fra]->rigid_body[no_obj]->POINT_NUMBER;
    temp->TETRAHEDRON_NUMBER=sub[no_fra]->rigid_body[no_obj]->TETRAHEDRON_NUMBER;
    temp->TRIANGLE_NUMBER=sub[no_fra]->rigid_body[no_obj]->TRIANGLE_NUMBER;
    temp->INNER_POINT_NUMBER=sub[no_fra]->rigid_body[no_obj]->INNER_POINT_NUMBER;
    temp->OUTER_POINT_NUMBER=sub[no_fra]->rigid_body[no_obj]->OUTER_POINT_NUMBER;

    memory_allocate_for_rigid_body(temp);

    cusp::array2d<double,cusp::device_memory,cusp::column_major> T_t(3,3), Omega(3,3), Omega_T(3,3), Omega_Omega_T(3,3), dOmega(3,3), dOmega_T(3,3), Omega_Omega_T_plus_dOmega_T(3,3);
//    double T_t[3][3],Omega[3][3],Omega_T[3][3],Omega_Omega_T[3][3],dOmega[3][3],dOmega_T[3][3];
/*
    T_t[0][0]=sub[no_fra]->orientation[0][0];
    T_t[0][1]=sub[no_fra]->orientation[0][1];
    T_t[0][2]=sub[no_fra]->orientation[0][2];
    T_t[1][0]=sub[no_fra]->orientation[1][0];
    T_t[1][1]=sub[no_fra]->orientation[1][1];
    T_t[1][2]=sub[no_fra]->orientation[1][2];
    T_t[2][0]=sub[no_fra]->orientation[2][0];
    T_t[2][1]=sub[no_fra]->orientation[2][1];
    T_t[2][2]=sub[no_fra]->orientation[2][2];

    Omega[0][0]=0;
    Omega[0][1]=-sub[no_fra]->angular_velocity[2];
    Omega[0][2]=sub[no_fra]->angular_velocity[1];
    Omega[1][0]=sub[no_fra]->angular_velocity[2];
    Omega[1][1]=0;
    Omega[1][2]=-sub[no_fra]->angular_velocity[0];
    Omega[2][0]=-sub[no_fra]->angular_velocity[1]; 
    Omega[2][1]=sub[no_fra]->angular_velocity[0];   
    Omega[2][2]=0;

    Multiply(Omega,T_t,Omega_T);
    Multiply(Omega,Omega_T,Omega_Omega_T);

    dOmega[0][0]=0;
    dOmega[0][1]=-sub[no_fra]->angular_acceleration[2];  
    dOmega[0][2]=sub[no_fra]->angular_acceleration[1];
    dOmega[1][0]=sub[no_fra]->angular_acceleration[2];  
    dOmega[1][1]=0;                                      
    dOmega[1][2]=-sub[no_fra]->angular_acceleration[0];
    dOmega[2][0]=-sub[no_fra]->angular_acceleration[1]; 
    dOmega[2][1]=sub[no_fra]->angular_acceleration[0];
    dOmega[2][2]=0;

    Multiply(dOmega,T_t,dOmega_T);
*/
    T_t(0,0)=sub[no_fra]->orientation[0][0];
    T_t(0,1)=sub[no_fra]->orientation[0][1];
    T_t(0,2)=sub[no_fra]->orientation[0][2];
    T_t(1,0)=sub[no_fra]->orientation[1][0];
    T_t(1,1)=sub[no_fra]->orientation[1][1];
    T_t(1,2)=sub[no_fra]->orientation[1][2];
    T_t(2,0)=sub[no_fra]->orientation[2][0];
    T_t(2,1)=sub[no_fra]->orientation[2][1];
    T_t(2,2)=sub[no_fra]->orientation[2][2];

    Omega(0,0)=0;
    Omega(0,1)=-sub[no_fra]->angular_velocity[2];
    Omega(0,2)=sub[no_fra]->angular_velocity[1];
    Omega(1,0)=sub[no_fra]->angular_velocity[2];
    Omega(1,1)=0;
    Omega(1,2)=-sub[no_fra]->angular_velocity[0];
    Omega(2,0)=-sub[no_fra]->angular_velocity[1]; 
    Omega(2,1)=sub[no_fra]->angular_velocity[0];   
    Omega(2,2)=0;

    cusp::blas::gemm(Omega,T_t,Omega_T);
    cusp::blas::gemm(Omega,Omega_T,Omega_Omega_T);
    
    dOmega(0,0)=0;
    dOmega(0,1)=-sub[no_fra]->angular_acceleration[2];  
    dOmega(0,2)=sub[no_fra]->angular_acceleration[1];
    dOmega(1,0)=sub[no_fra]->angular_acceleration[2];  
    dOmega(1,1)=0;                                      
    dOmega(1,2)=-sub[no_fra]->angular_acceleration[0];
    dOmega(2,0)=-sub[no_fra]->angular_acceleration[1]; 
    dOmega(2,1)=sub[no_fra]->angular_acceleration[0];
    dOmega(2,2)=0;

    cusp::blas::gemm(dOmega,T_t,dOmega_T);
    
    
    cusp::array2d<double, cusp::device_memory, cusp::column_major> sub_position(3,temp->POINT_NUMBER), sub_velocity(3,temp->POINT_NUMBER), sub_acceleration(3,temp->POINT_NUMBER), TEMP(3,temp->POINT_NUMBER);
    
    for(int i = 0; i < 3; i++){
	thrust::fill(sub_position.row(i).begin(), sub_position.row(i).end(), sub[no_fra]->position[i]);
	thrust::fill(sub_velocity.row(i).begin(), sub_velocity.row(i).end(), sub[no_fra]->velocity[i]);
	thrust::fill(sub_acceleration.row(i).begin(), sub_acceleration.row(i).end(), sub[no_fra]->acceleration[i]);
    }

    
    
    cusp::blas::gemm(T_t,sub[no_fra]->rigid_body[no_obj]->XYZ,temp->XYZ);
    cusp::blas::geam(temp->XYZ,sub_position,temp->XYZ,1.0,1.0);
    
    cusp::blas::gemm(T_t,sub[no_fra]->rigid_body[no_obj]->UVW,temp->UVW);
    cusp::blas::gemm(Omega_T,sub[no_fra]->rigid_body[no_obj]->XYZ,TEMP);
    cusp::blas::geam(temp->UVW,TEMP,temp->UVW,1.0,1.0);
    cusp::blas::geam(temp->UVW,sub_velocity,temp->UVW,1.0,1.0);
    
    cusp::blas::geam(Omega_Omega_T,dOmega_T,Omega_Omega_T_plus_dOmega_T,1.0,1.0);
    cusp::blas::gemm(T_t,sub[no_fra]->rigid_body[no_obj]->ACC,temp->ACC);
    cusp::blas::gemm(Omega_T,sub[no_fra]->rigid_body[no_obj]->UVW,TEMP);
    cusp::blas::geam(temp->ACC,TEMP,temp->ACC,1.0,2.0);
    cusp::blas::gemm(Omega_Omega_T_plus_dOmega_T,sub[no_fra]->rigid_body[no_obj]->XYZ,TEMP);
    cusp::blas::geam(temp->ACC,TEMP,temp->ACC,1.0,1.0);
    cusp::blas::geam(temp->ACC,sub_acceleration,temp->ACC,1.0,1.0);
    
    cusp::blas::gemm(T_t,sub[no_fra]->rigid_body[no_obj]->OUTER_NORMAL_VECTOR,temp->OUTER_NORMAL_VECTOR);
    
    
    temp->AREA = sub[no_fra]->rigid_body[no_obj]->AREA;
    temp->INNERMARK = sub[no_fra]->rigid_body[no_obj]->INNERMARK;
    temp->OUTERMARK = sub[no_fra]->rigid_body[no_obj]->OUTERMARK;
    
    temp->TETRAHEDRON = sub[no_fra]->rigid_body[no_obj]->TETRAHEDRON;
    temp->TRIANGLE = sub[no_fra]->rigid_body[no_obj]->TRIANGLE;
    
/*    
    for (int s=0; s<temp->POINT_NUMBER; s++)
	{
		temp->XYZ[0][s] = sub[no_fra]->position[0]
                        + T_t[0][0]*sub[no_fra]->rigid_body[no_obj]->XYZ[0][s]
                        + T_t[0][1]*sub[no_fra]->rigid_body[no_obj]->XYZ[1][s]
                        + T_t[0][2]*sub[no_fra]->rigid_body[no_obj]->XYZ[2][s];
		temp->XYZ[1][s] = sub[no_fra]->position[1]
                        + T_t[1][0]*sub[no_fra]->rigid_body[no_obj]->XYZ[0][s]
                        + T_t[1][1]*sub[no_fra]->rigid_body[no_obj]->XYZ[1][s]
                        + T_t[1][2]*sub[no_fra]->rigid_body[no_obj]->XYZ[2][s];
		temp->XYZ[2][s] = sub[no_fra]->position[2]
                        + T_t[2][0]*sub[no_fra]->rigid_body[no_obj]->XYZ[0][s]
                        + T_t[2][1]*sub[no_fra]->rigid_body[no_obj]->XYZ[1][s]
                        + T_t[2][2]*sub[no_fra]->rigid_body[no_obj]->XYZ[2][s];
		
		temp->UVW[0][s] = sub[no_fra]->velocity[0]
                        + T_t[0][0]*sub[no_fra]->rigid_body[no_obj]->UVW[0][s]
                        + T_t[0][1]*sub[no_fra]->rigid_body[no_obj]->UVW[1][s]
                        + T_t[0][2]*sub[no_fra]->rigid_body[no_obj]->UVW[2][s]
                        + Omega_T[0][0]*sub[no_fra]->rigid_body[no_obj]->XYZ[0][s]
                        + Omega_T[0][1]*sub[no_fra]->rigid_body[no_obj]->XYZ[1][s]
                        + Omega_T[0][2]*sub[no_fra]->rigid_body[no_obj]->XYZ[2][s];
		temp->UVW[1][s] = sub[no_fra]->velocity[1]
                        + T_t[1][0]*sub[no_fra]->rigid_body[no_obj]->UVW[0][s]
                        + T_t[1][1]*sub[no_fra]->rigid_body[no_obj]->UVW[1][s]
                        + T_t[1][2]*sub[no_fra]->rigid_body[no_obj]->UVW[2][s]
                        + Omega_T[1][0]*sub[no_fra]->rigid_body[no_obj]->XYZ[0][s]
                        + Omega_T[1][1]*sub[no_fra]->rigid_body[no_obj]->XYZ[1][s]
                        + Omega_T[1][2]*sub[no_fra]->rigid_body[no_obj]->XYZ[2][s];
		temp->UVW[2][s] = sub[no_fra]->velocity[2]
                        + T_t[2][0]*sub[no_fra]->rigid_body[no_obj]->UVW[0][s]
                        + T_t[2][1]*sub[no_fra]->rigid_body[no_obj]->UVW[1][s]
                        + T_t[2][2]*sub[no_fra]->rigid_body[no_obj]->UVW[2][s]
                        + Omega_T[2][0]*sub[no_fra]->rigid_body[no_obj]->XYZ[0][s]
                        + Omega_T[2][1]*sub[no_fra]->rigid_body[no_obj]->XYZ[1][s]
                        + Omega_T[2][2]*sub[no_fra]->rigid_body[no_obj]->XYZ[2][s];
		
		temp->ACC[0][s] = sub[no_fra]->acceleration[0]
                        + T_t[0][0]*sub[no_fra]->rigid_body[no_obj]->ACC[0][s]
                        + T_t[0][1]*sub[no_fra]->rigid_body[no_obj]->ACC[1][s]
                        + T_t[0][2]*sub[no_fra]->rigid_body[no_obj]->ACC[2][s]
                        + 2*Omega_T[0][0]*sub[no_fra]->rigid_body[no_obj]->UVW[0][s]
                        + 2*Omega_T[0][1]*sub[no_fra]->rigid_body[no_obj]->UVW[1][s]
                        + 2*Omega_T[0][2]*sub[no_fra]->rigid_body[no_obj]->UVW[2][s]
                        + (Omega_Omega_T[0][0]+dOmega_T[0][0])*sub[no_fra]->rigid_body[no_obj]->XYZ[0][s]
                        + (Omega_Omega_T[0][1]+dOmega_T[0][1])*sub[no_fra]->rigid_body[no_obj]->XYZ[1][s]
                        + (Omega_Omega_T[0][2]+dOmega_T[0][2])*sub[no_fra]->rigid_body[no_obj]->XYZ[2][s];	    
		temp->ACC[1][s] = sub[no_fra]->acceleration[1]
                        + T_t[1][0]*sub[no_fra]->rigid_body[no_obj]->ACC[0][s]
                        + T_t[1][1]*sub[no_fra]->rigid_body[no_obj]->ACC[1][s]
                        + T_t[1][2]*sub[no_fra]->rigid_body[no_obj]->ACC[2][s]
                        + 2*Omega_T[1][0]*sub[no_fra]->rigid_body[no_obj]->UVW[0][s]
                        + 2*Omega_T[1][1]*sub[no_fra]->rigid_body[no_obj]->UVW[1][s]
                        + 2*Omega_T[1][2]*sub[no_fra]->rigid_body[no_obj]->UVW[2][s]
                        + (Omega_Omega_T[1][0]+dOmega_T[1][0])*sub[no_fra]->rigid_body[no_obj]->XYZ[0][s]
                        + (Omega_Omega_T[1][1]+dOmega_T[1][1])*sub[no_fra]->rigid_body[no_obj]->XYZ[1][s]
                        + (Omega_Omega_T[1][2]+dOmega_T[1][2])*sub[no_fra]->rigid_body[no_obj]->XYZ[2][s];
		temp->ACC[2][s] = sub[no_fra]->acceleration[2]
                        + T_t[2][0]*sub[no_fra]->rigid_body[no_obj]->ACC[0][s]
                        + T_t[2][1]*sub[no_fra]->rigid_body[no_obj]->ACC[1][s]
                        + T_t[2][2]*sub[no_fra]->rigid_body[no_obj]->ACC[2][s]
                        + 2*Omega_T[2][0]*sub[no_fra]->rigid_body[no_obj]->UVW[0][s]
                        + 2*Omega_T[2][1]*sub[no_fra]->rigid_body[no_obj]->UVW[1][s]
                        + 2*Omega_T[2][2]*sub[no_fra]->rigid_body[no_obj]->UVW[2][s]
                        + (Omega_Omega_T[2][0]+dOmega_T[2][0])*sub[no_fra]->rigid_body[no_obj]->XYZ[0][s]
                        + (Omega_Omega_T[2][1]+dOmega_T[2][1])*sub[no_fra]->rigid_body[no_obj]->XYZ[1][s]
                        + (Omega_Omega_T[2][2]+dOmega_T[2][2])*sub[no_fra]->rigid_body[no_obj]->XYZ[2][s];
		
		temp->OUTER_NORMAL_VECTOR[0][s] = 
                          T_t[0][0]*sub[no_fra]->rigid_body[no_obj]->OUTER_NORMAL_VECTOR[0][s]
                        + T_t[0][1]*sub[no_fra]->rigid_body[no_obj]->OUTER_NORMAL_VECTOR[1][s]
                        + T_t[0][2]*sub[no_fra]->rigid_body[no_obj]->OUTER_NORMAL_VECTOR[2][s];
		temp->OUTER_NORMAL_VECTOR[1][s] = 
                          T_t[1][0]*sub[no_fra]->rigid_body[no_obj]->OUTER_NORMAL_VECTOR[0][s]
                        + T_t[1][1]*sub[no_fra]->rigid_body[no_obj]->OUTER_NORMAL_VECTOR[1][s]
                        + T_t[1][2]*sub[no_fra]->rigid_body[no_obj]->OUTER_NORMAL_VECTOR[2][s];
		temp->OUTER_NORMAL_VECTOR[2][s] = 
                          T_t[2][0]*sub[no_fra]->rigid_body[no_obj]->OUTER_NORMAL_VECTOR[0][s]
                        + T_t[2][1]*sub[no_fra]->rigid_body[no_obj]->OUTER_NORMAL_VECTOR[1][s]
                        + T_t[2][2]*sub[no_fra]->rigid_body[no_obj]->OUTER_NORMAL_VECTOR[2][s];
		
		temp->AREA[s]=sub[no_fra]->rigid_body[no_obj]->AREA[s];
		
		temp->INNERMARK[s]=sub[no_fra]->rigid_body[no_obj]->INNERMARK[s];
		temp->OUTERMARK[s]=sub[no_fra]->rigid_body[no_obj]->OUTERMARK[s];
	}
	
	
	for (int j=0; j<temp->TETRAHEDRON_NUMBER; j++)
	{
		temp->TETRAHEDRON[0][j]=sub[no_fra]->rigid_body[no_obj]->TETRAHEDRON[0][j];
		temp->TETRAHEDRON[1][j]=sub[no_fra]->rigid_body[no_obj]->TETRAHEDRON[1][j];
		temp->TETRAHEDRON[2][j]=sub[no_fra]->rigid_body[no_obj]->TETRAHEDRON[2][j];
		temp->TETRAHEDRON[3][j]=sub[no_fra]->rigid_body[no_obj]->TETRAHEDRON[3][j];
	}
	
	for (int j=0; j<temp->TRIANGLE_NUMBER; j++)
	{
		temp->TRIANGLE[0][j]=sub[no_fra]->rigid_body[no_obj]->TRIANGLE[0][j];
		temp->TRIANGLE[1][j]=sub[no_fra]->rigid_body[no_obj]->TRIANGLE[1][j];
		temp->TRIANGLE[2][j]=sub[no_fra]->rigid_body[no_obj]->TRIANGLE[2][j];
	}
*/	
	
	return temp;
}
/*
void Ref_Frame :: change_ref_frame()
{
    if(subframe_number==0){cout<<"No sub frame exist!"<<endl;}
    if(obj_number==0){cout<<"No object exist in current frame!"<<endl;}
    double T_t[3][3],Omega[3][3],Omega_T[3][3],Omega_Omega_T[3][3],dOmega[3][3],dOmega_T[3][3];
	
    for(int k=0;k<obj_number;k++)
    {
	for(int i=0;i<subframe_number;i++)
	{
	    if(sub[i]->obj_number==0){cout<<"No object exist in sub frame!"<<endl;continue;}
	
	    T_t[0][0]=sub[i]->orientation[0][0];
	    T_t[0][1]=sub[i]->orientation[0][1];
	    T_t[0][2]=sub[i]->orientation[0][2];
	    T_t[1][0]=sub[i]->orientation[1][0];
	    T_t[1][1]=sub[i]->orientation[1][1];
	    T_t[1][2]=sub[i]->orientation[1][2];
	    T_t[2][0]=sub[i]->orientation[2][0];
	    T_t[2][1]=sub[i]->orientation[2][1];
	    T_t[2][2]=sub[i]->orientation[2][2];
			
	    Omega[0][0]=0;
	    Omega[0][1]=-sub[i]->angular_velocity[2];  
	    Omega[0][2]=sub[i]->angular_velocity[1];
	    Omega[1][0]=sub[i]->angular_velocity[2];  
	    Omega[1][1]=0;
	    Omega[1][2]=-sub[i]->angular_velocity[0];
	    Omega[2][0]=-sub[i]->angular_velocity[1]; 
	    Omega[2][1]=sub[i]->angular_velocity[0];   
	    Omega[2][2]=0;
			
	    Multiply(Omega,T_t,Omega_T);
			
	    Multiply(Omega,Omega_T,Omega_Omega_T);
			
	    dOmega[0][0]=0;
	    dOmega[0][1]=-sub[i]->angular_acceleration[2];  
	    dOmega[0][2]=sub[i]->angular_acceleration[1];
	    dOmega[1][0]=sub[i]->angular_acceleration[2];  
	    dOmega[1][1]=0;
	    dOmega[1][2]=-sub[i]->angular_acceleration[0];
	    dOmega[2][0]=-sub[i]->angular_acceleration[1]; 
	    dOmega[2][1]=sub[i]->angular_acceleration[0];
	    dOmega[2][2]=0;
			
	    Multiply(dOmega,T_t,dOmega_T);
			
	    for(int j=0;j<sub[i]->obj_number;j++)
	    {
		if((rigid_body[k]->IDENTITY==sub[i]->rigid_body[j]->IDENTITY)&&(rigid_body[k]->LEVEL==sub[i]->rigid_body[j]->LEVEL-1))
		{
		    for(int s=0;s<rigid_body[k]->POINT_NUMBER;s++)
		    {
			rigid_body[k]->XYZ[0][s] = sub[i]->position[0]
                        + T_t[0][0]*sub[i]->rigid_body[j]->XYZ[0][s]
                        + T_t[0][1]*sub[i]->rigid_body[j]->XYZ[1][s]
                        + T_t[0][2]*sub[i]->rigid_body[j]->XYZ[2][s];
			rigid_body[k]->XYZ[1][s] = sub[i]->position[1]
                        + T_t[1][0]*sub[i]->rigid_body[j]->XYZ[0][s]
                        + T_t[1][1]*sub[i]->rigid_body[j]->XYZ[1][s]
                        + T_t[1][2]*sub[i]->rigid_body[j]->XYZ[2][s];
			rigid_body[k]->XYZ[2][s] = sub[i]->position[2]
                        + T_t[2][0]*sub[i]->rigid_body[j]->XYZ[0][s]
                        + T_t[2][1]*sub[i]->rigid_body[j]->XYZ[1][s]
                        + T_t[2][2]*sub[i]->rigid_body[j]->XYZ[2][s];
						
			rigid_body[k]->UVW[0][s] = sub[i]->velocity[0]
                        + T_t[0][0]*sub[i]->rigid_body[j]->UVW[0][s]
                        + T_t[0][1]*sub[i]->rigid_body[j]->UVW[1][s]
                        + T_t[0][2]*sub[i]->rigid_body[j]->UVW[2][s]
                        + Omega_T[0][0]*sub[i]->rigid_body[j]->XYZ[0][s]
                        + Omega_T[0][1]*sub[i]->rigid_body[j]->XYZ[1][s]
                        + Omega_T[0][2]*sub[i]->rigid_body[j]->XYZ[2][s];
			rigid_body[k]->UVW[1][s] = sub[i]->velocity[1]
                        + T_t[1][0]*sub[i]->rigid_body[j]->UVW[0][s]
                        + T_t[1][1]*sub[i]->rigid_body[j]->UVW[1][s]
                        + T_t[1][2]*sub[i]->rigid_body[j]->UVW[2][s]
                        + Omega_T[1][0]*sub[i]->rigid_body[j]->XYZ[0][s]
                        + Omega_T[1][1]*sub[i]->rigid_body[j]->XYZ[1][s]
                        + Omega_T[1][2]*sub[i]->rigid_body[j]->XYZ[2][s];
			rigid_body[k]->UVW[2][s] = sub[i]->velocity[2]
                        + T_t[2][0]*sub[i]->rigid_body[j]->UVW[0][s]
                        + T_t[2][1]*sub[i]->rigid_body[j]->UVW[1][s]
                        + T_t[2][2]*sub[i]->rigid_body[j]->UVW[2][s]
                        + Omega_T[2][0]*sub[i]->rigid_body[j]->XYZ[0][s]
                        + Omega_T[2][1]*sub[i]->rigid_body[j]->XYZ[1][s]
                        + Omega_T[2][2]*sub[i]->rigid_body[j]->XYZ[2][s];
						
			rigid_body[k]->ACC[0][s] = sub[i]->acceleration[0]
                        + T_t[0][0]*sub[i]->rigid_body[j]->ACC[0][s]
                        + T_t[0][1]*sub[i]->rigid_body[j]->ACC[1][s]
                        + T_t[0][2]*sub[i]->rigid_body[j]->ACC[2][s]
                        + 2*Omega_T[0][0]*sub[i]->rigid_body[j]->UVW[0][s]
                        + 2*Omega_T[0][1]*sub[i]->rigid_body[j]->UVW[1][s]
                        + 2*Omega_T[0][2]*sub[i]->rigid_body[j]->UVW[2][s]
                        + (Omega_Omega_T[0][0]+dOmega_T[0][0])*sub[i]->rigid_body[j]->XYZ[0][s]
                        + (Omega_Omega_T[0][1]+dOmega_T[0][1])*sub[i]->rigid_body[j]->XYZ[1][s]
                        + (Omega_Omega_T[0][2]+dOmega_T[0][2])*sub[i]->rigid_body[j]->XYZ[2][s];
			rigid_body[k]->ACC[1][s] = sub[i]->acceleration[1]
                        +T_t[1][0]*sub[i]->rigid_body[j]->ACC[0][s]
                        +T_t[1][1]*sub[i]->rigid_body[j]->ACC[1][s]
                        +T_t[1][2]*sub[i]->rigid_body[j]->ACC[2][s]
                        +2*Omega_T[1][0]*sub[i]->rigid_body[j]->UVW[0][s]
                        +2*Omega_T[1][1]*sub[i]->rigid_body[j]->UVW[1][s]
                        +2*Omega_T[1][2]*sub[i]->rigid_body[j]->UVW[2][s]
                        +(Omega_Omega_T[1][0]+dOmega_T[1][0])*sub[i]->rigid_body[j]->XYZ[0][s]
                        +(Omega_Omega_T[1][1]+dOmega_T[1][1])*sub[i]->rigid_body[j]->XYZ[1][s]
                        +(Omega_Omega_T[1][2]+dOmega_T[1][2])*sub[i]->rigid_body[j]->XYZ[2][s];
			rigid_body[k]->ACC[2][s]=sub[i]->acceleration[2]
                        +T_t[2][0]*sub[i]->rigid_body[j]->ACC[0][s]
                        +T_t[2][1]*sub[i]->rigid_body[j]->ACC[1][s]
                        +T_t[2][2]*sub[i]->rigid_body[j]->ACC[2][s]
                        +2*Omega_T[2][0]*sub[i]->rigid_body[j]->UVW[0][s]
                        +2*Omega_T[2][1]*sub[i]->rigid_body[j]->UVW[1][s]
                        +2*Omega_T[2][2]*sub[i]->rigid_body[j]->UVW[2][s]
                        +(Omega_Omega_T[2][0]+dOmega_T[2][0])*sub[i]->rigid_body[j]->XYZ[0][s]
                        +(Omega_Omega_T[2][1]+dOmega_T[2][1])*sub[i]->rigid_body[j]->XYZ[1][s]
                        +(Omega_Omega_T[2][2]+dOmega_T[2][2])*sub[i]->rigid_body[j]->XYZ[2][s];
													  
			rigid_body[k]->OUTER_NORMAL_VECTOR[0][s] = 
                          T_t[0][0]*sub[i]->rigid_body[j]->OUTER_NORMAL_VECTOR[0][s]
                        + T_t[0][1]*sub[i]->rigid_body[j]->OUTER_NORMAL_VECTOR[1][s]
                        + T_t[0][2]*sub[i]->rigid_body[j]->OUTER_NORMAL_VECTOR[2][s];
			rigid_body[k]->OUTER_NORMAL_VECTOR[1][s] = 
                          T_t[1][0]*sub[i]->rigid_body[j]->OUTER_NORMAL_VECTOR[0][s]
                        + T_t[1][1]*sub[i]->rigid_body[j]->OUTER_NORMAL_VECTOR[1][s]
                        + T_t[1][2]*sub[i]->rigid_body[j]->OUTER_NORMAL_VECTOR[2][s];
			rigid_body[k]->OUTER_NORMAL_VECTOR[2][s] = 
			              T_t[2][0]*sub[i]->rigid_body[j]->OUTER_NORMAL_VECTOR[0][s]
			            + T_t[2][1]*sub[i]->rigid_body[j]->OUTER_NORMAL_VECTOR[1][s]
			            + T_t[2][2]*sub[i]->rigid_body[j]->OUTER_NORMAL_VECTOR[2][s];

		    }
		}
	    }
	}
    }
}
*/

void REF_FRAME :: change_ref_frame()
{
    if(subframe_number==0){std::cout<<"No sub frame exist!"<<std::endl;}
    if(obj_number==0){std::cout<<"No object exist in current frame!"<<std::endl;}

    cusp::array2d<double,cusp::device_memory,cusp::column_major> T_t(3,3), Omega(3,3), Omega_T(3,3), Omega_Omega_T(3,3), dOmega(3,3), dOmega_T(3,3), Omega_Omega_T_plus_dOmega_T(3,3);
    
    for(int k=0;k<obj_number;k++)
    {
	for(int i=0;i<subframe_number;i++)
	{
	    if(sub[i]->obj_number==0){std::cout<<"No object exist in sub frame!"<<std::endl;continue;}

	    T_t(0,0)=sub[i]->orientation[0][0];
	    T_t(0,1)=sub[i]->orientation[0][1];
	    T_t(0,2)=sub[i]->orientation[0][2];
	    T_t(1,0)=sub[i]->orientation[1][0];
	    T_t(1,1)=sub[i]->orientation[1][1];
	    T_t(1,2)=sub[i]->orientation[1][2];
	    T_t(2,0)=sub[i]->orientation[2][0];
	    T_t(2,1)=sub[i]->orientation[2][1];
	    T_t(2,2)=sub[i]->orientation[2][2];

	    Omega(0,0)=0;
	    Omega(0,1)=-sub[i]->angular_velocity[2];
	    Omega(0,2)=sub[i]->angular_velocity[1];
	    Omega(1,0)=sub[i]->angular_velocity[2];
	    Omega(1,1)=0;
	    Omega(1,2)=-sub[i]->angular_velocity[0];
	    Omega(2,0)=-sub[i]->angular_velocity[1]; 
	    Omega(2,1)=sub[i]->angular_velocity[0];   
	    Omega(2,2)=0;

	    cusp::blas::gemm(Omega,T_t,Omega_T);
	    cusp::blas::gemm(Omega,Omega_T,Omega_Omega_T);

	    dOmega(0,0)=0;
	    dOmega(0,1)=-sub[i]->angular_acceleration[2];  
	    dOmega(0,2)=sub[i]->angular_acceleration[1];
	    dOmega(1,0)=sub[i]->angular_acceleration[2];  
	    dOmega(1,1)=0;                                      
	    dOmega(1,2)=-sub[i]->angular_acceleration[0];
	    dOmega(2,0)=-sub[i]->angular_acceleration[1]; 
	    dOmega(2,1)=sub[i]->angular_acceleration[0];
	    dOmega(2,2)=0;
	    
	    cusp::blas::gemm(dOmega,T_t,dOmega_T);
	    cusp::array2d<double, cusp::device_memory, cusp::column_major> sub_position(3,rigid_body[k]->POINT_NUMBER), sub_velocity(3,rigid_body[k]->POINT_NUMBER), sub_acceleration(3,rigid_body[k]->POINT_NUMBER), TEMP(3,rigid_body[k]->POINT_NUMBER);

	    for(int j=0;j<sub[i]->obj_number;j++)
	    {
		if((rigid_body[k]->IDENTITY==sub[i]->rigid_body[j]->IDENTITY)&&(rigid_body[k]->LEVEL==sub[i]->rigid_body[j]->LEVEL-1))
		{
		    for(int p = 0; p < 3; p++){
			thrust::fill(sub_position.row(p).begin(), sub_position.row(p).end(), sub[i]->position[p]);
			thrust::fill(sub_velocity.row(p).begin(), sub_velocity.row(p).end(), sub[i]->velocity[p]);
			thrust::fill(sub_acceleration.row(p).begin(), sub_acceleration.row(p).end(), sub[i]->acceleration[p]);
		    }
		    
		    cusp::blas::gemm(T_t,sub[i]->rigid_body[j]->XYZ,rigid_body[k]->XYZ);
		    cusp::blas::geam(rigid_body[k]->XYZ,sub_position,rigid_body[k]->XYZ,1.0,1.0);
		    
		    cusp::blas::gemm(T_t,sub[i]->rigid_body[j]->UVW,rigid_body[k]->UVW);
		    cusp::blas::gemm(Omega_T,sub[i]->rigid_body[j]->XYZ,TEMP);
		    cusp::blas::geam(rigid_body[k]->UVW,TEMP,rigid_body[k]->UVW,1.0,1.0);
		    cusp::blas::geam(rigid_body[k]->UVW,sub_velocity,rigid_body[k]->UVW,1.0,1.0);
		    
		    cusp::blas::geam(Omega_Omega_T,dOmega_T,Omega_Omega_T_plus_dOmega_T,1.0,1.0);
		    cusp::blas::gemm(T_t,sub[i]->rigid_body[j]->ACC,rigid_body[k]->ACC);
		    cusp::blas::gemm(Omega_T,sub[i]->rigid_body[j]->UVW,TEMP);
		    cusp::blas::geam(rigid_body[k]->ACC,TEMP,rigid_body[k]->ACC,1.0,2.0);
		    cusp::blas::gemm(Omega_Omega_T_plus_dOmega_T,sub[i]->rigid_body[j]->XYZ,TEMP);
		    cusp::blas::geam(rigid_body[k]->ACC,TEMP,rigid_body[k]->ACC,1.0,1.0);
		    cusp::blas::geam(rigid_body[k]->ACC,sub_acceleration,rigid_body[k]->ACC,1.0,1.0);
		    
		    cusp::blas::gemm(T_t,sub[i]->rigid_body[j]->OUTER_NORMAL_VECTOR,rigid_body[k]->OUTER_NORMAL_VECTOR);
		}
	    }
	}
    }
}

