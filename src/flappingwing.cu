#include <iostream>
#include <fstream>
#include <omp.h>

#include <cusp/print.h>

#include "../include/flow_field.h"
#include "../include/timer.h"
#include "../include/mesh.h"
#include "../include/linklist.h"

#include "../include/frame_structure.h"
#include "../include/control_module.h"
#include "../include/flap_pattern.h"
#include "../include/insect_parameters.h"
#include "../include/solver.h"
#include "../include/write_file.h"
#include "../include/global.h"

bool GPU_detection();
void CPU_detection();


void Initialize_Ref_Frame(SOLVER* solver, REF_FRAME* lab, REF_FRAME* insect, REF_FRAME* wingplane_1, REF_FRAME* wingplane_2, REF_FRAME* wing_1, REF_FRAME* wing_2)
{
	lab->level=0;
	insect->level=1;
	wingplane_1->level=2;
	wingplane_2->level=2;
	wing_1->level=3;
	wing_2->level=3;
	
	lab->super=NULL;
	insect->super=lab;
	wingplane_1->super=insect;
	wingplane_2->super=insect;
	wing_1->super=wingplane_1;
	wing_2->super=wingplane_2;
	
	lab->subframe_number=1;
	insect->subframe_number=2;
	wingplane_1->subframe_number=1;
	wingplane_2->subframe_number=1;
	wing_1->subframe_number=0;
	wing_2->subframe_number=0;
	
	lab->sub            = new REF_FRAME* [lab->subframe_number];
	insect->sub         = new REF_FRAME* [insect->subframe_number];
	wingplane_1->sub    = new REF_FRAME* [wingplane_1->subframe_number];
	wingplane_2->sub    = new REF_FRAME* [wingplane_2->subframe_number];
	wing_1->sub         = NULL;
	wing_2->sub         = NULL;
	lab->sub[0]         = insect;
	insect->sub[0]      = wingplane_1;
	insect->sub[1]      = wingplane_2;
	wingplane_1->sub[0] = wing_1;
	wingplane_2->sub[0] = wing_2;
	
	lab->position[0]=0;lab->position[1]=0;lab->position[2]=0;
	lab->velocity[0]=0;lab->velocity[1]=0;lab->velocity[2]=0;
	lab->angle[0]=0;lab->angle[1]=0;lab->angle[2]=0;
	lab->angular_velocity[0]=0;lab->angular_velocity[1]=0;lab->angular_velocity[2]=0;

	
	insect_parameters.Insect_Parameters_Initialize();
	insect_parameters.Insect_Parameters_Nondimensionalize();
	
	for(int i=0;i<3;i++) {
//		lastbodycentre[i]=8.0;
		CoM[i]=0;  // in flying frame
		Wing_Root[i] = -insect_parameters.Nondimensional_CoM[i];
	}
	for(int i=0;i<3;i++) {insect->momentum[i]=insect_parameters.Nondimensional_Mass*insect->velocity[i];}
	//input  local Inertia
	for(int i=0;i<3;i++) {
	    for(int j=0;j<3;j++) {
		insect->local_inertia[i][j]=insect_parameters.Nondimensional_MOI[i][j];
	    }
	}
	insect->local_inertia[0][0] -= insect_parameters.Nondimensional_Mass
		* (  insect_parameters.Nondimensional_CoM[1]*insect_parameters.Nondimensional_CoM[1]
		   + insect_parameters.Nondimensional_CoM[2]*insect_parameters.Nondimensional_CoM[2] );
	insect->local_inertia[0][1] -= insect_parameters.Nondimensional_Mass 
	    * ( -insect_parameters.Nondimensional_CoM[0]*insect_parameters.Nondimensional_CoM[1] );
	insect->local_inertia[0][2] -= insect_parameters.Nondimensional_Mass
	    * ( -insect_parameters.Nondimensional_CoM[0]*insect_parameters.Nondimensional_CoM[2] );
	insect->local_inertia[1][0] -= insect_parameters.Nondimensional_Mass
	    * ( -insect_parameters.Nondimensional_CoM[1]*insect_parameters.Nondimensional_CoM[0] );
	insect->local_inertia[1][1] -= insect_parameters.Nondimensional_Mass
	    * (  insect_parameters.Nondimensional_CoM[0]*insect_parameters.Nondimensional_CoM[0]
	       + insect_parameters.Nondimensional_CoM[2]*insect_parameters.Nondimensional_CoM[2] );
	insect->local_inertia[1][2] -= insect_parameters.Nondimensional_Mass
	    * ( -insect_parameters.Nondimensional_CoM[1]*insect_parameters.Nondimensional_CoM[2] );
	insect->local_inertia[2][0] -= insect_parameters.Nondimensional_Mass
	    * ( -insect_parameters.Nondimensional_CoM[2]*insect_parameters.Nondimensional_CoM[0] );
	insect->local_inertia[2][1] -= insect_parameters.Nondimensional_Mass
	    * ( -insect_parameters.Nondimensional_CoM[2]*insect_parameters.Nondimensional_CoM[1] );
	insect->local_inertia[2][2] -= insect_parameters.Nondimensional_Mass
	    * (  insect_parameters.Nondimensional_CoM[1]*insect_parameters.Nondimensional_CoM[1]
	       + insect_parameters.Nondimensional_CoM[0]*insect_parameters.Nondimensional_CoM[0] );
	
	double O_tran[3][3],temp[3][3];
	Transpose(insect->orientation,O_tran);
	Multiply(insect->orientation,insect->local_inertia,temp);
	Multiply(temp,O_tran,insect->inertia);
	for(int i=0;i<3;i++)
	{
		insect->angular_momentum[i] = insect->inertia[i][0]*insect->angular_velocity[0]
		                             + insect->inertia[i][1]*insect->angular_velocity[1]
		                             + insect->inertia[i][2]*insect->angular_velocity[2];
	}

	
//	Assign wing motion functions
	my_flap_pattern.Set_Flap_Pattern( Wing_Root, bodycentre );
	my_flap_pattern.solver = solver;
	
	insect->compute_x          = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: insect_X);
	insect->compute_y          = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: insect_Y);
	insect->compute_z          = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: insect_Z);
	insect->compute_phi        = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: insect_Phi);
	insect->compute_theta      = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: insect_Theta);
	insect->compute_psi        = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: insect_Psi);
	
	wingplane_1->compute_x     = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: wingplane1_X);
	wingplane_1->compute_y     = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: wingplane1_Y);
	wingplane_1->compute_z     = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: wingplane1_Z);
	wingplane_1->compute_phi   = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: wingplane1_Phi);
	wingplane_1->compute_theta = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: wingplane1_Theta);
	wingplane_1->compute_psi   = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: wingplane1_Psi);
	
	wingplane_2->compute_x     = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: wingplane2_X);
	wingplane_2->compute_y     = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: wingplane2_Y);
	wingplane_2->compute_z     = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: wingplane2_Z);
	wingplane_2->compute_phi   = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: wingplane2_Phi);
	wingplane_2->compute_theta = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: wingplane2_Theta);
	wingplane_2->compute_psi   = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: wingplane2_Psi);
	
	wing_1->compute_x     	   = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: wing1_X);
	wing_1->compute_y     	   = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: wing1_Y);
	wing_1->compute_z     	   = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: wing1_Z);
	wing_1->compute_phi        = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: wing1_Phi);
	wing_1->compute_theta      = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: wing1_Theta);
	wing_1->compute_psi        = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: wing1_Psi);
	
	wing_2->compute_x          = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: wing2_X);
	wing_2->compute_y          = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: wing2_Y);
	wing_2->compute_z          = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: wing2_Z);
	wing_2->compute_phi        = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: wing2_Phi);
	wing_2->compute_theta      = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: wing2_Theta);
	wing_2->compute_psi        = AttitudeEvaluator( &my_flap_pattern, &FLAPPATTERN :: wing2_Psi);
	
	lab->Time=0;
	insect->Time=0;
	wingplane_1->Time=0;
	wingplane_2->Time=0;
	wing_1->Time=0;
	wing_2->Time=0;
	
	insect->update_all();
	wingplane_1->update_all();
	wingplane_2->update_all();
	wing_1->update_all();
	wing_2->update_all();

	
	lab->obj_number=3;
	insect->obj_number=3;
	wingplane_1->obj_number=1;
	wingplane_2->obj_number=1;
	wing_1->obj_number=1;
	wing_2->obj_number=1;
	
	lab->rigid_body           = new RIGID_OBJ* [lab->obj_number];
	insect->rigid_body        = new RIGID_OBJ* [insect->obj_number];
	wingplane_1->rigid_body   = new RIGID_OBJ* [wingplane_1->obj_number];
	wingplane_2->rigid_body   = new RIGID_OBJ* [wingplane_2->obj_number];
	wing_1->rigid_body        = new RIGID_OBJ* [wing_1->obj_number];
	wing_2->rigid_body        = new RIGID_OBJ* [wing_2->obj_number];
	
	wing_1->rigid_body[0]=wing_1->get_rigid_body_from_file("../mesh/realwing1-wingthick0.02-cloudthick0.08-s0.015-g1.2-m0.025.dat");
	wing_2->rigid_body[0]=wing_2->get_rigid_body_from_file("../mesh/realwing2-wingthick0.02-cloudthick0.08-s0.015-g1.2-m0.025.dat");
		
	wingplane_1->rigid_body[0]=wingplane_1->get_rigid_body_from_subframe(0,0);
	wingplane_2->rigid_body[0]=wingplane_2->get_rigid_body_from_subframe(0,0);
	
	insect->rigid_body[0]=insect->get_rigid_body_from_file("../mesh/y-rotate-30degrees-fruitfly-cloudthick0.04-meshsize0.02.dat");
	

/*   need to invoke the thrust functions   */
    for(int i = 0; i < 3; i++)
    {
	thrust::transform(insect->rigid_body[0]->XYZ.row(i).begin(), insect->rigid_body[0]->XYZ.row(i).end(),
			  thrust::make_constant_iterator(-insect_parameters.Nondimensional_CoM[i]),
			  insect->rigid_body[0]->XYZ.row(i).begin(), thrust::plus<double>() );
    }
/*
	for (int i=0; i<insect->rigid_body[0]->POINT_NUMBER; i++) {
		insect->rigid_body[0]->XYZ(0,i) -= insect_parameters.Nondimensional_CoM[0];
		insect->rigid_body[0]->XYZ(1,i) -= insect_parameters.Nondimensional_CoM[1];
		insect->rigid_body[0]->XYZ(2,i) -= insect_parameters.Nondimensional_CoM[2];
	}
*/	
	insect->rigid_body[1]=insect->get_rigid_body_from_subframe(0,0);
	insect->rigid_body[2]=insect->get_rigid_body_from_subframe(1,0);
	
	lab->rigid_body[0]=lab->get_rigid_body_from_subframe(0,0);
	lab->rigid_body[1]=lab->get_rigid_body_from_subframe(0,1);
	lab->rigid_body[2]=lab->get_rigid_body_from_subframe(0,2);
/*	
	Point_Meshless=0;
	for(int i=0;i<lab->obj_number;i++) {
		Point_Meshless += lab->rigid_body[i]->POINT_NUMBER;
	}
	
	Re=insect_parameters.Re;
*/	
}

void Motion(int it, SOLVER& solver, MESH& mesh, REF_FRAME& lab, REF_FRAME& insect, REF_FRAME& wingplane_1, REF_FRAME& wingplane_2, REF_FRAME& wing_1, REF_FRAME& wing_2)
{
    
    my_flap_pattern.Update_Flap_Pattern(solver.Time, it, bodycentre, insect.angle);
	
    lab.Time=solver.Time;
    insect.Time=solver.Time;
    wingplane_1.Time=solver.Time;
    wingplane_2.Time=solver.Time;
    wing_1.Time=solver.Time;
    wing_2.Time=solver.Time;

    insect.update_all();
    wingplane_1.update_all();
    wingplane_2.update_all();
    wing_1.update_all();
    wing_2.update_all();
	
    wingplane_1.change_ref_frame();
    wingplane_2.change_ref_frame();
    insect.change_ref_frame();
    lab.change_ref_frame();
    
    
    mesh.UPDATE_MESHLESS(lab);
}


int main(int argc, char *argv[])
{
    std::cout << "Task Starting......" << std::endl << std::endl;
    
//    GPU_detection();
//    CPU_detection();
    
//    #pragma omp parallel
//    {
//	std::cout<<" Threads id: "<< omp_get_thread_num() << std::endl;
//    }

    
    timer          time;
    SOLVER         solver;
    Initialize_Ref_Frame(&solver, &lab, &insect, &wingplane_1, &wingplane_2, &wing_1, &wing_2);
    
    MESH_CARTESIAN cartesian;
    MESH_LESS      meshless(lab);
    MESH           mesh(cartesian, meshless);
    FLOW_FIELD     flow_field(mesh);
    LINKLIST       linklist;
    WRITE_FILE     write_file(lab,mesh,flow_field,solver,my_flap_pattern);
/*    
    linklist.Insert_point(3);
    linklist.Insert_point(6);
    linklist.Insert_point(7);
    linklist.Insert_point(30294);
    linklist.Insert_point(31286);
    
    linklist.Delete_point(cartesian,flow_field);
    
    std::cout<<linklist.LINKLISTPOINT.size()<<std::endl;
  
    std::cout<<" it contains: "<< std::endl;
    for(thrust::device_vector<LINKLIST_MEMBER>::iterator iter = linklist.LINKLISTPOINT.begin(); iter != linklist.LINKLISTPOINT.end(); iter++)  
    {
	std::cout << (static_cast<LINKLIST_MEMBER>(*iter)).Meshless_Ind << std::endl;
    }
    
    
    cartesian.POINTTYPE[1][30294]=4;
    linklist.Fresh_point(mesh,cartesian,flow_field);
*/

    mesh.SEARCH_TYPE(true,linklist,flow_field);
    
    for(int it=0;it<500;it++)
    {
	std::cout<<"Time = "<<it*0.002<<std::endl;
	solver.Time=it*0.002;
      
	mesh.UPDATE_POINTTYPE();
	
	for(int FSI=0; FSI<5; FSI++)
	{
	    mesh.UPDATE_POINTTYPE_IMPLICIT();
	
	    Motion(it,solver, mesh, lab, insect, wingplane_1, wingplane_2, wing_1, wing_2);
	
	    
	
	    mesh.SEARCH_TYPE(false,linklist,flow_field);
	    
	    mesh.TEST_FUNCTION();
	    linklist.Delete_point(mesh,flow_field);

	}

	write_file.OUTPUT_RESULT(it);
    }
    std::cout << "Elapsed time: " << time.milliseconds_elapsed() << "ms" << std::endl;

}
