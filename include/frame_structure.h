#ifndef FRAME_STRUCTURE_H
#define FRAME_STRUCTURE_H

#include <string>

#include "flap_pattern.h"
#include <cusp/array1d.h>
#include <cusp/array2d.h>

struct RIGID_OBJ
{
    int   IDENTITY;
    int   LEVEL;
    bool  POINTMODE;
    
    int   POINT_NUMBER, TETRAHEDRON_NUMBER, TRIANGLE_NUMBER;
    int   INNER_POINT_NUMBER, OUTER_POINT_NUMBER;
    
    cusp::array2d<double,cusp::device_memory,cusp::column_major>    XYZ;
    cusp::array2d<double,cusp::device_memory,cusp::column_major>    UVW;
    cusp::array2d<double,cusp::device_memory,cusp::column_major>    ACC;
    
    cusp::array2d<double,cusp::device_memory,cusp::column_major>    OUTER_NORMAL_VECTOR;
    cusp::array1d<double,cusp::device_memory>                       AREA;
    
    cusp::array2d<int,cusp::device_memory,cusp::column_major>       TETRAHEDRON;
    cusp::array2d<int,cusp::device_memory,cusp::column_major>       TRIANGLE;
    
    cusp::array1d<int,cusp::device_memory>                          INNERMARK;
    cusp::array1d<int,cusp::device_memory>                          OUTERMARK;
};


// A functor to update relative XYZ and angles from wing kinematics
// use template for this?
class AttitudeEvaluator 
{
  public:
	AttitudeEvaluator(FLAPPATTERN* _pt2object, double (FLAPPATTERN::*_pt2func)(double Time, int i))
	{
		pt2object = _pt2object;  pt2func = _pt2func;
	}
	AttitudeEvaluator() {}
	
	double operator()(double Time, int i) // call using operator
	{
		return (*pt2object.*pt2func)(Time, i);
	}
	
  private:
	double (FLAPPATTERN::*pt2func)(double Time, int i);   // pointer to member function
	FLAPPATTERN *pt2object;                  // pointer to object
};

class REF_FRAME
{
  public:
	double        Time;
  
	int           level;
	REF_FRAME     *super;
	int           subframe_number;
	REF_FRAME     **sub;
	double        local_inertia[3][3];
	
	double        position[3];            //X
	double        velocity[3];            //V
	double        acceleration[3];
	double        angle[3];
	double        angular_velocity[3];    //omega
	double        angular_acceleration[3];
	double        orientation[3][3];      //Rc
	double        force[3];               //F
	double        torque[3];              //T
	double        momentum[3];            //Mv
	double        inertia[3][3];          //I
	double        angular_momentum[3];    //L
	
	double        power;
	double        energy;
	
	double        position_old[3];            //X
	double        velocity_old[3];            //V
	double        acceleration_old[3];
	double        angle_old[3];
	double        angular_velocity_old[3];    //omega
	double        angular_acceleration_old[3];
	double        orientation_old[3][3];      //Rc
	double        force_old[3];               //F
	double        torque_old[3];              //T
	double        momentum_old[3];            // Mv
	double        inertia_old[3][3];         //I
	double        angular_momentum_old[3];    //L
	
	double        position_tmp[3];            //X
	double        velocity_tmp[3];            //V
	double        acceleration_tmp[3];
	double        angle_tmp[3];
	double        angular_velocity_tmp[3];    //omega
	double        angular_acceleration_tmp[3];
	double        orientation_tmp[3][3];      //Rc
	double        force_tmp[3];               //F
	double        torque_tmp[3];              //T
	double        momentum_tmp[3];            // Mv
	double        inertia_tmp[3][3];         //I
	double        angular_momentum_tmp[3];    //L
	
	int           obj_number;
	RIGID_OBJ     **rigid_body;

	RIGID_OBJ     *get_rigid_body_from_file(const char*);//input the file, allocate memory.
	RIGID_OBJ     *get_rigid_body_from_subframe(int,int);//input from other object, allocate memory.
	void          change_ref_frame();//modify the attributes of a body into the current frame;
	void          compute_orientation();
	void          compute_orientation_tmp();
	void          copy_new2old();
	
	AttitudeEvaluator compute_x; // functor version of -- double (*compute_x)(double Time, int i);
	AttitudeEvaluator compute_y;
	AttitudeEvaluator compute_z;
	AttitudeEvaluator compute_phi;
	AttitudeEvaluator compute_theta;
	AttitudeEvaluator compute_psi;
	
	void          update_all();
	
  private:
	void          memory_allocate_for_rigid_body(RIGID_OBJ*);
	double        compute_angular_velocity1(double Time);
	double        compute_angular_velocity2(double Time);
	double        compute_angular_velocity3(double Time);
	double        compute_angular_acceleration1(double Time);
	double        compute_angular_acceleration2(double Time);
	double        compute_angular_acceleration3(double Time);
};


#endif /*FRAME_STRUCTURE_H*/
