#ifndef FLAP_PATTERN_H
#define FLAP_PATTERN_H

#include <fstream>

#include "control_module.h"
#include "solver.h"

class FLAPPATTERN 
{
    public:
	/* ***** Controlling Angles ***** */
	double Angle_Current_Phi, Angle_Previous_Phi, Angle_Next_Phi;
	double Angle_Current_Theta, Angle_Previous_Theta, Angle_Next_Theta;
	double Angle_Current_Psi, Angle_Previous_Psi, Angle_Next_Psi;
	
	double Angle_Current_Stroke_Delta, Angle_Previous_Stroke_Delta, Angle_Next_Stroke_Delta;
	double Angle_Current_AoA_Delta, Angle_Previous_AoA_Delta, Angle_Next_AoA_Delta;
	/* ****************************** */
	
	PID_CONTROL_MODULE X_controller, Y_controller, Z_controller, V_controller;
	PID_CONTROL_MODULE Orientation_controller, Roll_controller, Yaw_controller;
		
	int check_intvl_Re;
	int check_intvl_Angle;
	short *modify_mark;
	
	double wing1_X(double Time, int i);
	double wing1_Y(double Time, int i);
	double wing1_Z(double Time, int i);
	double wing2_X(double Time, int i);
	double wing2_Y(double Time, int i);
	double wing2_Z(double Time, int i);

	double wingplane1_X(double Time, int i);
	double wingplane1_Y(double Time, int i);
	double wingplane1_Z(double Time, int i);
	double wingplane2_X(double Time, int i);
	double wingplane2_Y(double Time, int i);
	double wingplane2_Z(double Time, int i);

	double insect_X(double Time, int i);
	double insect_Y(double Time, int i);
	double insect_Z(double Time, int i);
	
	double wing1_Phi(double Time, int i);
	double wing1_Theta(double Time, int i);
	double wing1_Psi(double Time, int i);
	double wing2_Phi(double Time, int i);
	double wing2_Theta(double Time, int i);
	double wing2_Psi(double Time, int i);

	double wingplane1_Phi(double Time, int i);
	double wingplane1_Theta(double Time, int i);
	double wingplane1_Psi(double Time, int i);
	double wingplane2_Phi(double Time, int i);
	double wingplane2_Theta(double Time, int i);
	double wingplane2_Psi(double Time, int i);
	
	double insect_Phi(double Time, int i);
	double insect_Theta(double Time, int i);
	double insect_Psi(double Time, int i);
	
	void   Set_Flap_Pattern(double root[3], double body_centre[3]);
	void   Update_Flap_Pattern(double mtime, int it, double body_centre[3], double body_angle[3]);
	void   Update_Flap_Kinematics(double mtime,int it);
	void   Get_Current_Controlling_Angles(double mtime, double wingplane_angle[3], double wing1_angle[3], double wing2_angle[3]);
	
	int    Record_Flap_Pattern(std::ofstream &file);
	int    Load_Flap_Pattern(std::ifstream &file);
	int    Record_Flap_Pattern_History(double mtime,
					      /*for test*/ double wp_a[3], double w1_a[3], double w2_a[3], double freq);
	
	SOLVER *solver;
	double bodycentre[3];
	double lastbodycentre[3];
	
	FLAPPATTERN() 
	{
		modify_mark = new short [1000];
		//flap_history.open("flap_history");
	}
	~FLAPPATTERN() 
	{ 
		delete[] modify_mark;
		//flap_history.close();
	}
	
  private:
	double wingplane_root[3];
	double mid_aoa, stroke_amplitude;

	// wing plane motion parameters
	double cc_wingplane_phi[6];
	double cc_wingplane_theta[6];
	double cc_wingplane_psi[6];

	// wing motion parameters
	double cc_wing_phi[6], cc_wing_phi_ini[7];
	double cc_wing_psi_ini[7];
	double cc_wing_psi_a1[8];
	double cc_wing_psi_a2[8];
	double cc_wing_psi_b1[9];
	double cc_wing_psi_b2[9];
	double cc_wing_psi_c1[9];
	double cc_wing_psi_c2[9];
	
	
	std::ofstream flap_history;
	
	void   Init_Flap_Kinematics();
	void   Init_Flap_Controllers();
	void   Update_Flap_Controllers(double body_centre[3], double body_angle[3], double body_centre_old[3]);
	void   Cal_Poly_Connection(double x1, int order1, double x2, int order2, double x_ctrl, int ctrl_number, double *right);
};

#endif /* FLAP_PATTERN_H */
