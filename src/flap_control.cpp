#include <fstream>
#include <iostream>

#define  PI 3.1415926535897932

#include "../include/flap_pattern.h"

void FLAPPATTERN :: Init_Flap_Controllers()
{
    X_controller.name="X_controller";
    X_controller.set=0.0;
    X_controller.Kp=-12.5*PI/9.0;
    X_controller.Ki=-10;
    X_controller.Kd=-50;
    X_controller.epsilon=0.005;
    X_controller.Init=0;
    X_controller.Max=5*PI/180.0;
    X_controller.Min=-5*PI/180.0;
    
    
    Y_controller.name="Y_controller";
    Y_controller.set=0.0;
    Y_controller.Kp=-25*PI/9.0;
    Y_controller.Ki=-0; 
    Y_controller.Kd=-20;
    Y_controller.epsilon=0.005;
    Y_controller.Init=0*PI/180;
    Y_controller.Max=5*PI/180; // modified from 5 to 20
    Y_controller.Min=-5*PI/180;
    
    
    Z_controller.name="Z_controller";
    Z_controller.set=8.0;
    Z_controller.Kp=400;
    Z_controller.Ki=100;
    Z_controller.Kd=200;
    Z_controller.epsilon=0.005;
    Z_controller.Init=260;
    Z_controller.Max=300;
    Z_controller.Min=220;
    
    
    V_controller.name="V_controller";
    V_controller.Kp=-25*PI/9.0; 
    V_controller.Ki=-0;
    V_controller.Kd=-2;
    V_controller.epsilon=0.002;
    V_controller.Init=0;
    V_controller.Max=5*PI/180;
    V_controller.Min=-5*PI/180;

    
    Orientation_controller.name="Orientation_controller";
    Orientation_controller.set=0*PI/180;
    Orientation_controller.Kp=1.5;
    Orientation_controller.Ki=0.01;
    Orientation_controller.Kd=4.5;
    Orientation_controller.epsilon=3*PI/180;
    Orientation_controller.Init=-1.5*PI/180;
    Orientation_controller.Max=2.0*PI/180;
    Orientation_controller.Min=-5.0*PI/180;
	//Orientation_controller.Max_Change=3.0*PI/180;
	

    Roll_controller.name="Roll_controller";
    Roll_controller.set=0.0*PI/180;
    Roll_controller.Kp=0.4;
    Roll_controller.Ki=0.1;
    Roll_controller.Kd=1.0;
    Roll_controller.epsilon=2.0*PI/180;
    Roll_controller.Init=0*PI/180;
    Roll_controller.Max=3.0*PI/180;
    Roll_controller.Min=-3.0*PI/180;
    
    
    Yaw_controller.name="Yaw_controller";
    Yaw_controller.set=0.0*PI/180;
    Yaw_controller.Kp=-0.3;
    Yaw_controller.Ki=0;
    Yaw_controller.Kd=-0.4;
    Yaw_controller.epsilon=1.0*PI/180;
    Yaw_controller.Init=0*PI/180;
    Yaw_controller.Max=3.0*PI/180;
    Yaw_controller.Min=-3.0*PI/180;
}

void FLAPPATTERN :: Update_Flap_Controllers(double body_centre[3], double body_angle[3], double body_centre_old[3])
{
//
//  current input structure
//  body_centre[3]: 0 --> x,     1 --> y,           2 --> z
//  body_angle[3]:  0 --> yaw,   1 --> -1 * roll,   2 --> pitch
//

	double x_aim = 8.0, y_aim = 8.0;
	double x_delta = x_aim - body_centre[0], y_delta = y_aim - body_centre[1];
	double current_dist = sqrt(x_delta*x_delta + y_delta*y_delta);
	
	double x_input = -(x_delta*cos(-body_angle[0])-y_delta*sin(-body_angle[0]));
	double y_input = -(x_delta*sin(-body_angle[0])+y_delta*cos(-body_angle[0]));
	double v_input = (body_centre[0]-body_centre_old[0])*sin(-body_angle[0])
					+(body_centre[1]-body_centre_old[1])*cos(-body_angle[0]);
	/*
	double yaw_input = atan2(x_delta, y_delta) + insect->angle[0];
	if(current_dist < 20*maxD) {
		if(fabs(yaw_input) > 2.0/3.0*PI) yaw_input = (yaw_input > 0) ? (yaw_input - PI) : (yaw_input + PI);
		if(fabs(v_input) > maxV || current_dist < maxD) yaw_input = 0;
		if(current_dist<maxD && fabs(v_input)<maxV) yaw_input = insect->angle[0];
	}
	*/
	double yaw_input   = body_angle[0]; // double yaw_input = insect->angle[0];
	double roll_input  = body_angle[1];
	double pitch_input = body_angle[2];
	
//  longitudinal position and velocity controller, corresponding to wingplane angle PSI -- Angle_*_Psi
	double maxD=0.1, maxV=0.1;
	if(current_dist > 20*maxD) maxV = 0.2;
	if(y_input>maxD) {
		if(v_input < -maxV) {
			V_controller.set=0.0; Angle_Next_Psi = V_controller.Activate_controller(v_input);
			Y_controller.Activate_controller(y_input);
		}
		else {
			V_controller.set=-maxV; Angle_Next_Psi = V_controller.Activate_controller(v_input);
			Y_controller.Activate_controller(y_input);
		}
	}
	else if(y_input<-maxD) {
		if(v_input > maxV) {
			V_controller.set=0.0; Angle_Next_Psi=V_controller.Activate_controller(v_input);
			Y_controller.Activate_controller(y_input);
		}
		else {
			V_controller.set=maxV; Angle_Next_Psi=V_controller.Activate_controller(v_input);
			Y_controller.Activate_controller(y_input);
		}
	}
	else if(fabs(y_input)<=maxD) {
		if(fabs(v_input) > maxV) {
			V_controller.set=0.0; Angle_Next_Psi=V_controller.Activate_controller(v_input);
			Y_controller.Activate_controller(y_input);
		}
		else {
			V_controller.set=0.0; V_controller.Activate_controller(v_input);
			Angle_Next_Psi=Y_controller.Activate_controller(y_input);
		}
	}
	else {
		std::cout<<"Mistake: beyond control plan....."<<std::endl;
	}

//  lateral position controller, corresponding to wingplane angle THETA -- Angle_*_Theta
//	int x_ctrl_interval = 20;
//	if(X_controller.record.size() % x_ctrl_interval==0) {
//		if(x_input>0.05) {
//			Angle_Next_Theta = X_controller.Max;
//			X_controller.record.push_back(1);
//		}
//		else if(x_input<-0.05) {
//			Angle_Next_Theta = X_controller.Min;
//			X_controller.record.push_back(1);
//		}
//		else {
//			Angle_Next_Theta = 0;
//			X_controller.record.push_back(-1);
//		}
//	}
//	else if(X_controller.record.size() % x_ctrl_interval < (x_ctrl_interval/2)) {
//		Angle_Next_Theta = Angle_Current_Theta;
//		X_controller.record.push_back(1);
//	}
//	else {
//		Angle_Next_Theta = 0;
//		X_controller.record.push_back(-1);
//	}


//  pitch angle controller, corresponding to wingplane mean position angle PHI -- Angle_*_Phi
	Angle_Next_Phi=Orientation_controller.Activate_controller(pitch_input) - 0.3*Angle_Next_Psi;

//  yaw angle controller, corresponding to change in wing angle of attack -- Angle_*_AoA_Delta
	Angle_Next_AoA_Delta=Yaw_controller.Activate_controller(yaw_input);
	
//  roll angle controller, corresponding to change in wing stroke amplitude -- Angle_*_Stroke_Delta
	Angle_Next_Stroke_Delta=Roll_controller.Activate_controller(roll_input);
	
//  Manually set controlling angles
//	Angle_Next_Phi = Angle_Current_Phi;
//	Angle_Next_Psi = -15.0*PI/180.0 - insect->angle[2];
	Angle_Next_Theta = 0.0;
//	Angle_Next_AoA_Delta = 0.0;
//	Angle_Next_Stroke_Delta = 0.0;
	
}

void FLAPPATTERN :: Update_Flap_Pattern(double mtime, int it, double body_centre[3], double body_angle[3])
{
	if(it%check_intvl_Angle==0 && modify_mark[int(it/check_intvl_Angle)]==0)
	{
	    modify_mark[int(it/check_intvl_Angle)]=1;
	    
		Angle_Previous_Psi = Angle_Current_Psi;
		Angle_Previous_Theta = Angle_Current_Theta;
		Angle_Previous_Phi = Angle_Current_Phi;

		Angle_Previous_AoA_Delta = Angle_Current_AoA_Delta;
		Angle_Previous_Stroke_Delta = Angle_Current_Stroke_Delta;
	    
	    Update_Flap_Controllers(body_centre, body_angle, lastbodycentre);
		
		Update_Flap_Kinematics(mtime,it);
	    
	    for(int i=0;i<3;i++)
	    	lastbodycentre[i] = body_centre[i];
	}

	Angle_Current_Phi   = wingplane1_Phi(mtime, 0);
	Angle_Current_Theta = wingplane1_Theta(mtime, 0);
	Angle_Current_Psi   = wingplane1_Psi(mtime, 0);

	if(fabs(cos(2*PI*mtime)) > 0.001)
		Angle_Current_Stroke_Delta = (wing1_Phi(mtime, 0) + wing2_Phi(mtime, 0))/cos(2*PI*mtime)/2.0;
	
	Angle_Current_AoA_Delta = Angle_Next_AoA_Delta;
}

int FLAPPATTERN :: Load_Flap_Pattern(std::ifstream &in_info)
{
	if(!in_info.is_open()) {
		std::cout<<"Unable to open Motion Record File!"<<std::endl;
		return 1;
	}
	
	// controller
	in_info>>Angle_Previous_Phi         
	       >>Angle_Current_Phi         
	       >>Angle_Next_Phi;
	in_info>>Angle_Previous_Theta       
	       >>Angle_Current_Theta
	       >>Angle_Next_Theta;
	in_info>>Angle_Previous_Psi
	       >>Angle_Current_Psi
	       >>Angle_Next_Psi;
	in_info>>Angle_Previous_Stroke_Delta
	       >>Angle_Current_Stroke_Delta
	       >>Angle_Next_Stroke_Delta;
	in_info>>Angle_Previous_AoA_Delta   
	       >>Angle_Current_AoA_Delta
	       >>Angle_Next_AoA_Delta;
	
	int count=7; double tmp[count];
	for(int i=0;i<count;i++) in_info>>tmp[i];
	X_controller.Evoke_controller(tmp);
	for(int i=0;i<count;i++) in_info>>tmp[i];
	Y_controller.Evoke_controller(tmp);
	for(int i=0;i<count;i++) in_info>>tmp[i];
	Z_controller.Evoke_controller(tmp);
	for(int i=0;i<count;i++) in_info>>tmp[i];
	V_controller.Evoke_controller(tmp);
	for(int i=0;i<count;i++) in_info>>tmp[i];
	Orientation_controller.Evoke_controller(tmp);
	for(int i=0;i<count;i++) in_info>>tmp[i];
	Roll_controller.Evoke_controller(tmp);
	for(int i=0;i<count;i++) in_info>>tmp[i];
	Yaw_controller.Evoke_controller(tmp);
	
	in_info>>lastbodycentre[0]>>lastbodycentre[1]>>lastbodycentre[2];
	
	return 0;
}

int FLAPPATTERN :: Record_Flap_Pattern(std::ofstream &out_info)
{
	if(!out_info.is_open()) {
		std::cout<<"Unable to open Motion Record File!"<<std::endl;
		return 1;
	}
	// controller
	out_info<<Angle_Previous_Phi<<" "
	        <<Angle_Current_Phi<<" "
	        <<Angle_Next_Phi<<std::endl;
	out_info<<Angle_Previous_Theta<<" "
	        <<Angle_Current_Theta<<" "
	        <<Angle_Next_Theta<<std::endl;
	out_info<<Angle_Previous_Psi<<" "
	        <<Angle_Current_Psi<<" "
	        <<Angle_Next_Psi<<std::endl;
	out_info<<Angle_Previous_Stroke_Delta<<" "
	        <<Angle_Current_Stroke_Delta<<" "
	        <<Angle_Next_Stroke_Delta<<std::endl;
	out_info<<Angle_Previous_AoA_Delta<<" "
	        <<Angle_Current_AoA_Delta<<" "
	        <<Angle_Next_AoA_Delta<<std::endl;
	
	int count=7; double tmp[count];
	X_controller.Record_controller(tmp); 
	for(int i=0;i<count;i++) out_info<<tmp[i]<<" "; out_info<<std::endl;
	Y_controller.Record_controller(tmp); 
	for(int i=0;i<count;i++) out_info<<tmp[i]<<" "; out_info<<std::endl;
	Z_controller.Record_controller(tmp); 
	for(int i=0;i<count;i++) out_info<<tmp[i]<<" "; out_info<<std::endl;
	V_controller.Record_controller(tmp); 
	for(int i=0;i<count;i++) out_info<<tmp[i]<<" "; out_info<<std::endl;
	Orientation_controller.Record_controller(tmp); 
	for(int i=0;i<count;i++) out_info<<tmp[i]<<" "; out_info<<std::endl;
	Roll_controller.Record_controller(tmp); 
	for(int i=0;i<count;i++) out_info<<tmp[i]<<" "; out_info<<std::endl;
	Yaw_controller.Record_controller(tmp); 
	for(int i=0;i<count;i++) out_info<<tmp[i]<<" "; out_info<<std::endl;
	
	for(int i=0;i<3;i++) out_info<<lastbodycentre[i]<<" "; 
	out_info<<std::endl;
	
	return 0;
}


int FLAPPATTERN :: Record_Flap_Pattern_History(double mtime,
								  /*for test*/ double wp_a[3], double w1_a[3], double w2_a[3], double freq)
{
	if(!flap_history.is_open()) {
		std::cout<<"Unable to open Flap History File!"<<std::endl;
		return 1;
	}
	
	flap_history<<mtime<<" "
	            <<w1_a[0]<<" "<<w1_a[1]<<" "<<w1_a[2]<<" "
	            <<w2_a[0]<<" "<<w2_a[1]<<" "<<w2_a[2]<<" "
	            <<freq<<" "
	            <<Angle_Current_Stroke_Delta<<" "
	            <<Angle_Current_AoA_Delta<<" "
	            <<Angle_Previous_Phi<<" "
	            <<Angle_Previous_Theta<<" "
	            <<Angle_Previous_Psi<<" "
	            <<wp_a[0]<<" "<<wp_a[1]<<" "<<wp_a[2]<<" "
	            <<std::endl;
		
	return 0;
}

