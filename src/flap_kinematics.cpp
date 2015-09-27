#include <iostream>
#include <math.h>

using namespace std;

#define  PI 3.1415926535897932

#include "../include/basis.h"
#include "../include/flap_pattern.h"


void FLAPPATTERN :: Set_Flap_Pattern(double root[3], double body_centre[3]) 
{
	mid_aoa = 45/180.0*PI;
	stroke_amplitude = 7/18.0*PI;
	
	for(int i=0;i<3;i++) {
		wingplane_root[i] = root[i];
		bodycentre[i]     = body_centre[i];
		lastbodycentre[i] = body_centre[i];
	}
	
	Init_Flap_Kinematics();
	
	check_intvl_Re    = 500;
	check_intvl_Angle = 500;
	
	Init_Flap_Controllers();
	
	Angle_Current_Phi    = wingplane1_Phi(0, 0);
	Angle_Previous_Phi   = wingplane1_Phi(0, 0);
	Angle_Next_Phi       = wingplane1_Phi(0, 0);
	
	Angle_Current_Theta  = wingplane1_Theta(0, 0);
	Angle_Previous_Theta = wingplane1_Theta(0, 0);
	Angle_Next_Theta     = wingplane1_Theta(0, 0);
	
	Angle_Current_Psi    = wingplane1_Psi(0, 0);
	Angle_Previous_Psi   = wingplane1_Psi(0, 0);
	Angle_Next_Psi       = wingplane1_Psi(0, 0);
	
	Angle_Current_Stroke_Delta  = 0;
	Angle_Previous_Stroke_Delta = 0;
	Angle_Next_Stroke_Delta     = 0;
	
	Angle_Current_AoA_Delta     = 0;
	Angle_Previous_AoA_Delta    = 0;
	Angle_Next_AoA_Delta        = 0;
	
	for(int i=0;i<1000;i++) modify_mark[i]=0;
}

double FLAPPATTERN :: wing1_X(double Time, int i) 
{
	if(i==0)
		return 0.27;
	else if(i==1)
		return 0;
	else if(i==2)
		return 0;
	else {
		cout<<"Mistake in function wing1_X!"<<endl;
		return 0;
	}
}
double FLAPPATTERN :: wing1_Y(double Time, int i) 
{
	return 0;
}
double FLAPPATTERN :: wing1_Z(double Time, int i) 
{
	return 0;
}
double FLAPPATTERN :: wing2_X(double Time, int i) 
{
	if(i==0)
		return -0.27;
	else if(i==1)
		return 0;
	else if(i==2)
		return 0;
	else {
		cout<<"Mistake in function wing2_X!"<<endl;
		return 0;
	}
}
double FLAPPATTERN :: wing2_Y(double Time, int i) 
{
	return 0;
}
double FLAPPATTERN :: wing2_Z(double Time, int i) 
{
	return 0;
}

double FLAPPATTERN :: wingplane1_X(double Time, int i) 
{
	if(i==0)
		return wingplane_root[0];
	else if(i==1)
		return 0;
	else if(i==2)
		return 0;
	else {
		cout<<"Mistake in function wingplane1_X!"<<endl;
		return 0;
	}
}
double FLAPPATTERN :: wingplane1_Y(double Time, int i) 
{
    if(i==0)
		return wingplane_root[1];
	else if(i==1)
		return 0;
	else if(i==2)
		return 0;
	else {
		cout<<"Mistake in function wingplane1_Y!"<<endl;
		return 0;
	}
}
double FLAPPATTERN :: wingplane1_Z(double Time, int i) 
{
    if(i==0)
		return wingplane_root[2];
	else if(i==1)
		return 0;
	else if(i==2)
		return 0;
	else {
		cout<<"Mistake in function wingplane1_Z!"<<endl;
		return 0;
	}
}
double FLAPPATTERN :: wingplane2_X(double Time, int i) 
{
	if(i==0)
		return -wingplane_root[0];
	else if(i==1)
		return 0;
	else if(i==2)
		return 0;
	else {
		cout<<"Mistake in function wingplane2_X!"<<endl;
		return 0;
	}
}
double FLAPPATTERN :: wingplane2_Y(double Time, int i) 
{
    if(i==0)
		return wingplane_root[1];
	else if(i==1)
		return 0;
	else if(i==2)
		return 0;
	else {
		cout<<"Mistake in function wingplane2_Y!"<<endl;
		return 0;
	}
}
double FLAPPATTERN :: wingplane2_Z(double Time, int i) 
{
    if(i==0)
		return wingplane_root[2];
	else if(i==1)
		return 0;
	else if(i==2)
		return 0;
	else {
		cout<<"Mistake in function wingplane2_Z!"<<endl;
		return 0;
	}
}

double FLAPPATTERN :: insect_X(double Time, int i)
{
    if(i==0)
	return bodycentre[0];
    else if(i==1)
	return 0;
    else if(i==2)
	return 0;
    else{
	cout<<"Mistake in function insect_X!"<<endl;
	return 0;
    }
}

double FLAPPATTERN :: insect_Y(double Time, int i)
{
    if(i==0)
	return bodycentre[1];
    else if(i==1)
	return 0;
    else if(i==2)
	return 0;
    else{
	cout<<"Mistake in function insect_Y!"<<endl;
	return 0;
    }
}

double FLAPPATTERN :: insect_Z(double Time, int i)
{
    if(i==0)
	return bodycentre[2];
    else if(i==1)
	return 0;
    else if(i==2)
	return 0;
    else{
	cout<<"Mistake in function insect_Z!"<<endl;
	return 0;
    }
}


double FLAPPATTERN :: wing1_Phi(double Time, int i) 
{
/*
	double phi_a=-6.1110694731930507e+03;
	double phi_b=6.1941704562879631e+03;
	double phi_c=-2.3982617331604333e+03;
	double phi_d=3.8610599078009807e+02;
	double phi_e=0.0000000000000000e+00;
	double phi_f=0.0000000000000000e+00;
	double phi_g=-1.2217304763960324e+00;
*/
	double tmp_wing_phi;
	int order;
	
	if(Time<0.25) {
		tmp_wing_phi = 0;
		order = (sizeof(cc_wing_phi_ini) / sizeof(cc_wing_phi_ini[0]));
		if(i==0) {
			for(int k=0;k<order;k++)
				tmp_wing_phi += cc_wing_phi_ini[k]*pow(Time,k);
			return tmp_wing_phi;
			//return phi_a*pow(Time,6)+phi_b*pow(Time,5)+phi_c*pow(Time,4)+phi_d*pow(Time,3)+phi_e*pow(Time,2)+phi_f*Time+phi_g;
	    }
		if(i==1) {
			for(int k=1;k<order;k++)
				tmp_wing_phi += k*cc_wing_phi_ini[k]*pow(Time,k-1);
			return tmp_wing_phi;
			//return 6*phi_a*pow(Time,5)+5*phi_b*pow(Time,4)+4*phi_c*pow(Time,3)+3*phi_d*pow(Time,2)+2*phi_e*Time+phi_f;
	    }
		if(i==2) {
			for(int k=2;k<order;k++)
				tmp_wing_phi += k*(k-1)*cc_wing_phi_ini[k]*pow(Time,k-2);
			return tmp_wing_phi;
			//return 30*phi_a*pow(Time,4)+20*phi_b*pow(Time,3)+12*phi_c*pow(Time,2)+6*phi_d*Time+2*phi_e;
	    }
		else {
			cout<<"Mistake in function wing_Phi!"<<endl;
			return 0;
		}
	}
	
	double Time_p = Time-floor(Time);
	if(Time_p>0.9999) Time_p-=1.0;
	
	double delta_stroke_amplitude[3] = {0};
	if(Time>0.9999) {
		order = (sizeof(cc_wing_phi) / sizeof(cc_wing_phi[0]));
		if(Time_p<0.5) {
			for(int k=0;k<order;k++) 
				delta_stroke_amplitude[0] += cc_wing_phi[k]*pow(Time_p,k);
			for(int k=1;k<order;k++) 
				delta_stroke_amplitude[1] += k*cc_wing_phi[k]*pow(Time_p,k-1);
			for(int k=2;k<order;k++) 
				delta_stroke_amplitude[2] += k*(k-1)*cc_wing_phi[k]*pow(Time_p,k-2);
		}
		else {
			for(int k=0;k<order;k++) 
				delta_stroke_amplitude[0] += cc_wing_phi[k]*pow(0.5,k);
			delta_stroke_amplitude[1] = 0;
			delta_stroke_amplitude[2] = 0;
		}
	}
	
	if(i==0) {
		return -(stroke_amplitude - delta_stroke_amplitude[0])*cos(2*PI*Time_p);
	}
	else if(i==1) {
		return stroke_amplitude*2*PI*sin(2*PI*Time_p)
		    + ( - delta_stroke_amplitude[0]*2*PI*sin(2*PI*Time_p) 
		        + delta_stroke_amplitude[1]*cos(2*PI*Time_p) );
	}
	else if(i==2) {
		return stroke_amplitude*4*PI*PI*cos(2*PI*Time_p) 
			+ ( - delta_stroke_amplitude[0]*4*PI*PI*cos(2*PI*Time_p)
			    - delta_stroke_amplitude[1]*4*PI*sin(2*PI*Time_p) 
			    + delta_stroke_amplitude[2]*cos(2*PI*Time_p) );
	}
	else {
		cout<<"Mistake in function wing1_Phi!"<<endl;
		return 0;
	}
}

double FLAPPATTERN :: wing1_Theta(double Time, int i) 
{
	if(i==0)
		return 0.0*PI/180.0;
	else if(i==1)
		return 0;
	else if(i==2)
		return 0;
	else {
		cout<<"Mistake in function wingplane1_X!"<<endl;
		return 0;
	}
}
double FLAPPATTERN :: wing1_Psi(double Time, int i) {
/*
	double psi_a=2.0263498527524524e+04;
	double psi_b=-1.9030909357787743e+04;
	double psi_c=6.3192344944722654e+03;
	double psi_d=-7.5725943570633979e+02;
	double psi_e=0.0000000000000000e+00;
	double psi_f=0.0000000000000000e+00;
	double psi_g=1.5707963267948974e+00;
*/
	double temp;
	int psi_a_order, psi_b_order, psi_c_order;
	
	if(Time<0.25) {
		temp = 0;
		psi_a_order = (sizeof(cc_wing_psi_ini) / sizeof(cc_wing_psi_ini[0]));
		if(i==0) {
			for(int k=0;k<psi_a_order;k++)
				temp += cc_wing_psi_ini[k]*pow(Time,k);
			return temp;
			//return psi_a*pow(Time,6)+psi_b*pow(Time,5)+psi_c*pow(Time,4)+psi_d*pow(Time,3)+psi_e*pow(Time,2)+psi_f*Time+psi_g;
	    }
		else if(i==1) {
			for(int k=1;k<psi_a_order;k++)
				temp += k*cc_wing_psi_ini[k]*pow(Time,k-1);
			return temp;
			//return 6*psi_a*pow(Time,5)+5*psi_b*pow(Time,4)+4*psi_c*pow(Time,3)+3*psi_d*pow(Time,2)+2*psi_e*Time+psi_f;
	    }
		else if(i==2) {
			for(int k=2;k<psi_a_order;k++)
				temp += k*(k-1)*cc_wing_psi_ini[k]*pow(Time,k-2);
			return temp;
			//return 30*psi_a*pow(Time,4)+20*psi_b*pow(Time,3)+12*psi_c*pow(Time,2)+6*psi_d*Time+2*psi_e;
	    }
		else {
			cout<<"Mistake in function wing_Psi!"<<endl;
			return 0;
		}	
	}
	
	double Time_p = Time-floor(Time);
	if(Time_p>0.9999) Time_p -= 1.0;
	
	temp = 0;
	psi_a_order = (sizeof(cc_wing_psi_a1) / sizeof(cc_wing_psi_a1[0]));
	psi_b_order = (sizeof(cc_wing_psi_b1) / sizeof(cc_wing_psi_b1[0]));
	psi_c_order = (sizeof(cc_wing_psi_c1) / sizeof(cc_wing_psi_c1[0]));
	if(i==0) {
		// return PI/2.0-Angle_mid_aoa*sin(2*PI*Time);
		if(Time_p<=0.1) {
			for(int k=0;k<psi_a_order;k++) 
				temp += cc_wing_psi_a1[k]*pow(Time_p,k);
			return temp;
	    }
		else if(0.1<Time_p&&Time_p<=0.4) {
		    return (mid_aoa + Angle_Current_AoA_Delta)*sin(2*PI*Time_p) + PI/2.0*(1-sin(2*PI*Time_p));
	    }
	    else if(0.4<Time_p&&Time_p<=0.6) {
			for(int k=0;k<psi_b_order;k++) 
				temp += cc_wing_psi_b1[k]*pow(Time_p,k);
			return temp;
	    }  
	    else if(0.6<Time_p&&Time_p<=0.9) {
			return (mid_aoa - Angle_Current_AoA_Delta)*sin(2*PI*Time_p) + PI/2.0*(1-sin(2*PI*Time_p));
	    }
	    else {
			temp = 0;
			for(int k=0;k<psi_c_order;k++) 
				temp += cc_wing_psi_c1[k]*pow(Time_p,k);
			return temp;
	    }
	}
	else if(i==1) {
		// return -2*PI*Angle_mid_aoa*cos(2*PI*Time);
		if(Time_p<=0.1) {
			for(int k=1;k<psi_a_order;k++) 
				temp += k*cc_wing_psi_a1[k]*pow(Time_p,k-1);
			return temp;
	    }
		else if(0.1<Time_p&&Time_p<=0.4) {
		    return 2.0*PI*(mid_aoa + Angle_Current_AoA_Delta)*cos(2*PI*Time_p) - PI*PI*cos(2*PI*Time_p);
	    }
	    else if(0.4<Time_p&&Time_p<=0.6) {
			for(int k=1;k<psi_b_order;k++) 
				temp += k*cc_wing_psi_b1[k]*pow(Time_p,k-1);
			return temp;
	    }  
	    else if(0.6<Time_p&&Time_p<=0.9) {
			return 2.0*PI*(mid_aoa - Angle_Current_AoA_Delta)*cos(2*PI*Time_p) - PI*PI*cos(2*PI*Time_p);
	    }
	    else {
			for(int k=1;k<psi_c_order;k++) 
				temp += k*cc_wing_psi_c1[k]*pow(Time_p,k-1);
			return temp;
	    }
	}
	else if(i==2) {
		// return 4*PI*PI*Angle_mid_aoa*sin(2*PI*Time);
		if(Time_p<=0.1) {
			for(int k=2;k<psi_a_order;k++) 
				temp += k*(k-1)*cc_wing_psi_a1[k]*pow(Time_p,k-2);
			return temp;
	    }
		else if(0.1<Time_p&&Time_p<=0.4) {
		    return -4.0*PI*PI*(mid_aoa + Angle_Current_AoA_Delta)*sin(2*PI*Time_p) + 2*PI*PI*PI*sin(2*PI*Time_p);
	    }
	    else if(0.4<Time_p&&Time_p<=0.6) {
			for(int k=2;k<psi_b_order;k++) 
				temp += k*(k-1)*cc_wing_psi_b1[k]*pow(Time_p,k-2);
			return temp;
	    }  
	    else if(0.6<Time_p&&Time_p<=0.9) {
			return -4.0*PI*PI*(mid_aoa - Angle_Current_AoA_Delta)*sin(2*PI*Time_p) + 2*PI*PI*PI*sin(2*PI*Time_p);
		}
	    else {
			for(int k=2;k<psi_c_order;k++) 
				temp += k*(k-1)*cc_wing_psi_c1[k]*pow(Time_p,k-2);
			return temp;
		}
	}
	else {
		cout<<"Mistake in function wing1_Psi!"<<endl;
		return 0;
	}
}

double FLAPPATTERN :: wing2_Phi(double Time, int i)
{
	if(Time<0.9999)
		return -wing1_Phi(Time,i);
	
	double Time_p = Time-floor(Time);
	if(Time_p>0.9999) Time_p-=1.0;
	
	double delta_stroke_amplitude[3] = {0};
	int order = (sizeof(cc_wing_phi) / sizeof(cc_wing_phi[0]));
	if(Time_p<0.5) {
		for(int k=0;k<order;k++) delta_stroke_amplitude[0] += cc_wing_phi[k]*pow(Time_p,k);
		for(int k=1;k<order;k++) delta_stroke_amplitude[1] += k*cc_wing_phi[k]*pow(Time_p,k-1);
		for(int k=2;k<order;k++) delta_stroke_amplitude[2] += k*(k-1)*cc_wing_phi[k]*pow(Time_p,k-2);
	}
	else {
		for(int k=0;k<order;k++) delta_stroke_amplitude[0] += cc_wing_phi[k]*pow(0.5,k);
		delta_stroke_amplitude[1] = 0;
		delta_stroke_amplitude[2] = 0;
	}
	
	if(i==0) {
		return (stroke_amplitude + delta_stroke_amplitude[0])*cos(2*PI*Time_p);
	}
	else if(i==1) {
		return -stroke_amplitude*2*PI*sin(2*PI*Time_p) 
		    + ( - delta_stroke_amplitude[0]*2*PI*sin(2*PI*Time_p)
		        + delta_stroke_amplitude[1]*cos(2*PI*Time_p) );
	}
	else if(i==2) {
		return -stroke_amplitude*4*PI*PI*cos(2*PI*Time_p) 
			+ ( - delta_stroke_amplitude[0]*4*PI*PI*cos(2*PI*Time_p)
			    - delta_stroke_amplitude[1]*4*PI*sin(2*PI*Time_p)
			    + delta_stroke_amplitude[2]*cos(2*PI*Time_p) );
	}
	else {
		cout<<"Mistake in function wing2_Phi!"<<endl;
		return 0;
	}
}
double FLAPPATTERN :: wing2_Theta(double Time, int i)
{
	return -wing1_Theta(Time,i);
}
double FLAPPATTERN :: wing2_Psi(double Time, int i)
{	
	if(Time<0.9999)
		return wing1_Psi(Time,i);
	
	double Time_p = Time-floor(Time);
	if(Time_p>0.9999) Time_p-=1.0;
	int psi_a_order = (sizeof(cc_wing_psi_a2) / sizeof(cc_wing_psi_a2[0]));
	int psi_b_order = (sizeof(cc_wing_psi_b2) / sizeof(cc_wing_psi_b2[0]));
	int psi_c_order = (sizeof(cc_wing_psi_c2) / sizeof(cc_wing_psi_c2[0]));
	double temp = 0;
	
	if(i==0) {
		// return PI/2.0-Angle_mid_aoa*sin(2*PI*Time);
		if(Time_p<=0.1) {
			for(int k=0;k<psi_a_order;k++) 
				temp += cc_wing_psi_a2[k]*pow(Time_p,k);
			return temp;
	    }
		else if(0.1<Time_p&&Time_p<=0.4) {
		    return (mid_aoa - Angle_Current_AoA_Delta)*sin(2*PI*Time_p) + PI/2.0*(1-sin(2*PI*Time_p));
	    }
	    else if(0.4<Time_p&&Time_p<=0.6) {
			for(int k=0;k<psi_b_order;k++) 
				temp += cc_wing_psi_b2[k]*pow(Time_p,k);
			return temp;
	    }  
	    else if(0.6<Time_p&&Time_p<=0.9) {
			return (mid_aoa + Angle_Current_AoA_Delta)*sin(2*PI*Time_p) + PI/2.0*(1-sin(2*PI*Time_p));
	    }
	    else {
			for(int k=0;k<psi_c_order;k++) 
				temp += cc_wing_psi_c2[k]*pow(Time_p,k);
			return temp;
	    }
	}
	else if(i==1) {
		// return -2*PI*Angle_mid_aoa*cos(2*PI*Time);
		if(Time_p<=0.1) {
			for(int k=1;k<psi_a_order;k++) 
				temp += k*cc_wing_psi_a2[k]*pow(Time_p,k-1);
			return temp;
	    }
		else if(0.1<Time_p&&Time_p<=0.4) {
		    return 2.0*PI*(mid_aoa - Angle_Current_AoA_Delta)*cos(2*PI*Time_p) - PI*PI*cos(2*PI*Time_p);
	    }
	    else if(0.4<Time_p&&Time_p<=0.6) {
			for(int k=1;k<psi_b_order;k++) 
				temp += k*cc_wing_psi_b2[k]*pow(Time_p,k-1);
			return temp;
	    }  
	    else if(0.6<Time_p&&Time_p<=0.9) {
			return 2.0*PI*(mid_aoa + Angle_Current_AoA_Delta)*cos(2*PI*Time_p) - PI*PI*cos(2*PI*Time_p);
	    }
	    else {
			for(int k=1;k<psi_c_order;k++) 
				temp += k*cc_wing_psi_c2[k]*pow(Time_p,k-1);
			return temp;
	    }
	}
	else if(i==2) {
		// return 4*PI*PI*Angle_mid_aoa*sin(2*PI*Time);
		if(Time_p<=0.1) {
			for(int k=2;k<psi_a_order;k++) 
				temp += k*(k-1)*cc_wing_psi_a2[k]*pow(Time_p,k-2);
			return temp;
	    }
		else if(0.1<Time_p&&Time_p<=0.4) {
		    return -4.0*PI*PI*(mid_aoa - Angle_Current_AoA_Delta)*sin(2*PI*Time_p) + 2*PI*PI*PI*sin(2*PI*Time_p);
	    }
	    else if(0.4<Time_p&&Time_p<=0.6) {
			for(int k=2;k<psi_b_order;k++) 
				temp += k*(k-1)*cc_wing_psi_b2[k]*pow(Time_p,k-2);
			return temp;
	    }
	    else if(0.6<Time_p&&Time_p<=0.9) {
			return -4.0*PI*PI*(mid_aoa + Angle_Current_AoA_Delta)*sin(2*PI*Time_p) + 2*PI*PI*PI*sin(2*PI*Time_p);
		}
	    else {
			for(int k=2;k<psi_c_order;k++) 
				temp += k*(k-1)*cc_wing_psi_c2[k]*pow(Time_p,k-2);
			return temp;
		}
	}
	else {
		cout<<"Mistake in function wing2_Psi!"<<endl;
		return 0;
	}
}

double FLAPPATTERN :: wingplane1_Phi(double Time, int i) {
	int order = (sizeof(cc_wingplane_phi) / sizeof(cc_wingplane_phi[0]));
	double temp=0;
	
	if(i==0) {
		if(Time<=1.0)
			return -1.5*PI/180.0;
		else {
			for(int k=0;k<order;k++) 
				temp += cc_wingplane_phi[k]*pow(Time,k);
			return temp;
		}
    }
    else if(i==1) {
		if(Time<=1.0)
			return 0;
		else {
			for(int k=1;k<order;k++) 
				temp += k*cc_wingplane_phi[k]*pow(Time,k-1);
			return temp;
		}
    }
    else if(i==2) {
		if(Time<=1.0)
			return 0;
		else {
			for(int k=2;k<order;k++) 
				temp += k*(k-1)*cc_wingplane_phi[k]*pow(Time,k-2);
			return temp;
		}
    }
    else {
		cout<<"Mistake takes place in Function wingplane_Psi."<<endl;
		return 0;
    }
}
double FLAPPATTERN :: wingplane1_Theta(double Time, int i) {
	int order = (sizeof(cc_wingplane_theta) / sizeof(cc_wingplane_theta[0]));
	double temp=0;
	
	if(i==0) {
		if(Time<=1.0)
			return 0;
		else {
			for(int k=0;k<order;k++) 
				temp += cc_wingplane_theta[k]*pow(Time,k);
			return temp;
		}
	}
	else if(i==1) {
		if(Time<=1.0)
			return 0;
		else {
			for(int k=1;k<order;k++) 
				temp += k*cc_wingplane_theta[k]*pow(Time,k-1);
			return temp;
		}
	}
	else if(i==2) {
		if(Time<=1.0)
			return 0;
		else {
			for(int k=2;k<order;k++) 
				temp += k*(k-1)*cc_wingplane_theta[k]*pow(Time,k-2);
			return temp;
		}
	}
	else {
		cout<<"Mistake takes place in Function wingplane_Theta."<<endl;
		return 0;
	}
}
double FLAPPATTERN :: wingplane1_Psi(double Time, int i) {
	int order = (sizeof(cc_wingplane_psi) / sizeof(cc_wingplane_psi[0]));
	double temp=0;
	
	if(i==0) {
		if(Time<=1.0)
			return -0*PI/180.0;
		else {
			for(int k=0;k<order;k++) 
				temp += cc_wingplane_psi[k]*pow(Time,k);
			return temp;
		}
	}
	else if(i==1) {
		if(Time<=1.0)
			return 0;
		else {
			for(int k=1;k<order;k++) 
				temp += k*cc_wingplane_psi[k]*pow(Time,k-1);
			return temp;
		}
	}
	else if(i==2) {
		if(Time<=1.0)
			return 0;
		else {
			for(int k=2;k<order;k++) 
				temp += k*(k-1)*cc_wingplane_psi[k]*pow(Time,k-2);
			return temp;
		}
	}
	else{
		cout<<"Mistake takes place in Function wingplane_Psi."<<endl;
		return 0;
	}
}

double FLAPPATTERN :: wingplane2_Phi(double Time, int i) {
	return -wingplane1_Phi(Time,i);
}
double FLAPPATTERN :: wingplane2_Theta(double Time, int i) {
	return wingplane1_Theta(Time,i);
}
double FLAPPATTERN :: wingplane2_Psi(double Time, int i) {
    return wingplane1_Psi(Time,i);
}

double FLAPPATTERN :: insect_Phi(double Time, int i) {
    
    if( i == 0 ) {
	return solver->angle[0];
    }
    else if( i == 1 ) {
	return solver->angular_velocity[0];
    }
    else if( i == 2 ) {
	return solver->angular_acceleration[0];
    }
    else {
	cout<<"Mistake takes place in Function insect_Phi."<<endl;
	return 0;
    }
}

double FLAPPATTERN :: insect_Theta(double Time, int i) {
    
    if( i == 0 ) {
	return solver->angle[1];
    }
    else if( i == 1 ) {
	return solver->angular_velocity[1];
    }
    else if( i == 2 ) {
	return solver->angular_acceleration[1];
    }
    else {
	cout<<"Mistake takes place in Function insect_Theta."<<endl;
	return 0;
    }
}

double FLAPPATTERN :: insect_Psi(double Time, int i) {
    
    if( i == 0 ) {
	return solver->angle[2];
    }
    else if( i == 1 ) {
	return solver->angular_velocity[2];
    }
    else if( i == 2 ) {
	return solver->angular_acceleration[2];
    }
    else {
	cout<<"Mistake takes place in Function insect_Phi."<<endl;
	return 0;
    }
}


void FLAPPATTERN :: Init_Flap_Kinematics() {
	double x1, x2;
	
	// Wing_Phi Initialization
	x1=0; x2=0.25;
	cc_wing_phi_ini[0]=-stroke_amplitude;
	cc_wing_phi_ini[1]=0;
	cc_wing_phi_ini[2]=0;
	cc_wing_phi_ini[3]=-stroke_amplitude*cos(2*PI*x2);
	cc_wing_phi_ini[4]=2*PI*stroke_amplitude*sin(2*PI*x2);
	cc_wing_phi_ini[5]=4*PI*PI*stroke_amplitude*cos(2*PI*x2);
	cc_wing_phi_ini[6]=-8*PI*PI*PI*stroke_amplitude*sin(2*PI*x2);
	Cal_Poly_Connection(x1, 3, x2, 4, 0, 0, cc_wing_phi_ini);
	
	
	// Wing_Psi Initialization
	// 		start section
	x1=0; x2=0.25;
	cc_wing_psi_ini[0] = PI/2.0;
	cc_wing_psi_ini[1] = 0;
	cc_wing_psi_ini[2] = 0;
	cc_wing_psi_ini[3] = mid_aoa*sin(2*PI*x2) + PI/2.0*(1-sin(2*PI*x2));
	cc_wing_psi_ini[4] = 2*PI*mid_aoa*cos(2*PI*x2) - PI*PI*cos(2*PI*x2);
	cc_wing_psi_ini[5] = -4*PI*PI*mid_aoa*sin(2*PI*x2) + 2*PI*PI*PI*sin(2*PI*x2);
	cc_wing_psi_ini[6] = -8*PI*PI*PI*mid_aoa*cos(2*PI*x2) + 4*PI*PI*PI*PI*cos(2*PI*x2);
	Cal_Poly_Connection(x1, 3, x2, 4, 0, 0, cc_wing_psi_ini);
	for(int i=0;i<7;i++) { cc_wing_psi_a1[i]=0; cc_wing_psi_a2[i]=0; }

	//		middle section (downstroke to upstroke)
	x1=0.4;x2=0.6;
	cc_wing_psi_b1[0] = mid_aoa*sin(2*PI*x1) + PI/2.0*(1-sin(2*PI*x1));
	cc_wing_psi_b1[1] = 2.0*PI*mid_aoa*cos(2*PI*x1) - PI*PI*cos(2*PI*x1);
	cc_wing_psi_b1[2] = -4.0*PI*PI*mid_aoa*sin(2*PI*x1) + 2*PI*PI*PI*sin(2*PI*x1);
	cc_wing_psi_b1[3] = -8.0*PI*PI*PI*mid_aoa*cos(2*PI*x1) + 4*PI*PI*PI*PI*cos(2*PI*x1);
	cc_wing_psi_b1[4] = mid_aoa*sin(2*PI*x2) + PI/2.0*(1-sin(2*PI*x2));
	cc_wing_psi_b1[5] = 2.0*PI*mid_aoa*cos(2*PI*x2) - PI*PI*cos(2*PI*x2);
	cc_wing_psi_b1[6] = -4.0*PI*PI*mid_aoa*sin(2*PI*x2) + 2*PI*PI*PI*sin(2*PI*x2);
	cc_wing_psi_b1[7] = -8.0*PI*PI*PI*mid_aoa*cos(2*PI*x2) + 4*PI*PI*PI*PI*cos(2*PI*x2);
	cc_wing_psi_b1[8] = PI/2.0; // As 1:1 stroke is needed, let curve pass (t=0.5, AoA=PI/2)
	Cal_Poly_Connection(x1, 4, x2, 4, 0.5, 1, cc_wing_psi_b1);
	for(int i=0;i<9;i++) cc_wing_psi_b2[i]=cc_wing_psi_b1[i];

	//		end section (upstroke to downstroke)
	x1=0.9;x2=1.1;
	cc_wing_psi_c1[0] = mid_aoa*sin(2*PI*x1) + PI/2.0*(1-sin(2*PI*x1));
	cc_wing_psi_c1[1] = 2.0*PI*mid_aoa*cos(2*PI*x1) - PI*PI*cos(2*PI*x1);
	cc_wing_psi_c1[2] = -4.0*PI*PI*mid_aoa*sin(2*PI*x1) + 2*PI*PI*PI*sin(2*PI*x1);
	cc_wing_psi_c1[3] = -8.0*PI*PI*PI*mid_aoa*cos(2*PI*x1) + 4*PI*PI*PI*PI*cos(2*PI*x1);
	cc_wing_psi_c1[4] = mid_aoa*sin(2*PI*x2) + PI/2.0*(1-sin(2*PI*x2));
	cc_wing_psi_c1[5] = 2.0*PI*mid_aoa*cos(2*PI*x2) - PI*PI*cos(2*PI*x2);
	cc_wing_psi_c1[6] = -4.0*PI*PI*mid_aoa*sin(2*PI*x2) + 2*PI*PI*PI*sin(2*PI*x2);
	cc_wing_psi_c1[7] = -8.0*PI*PI*PI*mid_aoa*cos(2*PI*x2) + 4*PI*PI*PI*PI*cos(2*PI*x2);
	cc_wing_psi_c1[8] = PI/2.0; // let curve pass (t=1.0, AoA=PI/2)
	Cal_Poly_Connection(x1, 4, x2, 4, 1.0, 1, cc_wing_psi_c1);
	for(int i=0;i<9;i++) cc_wing_psi_c2[i]=cc_wing_psi_c1[i];
}


void FLAPPATTERN :: Update_Flap_Kinematics(double mtime,int it) {
	double x1,x2;
	
	// x1=mtime;x2=mtime+Check_Frequency_Angle*Dt;
	x1=mtime;x2=mtime+1.0;
	// wingplane phi
	cc_wingplane_phi[0]=Angle_Current_Phi;
	cc_wingplane_phi[1]=0;
	cc_wingplane_phi[2]=0;
	cc_wingplane_phi[3]=Angle_Next_Phi;
	cc_wingplane_phi[4]=0;
	cc_wingplane_phi[5]=0;
	Cal_Poly_Connection(x1, 3, x2, 3, 0, 0, cc_wingplane_phi);
	//c5=right[5];c4=right[4];c3=right[3];c2=right[2];c1=right[1];c0=right[0];
	// wingplane theta
	cc_wingplane_theta[0]=Angle_Current_Theta;
	cc_wingplane_theta[1]=0;
	cc_wingplane_theta[2]=0;
	cc_wingplane_theta[3]=Angle_Next_Theta;
	cc_wingplane_theta[4]=0;
	cc_wingplane_theta[5]=0;
	Cal_Poly_Connection(x1, 3, x2, 3, 0, 0, cc_wingplane_theta);
	//a5=right[5];a4=right[4];a3=right[3];a2=right[2];a1=right[1];a0=right[0];
	// wingplane psi
	cc_wingplane_psi[0]=Angle_Current_Psi;
	cc_wingplane_psi[1]=0;
	cc_wingplane_psi[2]=0;
	cc_wingplane_psi[3]=Angle_Next_Psi;
	cc_wingplane_psi[4]=0;
	cc_wingplane_psi[5]=0;
	Cal_Poly_Connection(x1, 3, x2, 3, 0, 0, cc_wingplane_psi);
	//b5=right[5];b4=right[4];b3=right[3];b2=right[2];b1=right[1];b0=right[0];
	
	
	// Stroke amplitude control
	x1=0.0;x2=0.5;
	cc_wing_phi[0]=Angle_Current_Stroke_Delta; 
	cc_wing_phi[1]=0; 
	cc_wing_phi[2]=0; 
	cc_wing_phi[3]=Angle_Next_Stroke_Delta;
	cc_wing_phi[4]=0;	
	cc_wing_phi[5]=0;
	Cal_Poly_Connection(x1, 3, x2, 3, 0, 0, cc_wing_phi);
	//d5=right[5];d4=right[4];d3=right[3];d2=right[2];d1=right[1];d0=right[0];
	
	
	// AoA control
	// 		start section 
	//		using current psi_c, psi_c will change in next section
	x1=0.0;x2=0.1;
	// 			wing1
	cc_wing_psi_a1[0] = 0; for(int i=0;i<9;i++) cc_wing_psi_a1[0] += cc_wing_psi_c1[i];  // or: right[0] = PI/2.0;
	cc_wing_psi_a1[1] = 0; for(int i=1;i<9;i++) cc_wing_psi_a1[1] += i*cc_wing_psi_c1[i];
	cc_wing_psi_a1[2] = 0; for(int i=2;i<9;i++) cc_wing_psi_a1[2] += i*(i-1)*cc_wing_psi_c1[i];
	cc_wing_psi_a1[3] = 0; for(int i=3;i<9;i++) cc_wing_psi_a1[3] += i*(i-1)*(i-2)*cc_wing_psi_c1[i];
	cc_wing_psi_a1[4] = (mid_aoa + Angle_Next_AoA_Delta)*sin(2*PI*x2) + PI/2.0*(1-sin(2*PI*x2));
	cc_wing_psi_a1[5] = 2.0*PI*(mid_aoa + Angle_Next_AoA_Delta)*cos(2*PI*x2) - PI*PI*cos(2*PI*x2);
	cc_wing_psi_a1[6] = -4.0*PI*PI*(mid_aoa + Angle_Next_AoA_Delta)*sin(2*PI*x2) + 2*PI*PI*PI*sin(2*PI*x2);
	cc_wing_psi_a1[7] = -8.0*PI*PI*PI*(mid_aoa + Angle_Next_AoA_Delta)*cos(2*PI*x2) + 4*PI*PI*PI*PI*cos(2*PI*x2);
	Cal_Poly_Connection(x1, 4, x2, 4, 0, 0, cc_wing_psi_a1);
	// 			wing2
	cc_wing_psi_a2[0] = 0; for(int i=0;i<9;i++) cc_wing_psi_a2[0] += cc_wing_psi_c2[i];  // or: right[0] = PI/2.0;
	cc_wing_psi_a2[1] = 0; for(int i=1;i<9;i++) cc_wing_psi_a2[1] += i*cc_wing_psi_c2[i];
	cc_wing_psi_a2[2] = 0; for(int i=2;i<9;i++) cc_wing_psi_a2[2] += i*(i-1)*cc_wing_psi_c2[i];
	cc_wing_psi_a2[3] = 0; for(int i=3;i<9;i++) cc_wing_psi_a2[3] += i*(i-1)*(i-2)*cc_wing_psi_c2[i];
	cc_wing_psi_a2[4] = (mid_aoa - Angle_Next_AoA_Delta)*sin(2*PI*x2) + PI/2.0*(1-sin(2*PI*x2));
	cc_wing_psi_a2[5] = 2.0*PI*(mid_aoa - Angle_Next_AoA_Delta)*cos(2*PI*x2) - PI*PI*cos(2*PI*x2);
	cc_wing_psi_a2[6] = -4.0*PI*PI*(mid_aoa - Angle_Next_AoA_Delta)*sin(2*PI*x2) + 2*PI*PI*PI*sin(2*PI*x2);
	cc_wing_psi_a2[7] = -8.0*PI*PI*PI*(mid_aoa - Angle_Next_AoA_Delta)*cos(2*PI*x2) + 4*PI*PI*PI*PI*cos(2*PI*x2);
	Cal_Poly_Connection(x1, 4, x2, 4, 0, 0, cc_wing_psi_a2);
	
	//		middle section (downstroke to upstroke)
	x1=0.4;x2=0.6;
	//			wing1
	cc_wing_psi_b1[0] = (mid_aoa + Angle_Next_AoA_Delta)*sin(2*PI*x1) + PI/2.0*(1-sin(2*PI*x1));
	cc_wing_psi_b1[1] = 2.0*PI*(mid_aoa + Angle_Next_AoA_Delta)*cos(2*PI*x1) - PI*PI*cos(2*PI*x1);
	cc_wing_psi_b1[2] = -4.0*PI*PI*(mid_aoa + Angle_Next_AoA_Delta)*sin(2*PI*x1) + 2*PI*PI*PI*sin(2*PI*x1);
	cc_wing_psi_b1[3] = -8.0*PI*PI*PI*(mid_aoa + Angle_Next_AoA_Delta)*cos(2*PI*x1) + 4*PI*PI*PI*PI*cos(2*PI*x1);
	cc_wing_psi_b1[4] = (mid_aoa - Angle_Next_AoA_Delta)*sin(2*PI*x2) + PI/2.0*(1-sin(2*PI*x2));
	cc_wing_psi_b1[5] = 2.0*PI*(mid_aoa - Angle_Next_AoA_Delta)*cos(2*PI*x2) - PI*PI*cos(2*PI*x2);
	cc_wing_psi_b1[6] = -4.0*PI*PI*(mid_aoa - Angle_Next_AoA_Delta)*sin(2*PI*x2) + 2*PI*PI*PI*sin(2*PI*x2);
	cc_wing_psi_b1[7] = -8.0*PI*PI*PI*(mid_aoa - Angle_Next_AoA_Delta)*cos(2*PI*x2) + 4*PI*PI*PI*PI*cos(2*PI*x2);
	cc_wing_psi_b1[8] = PI/2.0; // if 1:1 stroke is needed, let curve pass (t=0.5, AoA=PI/2)
	Cal_Poly_Connection(x1, 4, x2, 4, 0.5, 1, cc_wing_psi_b1);
	//			wing2
	cc_wing_psi_b2[0] = (mid_aoa - Angle_Next_AoA_Delta)*sin(2*PI*x1) + PI/2.0*(1-sin(2*PI*x1));
	cc_wing_psi_b2[1] = 2.0*PI*(mid_aoa - Angle_Next_AoA_Delta)*cos(2*PI*x1) - PI*PI*cos(2*PI*x1);
	cc_wing_psi_b2[2] = -4.0*PI*PI*(mid_aoa - Angle_Next_AoA_Delta)*sin(2*PI*x1) + 2*PI*PI*PI*sin(2*PI*x1);
	cc_wing_psi_b2[3] = -8.0*PI*PI*PI*(mid_aoa - Angle_Next_AoA_Delta)*cos(2*PI*x1) + 4*PI*PI*PI*PI*cos(2*PI*x1);
	cc_wing_psi_b2[4] = (mid_aoa + Angle_Next_AoA_Delta)*sin(2*PI*x2) + PI/2.0*(1-sin(2*PI*x2));
	cc_wing_psi_b2[5] = 2.0*PI*(mid_aoa + Angle_Next_AoA_Delta)*cos(2*PI*x2) - PI*PI*cos(2*PI*x2);
	cc_wing_psi_b2[6] = -4.0*PI*PI*(mid_aoa + Angle_Next_AoA_Delta)*sin(2*PI*x2) + 2*PI*PI*PI*sin(2*PI*x2);
	cc_wing_psi_b2[7] = -8.0*PI*PI*PI*(mid_aoa + Angle_Next_AoA_Delta)*cos(2*PI*x2) + 4*PI*PI*PI*PI*cos(2*PI*x2);
	cc_wing_psi_b2[8] = PI/2.0; // if 1:1 stroke is needed, let curve pass (t=0.5, AoA=PI/2)
	Cal_Poly_Connection(x1, 4, x2, 4, 0.5, 1, cc_wing_psi_b2);
	
	//		end section (upstroke to downstroke)
	x1=0.9;x2=1.1;
	//			wing1
	cc_wing_psi_c1[0] = (mid_aoa - Angle_Next_AoA_Delta)*sin(2*PI*x1) + PI/2.0*(1-sin(2*PI*x1));
	cc_wing_psi_c1[1] = 2.0*PI*(mid_aoa - Angle_Next_AoA_Delta)*cos(2*PI*x1) - PI*PI*cos(2*PI*x1);
	cc_wing_psi_c1[2] = -4.0*PI*PI*(mid_aoa - Angle_Next_AoA_Delta)*sin(2*PI*x1) + 2*PI*PI*PI*sin(2*PI*x1);
	cc_wing_psi_c1[3] = -8.0*PI*PI*PI*(mid_aoa - Angle_Next_AoA_Delta)*cos(2*PI*x1) + 4*PI*PI*PI*PI*cos(2*PI*x1);
	cc_wing_psi_c1[4] = (mid_aoa + Angle_Next_AoA_Delta)*sin(2*PI*x2) + PI/2.0*(1-sin(2*PI*x2));
	cc_wing_psi_c1[5] = 2.0*PI*(mid_aoa + Angle_Next_AoA_Delta)*cos(2*PI*x2) - PI*PI*cos(2*PI*x2);
	cc_wing_psi_c1[6] = -4.0*PI*PI*(mid_aoa + Angle_Next_AoA_Delta)*sin(2*PI*x2) + 2*PI*PI*PI*sin(2*PI*x2);
	cc_wing_psi_c1[7] = -8.0*PI*PI*PI*(mid_aoa + Angle_Next_AoA_Delta)*cos(2*PI*x2) + 4*PI*PI*PI*PI*cos(2*PI*x2);
	cc_wing_psi_c1[8] = PI/2.0; // let curve pass (t=1.0, AoA=PI/2)
	Cal_Poly_Connection(x1, 4, x2, 4, 0.5, 1, cc_wing_psi_c1);
	//			wing2
	cc_wing_psi_c2[0] = (mid_aoa + Angle_Next_AoA_Delta)*sin(2*PI*x1) + PI/2.0*(1-sin(2*PI*x1));
	cc_wing_psi_c2[1] = 2.0*PI*(mid_aoa + Angle_Next_AoA_Delta)*cos(2*PI*x1) - PI*PI*cos(2*PI*x1);
	cc_wing_psi_c2[2] = -4.0*PI*PI*(mid_aoa + Angle_Next_AoA_Delta)*sin(2*PI*x1) + 2*PI*PI*PI*sin(2*PI*x1);
	cc_wing_psi_c2[3] = -8.0*PI*PI*PI*(mid_aoa + Angle_Next_AoA_Delta)*cos(2*PI*x1) + 4*PI*PI*PI*PI*cos(2*PI*x1);
	cc_wing_psi_c2[4] = (mid_aoa - Angle_Next_AoA_Delta)*sin(2*PI*x2) + PI/2.0*(1-sin(2*PI*x2));
	cc_wing_psi_c2[5] = 2.0*PI*(mid_aoa - Angle_Next_AoA_Delta)*cos(2*PI*x2) - PI*PI*cos(2*PI*x2);
	cc_wing_psi_c2[6] = -4.0*PI*PI*(mid_aoa - Angle_Next_AoA_Delta)*sin(2*PI*x2) + 2*PI*PI*PI*sin(2*PI*x2);
	cc_wing_psi_c2[7] = -8.0*PI*PI*PI*(mid_aoa - Angle_Next_AoA_Delta)*cos(2*PI*x2) + 4*PI*PI*PI*PI*cos(2*PI*x2);
	cc_wing_psi_c2[8] = PI/2.0; // let curve pass (t=1.0, AoA=PI/2)
	Cal_Poly_Connection(x1, 4, x2, 4, 0.5, 1, cc_wing_psi_c2);
}


void FLAPPATTERN :: Cal_Poly_Connection(double x1, int order1, double x2, int order2, double x_ctrl, int ctrl_number, double *right) {
	/*
	Using a polynomial to connect motion functions
		Only one middle control point can be handled, and derivatives of it are not controlled
		The order of continuous derivatives at both end need to be the same
	*/
	int i,j,k, order;
	double **temp;
	
	order=order1+order2+ctrl_number;
	temp=new double* [order]; for(i=0;i<order;i++) temp[i]=new double [order];

	// x1
	for(i=0;i<order1;i++) {
		for(j=0;j<order;j++) {	
			if(j<i) {
				temp[i][j]=0;
			}
			else {
				temp[i][j]=pow(x1,j-i);
				for(k=j;k>j-i;k--) {
					temp[i][j] *= k;
				}
			}
		}			
	}
	// x2
	for(i=0;i<order2;i++) {
		for(j=0;j<order;j++) {	
			if(j<i) {
				temp[i+order1][j]=0;
			}
			else {
				temp[i+order1][j]=pow(x2,j-i);
				for(k=j;k>j-i;k--) {
					temp[i+order1][j] *= k;
				}
			}
		}			
	}
	// control point
	if(ctrl_number>0) {
		for(i=order-ctrl_number;i<order;i++)
			for(j=0;j<order;j++)
				temp[i][j]=pow(x_ctrl,j);
	}

	gauss(temp,right,order); // using 'right' to store the results
	// c5=right[5];c4=right[4];c3=right[3];c2=right[2];c1=right[1];c0=right[0];
	for(i=0;i<order;i++) delete[] temp[i];
	delete[] temp;
}


