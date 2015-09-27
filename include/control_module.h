#ifndef CONTROL_MODULE_H
#define CONTROL_MODULE_H

#include <iostream>
#include <vector>
#include <math.h>

class PID_CONTROL_MODULE
{
    public:
	const char*           name;
	double                set;
	double                Kp,Ki,Kd,epsilon,Init;
	double                Max,Min,Max_Change;
	std::vector<double>   record;

	PID_CONTROL_MODULE() { Max_Change=-1; }

	double Activate_controller(double input) 
	{
	    double output;
	
	    std::cout<<"----------------- "<<name<<" is activated! -----------------"<<std::endl;
	 
	    process=input;
	
	    if(record.size()<1)
		pre_error=set-process;
	    else
		pre_error=error;
	
	    error=set-process;
	    Proportional=error;
	    if(fabs(error)<epsilon) {
		Integral=Integral+error;
	    }
	    Derivative=error-pre_error;
	
	    output = Kp*Proportional+Ki*Integral+Kd*Derivative+Init;
	
	    output = (output>Max) ? Max:output;
	    output = (output<Min) ? Min:output;
	
	    if(Max_Change>0) {
		if(record.size()<1) {
		    output = (output>(Init + Max_Change)) ? (Init + Max_Change):output;
		    output = (output<(Init - Max_Change)) ? (Init - Max_Change):output;
		}
		else {
		    int tmp_idx = record.size()-1;
		    output = (output>(record[tmp_idx]+Max_Change)) ? (record[tmp_idx] + Max_Change):output;
		    output = (output<(record[tmp_idx]-Max_Change)) ? (record[tmp_idx] - Max_Change):output;
		}
	    }
	
	    record.push_back(output);
	    return output;
	}

	
	void Record_controller(double protected_var[7])
	{
	    protected_var[0] = Proportional;
	    protected_var[1] = Integral;
	    protected_var[2] = Derivative;
	
	    protected_var[3] = error;
	    protected_var[4] = pre_error;
	    protected_var[5] = process;
	
	    if(record.size()<1)
		protected_var[6] = Init;
	    else
		protected_var[6] = record[record.size()-1];
	}
	
	void Evoke_controller(double protected_var[7]) 
	{
	    Proportional        = protected_var[0];
	    Integral            = protected_var[1];
	    Derivative          = protected_var[2];
	
	    error               = protected_var[3];
	    pre_error           = protected_var[4];
	    process             = protected_var[5];
	
	    record.push_back(protected_var[6]);
	}
	
    protected:       
	double           Proportional,Integral,Derivative;
	double           error,pre_error,process;
};

#endif /* CONTROL_MODULE_H */
