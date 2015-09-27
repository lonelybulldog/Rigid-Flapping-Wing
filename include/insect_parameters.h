#ifndef INSECT_PARA_H
#define INSECT_PARA_H

class INSECT_MODULE
{
    public:
	double Re,Frequency,AirDensity,AirViscosity;
	double Dimensional_Mass,               Nondimensional_Mass;
	double Dimensional_WingLength,         Nondimensional_WingLength;
	double Dimensional_MeanWingChord,      Nondimensional_MeanWingChord;
	double Dimensional_WingArea,           Nondimensional_WingArea;
	double Dimensional_Gravity,            Nondimensional_Gravity;
	double Nondimensional_MOI[3][3];
	double Nondimensional_CoM[3];
	
	void Insect_Parameters_Initialize();
	void Insect_Parameters_Nondimensionalize();
	void Insect_Parameters_Update();
};

void INSECT_MODULE :: Insect_Parameters_Initialize()
{
    // fruit fly data
    Dimensional_Mass           =  0.96*1e-6;
    Dimensional_WingLength     =  2.39*1e-3;
    Dimensional_MeanWingChord  =  0.874*1e-3;
    Dimensional_WingArea       =  2.088857849*1e-6;
    Dimensional_Gravity        =  9.81;
    
    //StrokeAmplitude            =  2.0*PI/3.0;
    Frequency                  =  260;
    
    /* temperature : 15 centigrade degrees  */
    AirDensity                 =  1.225;
    AirViscosity               =  1.78*1e-5;
}

void INSECT_MODULE :: Insect_Parameters_Nondimensionalize()
{
    Re=pow(Dimensional_WingLength,2)*Frequency*AirDensity/AirViscosity;
    
    Nondimensional_Mass=Dimensional_Mass/(AirDensity*pow(Dimensional_WingLength,3));
    Nondimensional_WingLength=1.0;
    Nondimensional_MeanWingChord=Dimensional_MeanWingChord/Dimensional_WingLength;
    Nondimensional_WingArea=Dimensional_WingArea/(pow(Dimensional_WingLength,2));
    Nondimensional_Gravity=Dimensional_Gravity/(Dimensional_WingLength*pow(Frequency,2));

    Nondimensional_MOI[0][0]=8.275507637166358; Nondimensional_MOI[0][1]=0; Nondimensional_MOI[0][2]=0;
    Nondimensional_MOI[1][0]=0; Nondimensional_MOI[1][1]=6.456649380801744; Nondimensional_MOI[1][2]=-2.118453249274909;
    Nondimensional_MOI[2][0]=0; Nondimensional_MOI[2][1]=-2.118453249274909;Nondimensional_MOI[2][2]=2.435802992578161;

    Nondimensional_CoM[0]=0; Nondimensional_CoM[1]=-2.5202073690360940e-02; Nondimensional_CoM[2]=-2.5384708611080431e-01;
}

void INSECT_MODULE :: Insect_Parameters_Update()
{
    Re=pow(Dimensional_WingLength,2)*Frequency*AirDensity/AirViscosity;
    Nondimensional_Gravity=Dimensional_Gravity/(Dimensional_WingLength*pow(Frequency,2));
}

#endif /* INSECT_PARA_H */
