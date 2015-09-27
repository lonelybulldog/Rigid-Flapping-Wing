#ifndef GLOBAL_H
#define GLOBAL_H

REF_FRAME        lab,insect,wing_1,wing_2, wingplane_1, wingplane_2;
INSECT_MODULE    insect_parameters;
FLAPPATTERN      my_flap_pattern;


double CoM[3],CoM_Global[3], Wing_Root[3];
double bodyforce[3]={0}, bodytorque[3]={0};
double wing1force[3]={0}, wing1torque[3]={0};
double wing2force[3]={0}, wing2torque[3]={0};

double bodycentre_old[3]={8,8,8};
double bodycentre[3]={8,8,8};
double bodycentre_tmp[3]={8,8,8};

void Transpose(double a[3][3],double b[3][3]);
void Multiply(double a[3][3],double b[3][3],double c[3][3]);

void Write_Result_DAT(int iter, MESH& mesh);

#endif /* GLOBAL_H */
