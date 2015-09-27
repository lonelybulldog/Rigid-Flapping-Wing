#ifndef SOLVER_H
#define SOLVER_H


class SOLVER
{
public:
    double        Time;
    
    double        position[3];
    double        velocity[3];
    double        acceleration[3];
    double        angle[3];
    double        angular_velocity[3];
    double        angular_acceleration[3];
    double        orientation[3][3];
    double        force[3];
    double        torque[3];
    double        momentum[3];
    double        inertia[3][3];
    double        angular_momentum[3];
    
    double        position_old[3];
    double        velocity_old[3];
    double        acceleration_old[3];
    double        angle_old[3];
    double        angular_velocity_old[3];
    double        angular_acceleration_old[3];
    double        orientation_old[3][3];
    double        force_old[3];
    double        torque_old[3];
    double        momentum_old[3];
    double        inertia_old[3][3];
    double        angular_momentum_old[3];
	
    double        position_tmp[3];
    double        velocity_tmp[3];
    double        acceleration_tmp[3];
    double        angle_tmp[3];
    double        angular_velocity_tmp[3];
    double        angular_acceleration_tmp[3];
    double        orientation_tmp[3][3];
    double        force_tmp[3];
    double        torque_tmp[3];
    double        momentum_tmp[3];
    double        inertia_tmp[3][3];
    double        angular_momentum_tmp[3];
    
    
    SOLVER();
    ~SOLVER();
};










#endif /* SOLVER_H */