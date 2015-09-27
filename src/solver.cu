#include <iostream>

#include "../include/solver.h"

SOLVER :: SOLVER()
{
    std::cout << "Solver is being initialized......";
    
    Time = 0;
    
    for ( int i = 0; i < 3; i++ )
    {
	position[i]             = 0;
	velocity[i]             = 0;
	acceleration[i]         = 0;
	angle[i]                = 0;
	angular_velocity[i]     = 0;
	angular_acceleration[i] = 0;
	
	force[i]                = 0;
	torque[i]               = 0;
	momentum[i]             = 0;
	angular_momentum[i]     = 0;
	
	for ( int j = 0; j < 3; j++ )
	{
	    orientation[i][j]   = 0;
	    inertia[i][j]       = 0;
	}
    }
    
    std::cout << "Initialization done." << std::endl;
}

SOLVER :: ~SOLVER()
{
    std::cout << "Solver is being deleted." << std::endl;
}