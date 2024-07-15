#include "utils.h"

void Robotarm::scheme()
{
    // get_torque();

    // float kp = 0.1;
    // float kd = 0.1;

    // for (int i=0; i<actuator_count; i++)
    // {
    //     // cur_d[i] = cur_p[i] + kp*(pos_d[i] - pos_p[i]) - kd*vel_p[i];
    //     std::cout << pos_p[i] << "  " << cur_p[i] << std::endl;
    // } std::cout << std::endl;
    
}

int main(int argc, char **argv)
{
    Robotarm kinova;
    kinova.connet();
    
    // kinova.move_to_home_position();
    
    auto isOk = kinova.current_control();
    if (!isOk)
    {
        std::cout << "There has been an unexpected error in example_cyclic_armbase() function." << std::endl;
    }
    // while (1)
    // {   
    //     kinova.scheme();
    //     usleep(100000);
        
    // }

    kinova.disconnet();
}