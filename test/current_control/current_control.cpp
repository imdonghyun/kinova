#include "utils.h"




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