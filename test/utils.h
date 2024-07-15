#include <iostream>
#include <string>
#include <vector>
#include <math.h>

#include <KDetailedException.h>

#include <BaseClientRpc.h>
#include <BaseCyclicClientRpc.h>
#include <ActuatorConfigClientRpc.h>
#include <SessionClientRpc.h>
#include <SessionManager.h>

#include <RouterClient.h>
#include <TransportClientTcp.h>
#include <TransportClientUdp.h>

#include <google/protobuf/util/json_util.h>

#if defined(_MSC_VER)
#include <Windows.h>
#else
#include <unistd.h>
#endif
#include <time.h>
#include "config.h"

namespace k_api = Kinova::Api;

#define PORT 10000
#define PORT_REAL_TIME 10001

#define DURATION 1             // Network timeout (seconds)

// Waiting time during actions
const auto ACTION_WAITING_TIME = std::chrono::seconds(1);

// Create closure to set finished to true after an END or an ABORT

class Robotarm
{
    private:
    std::string ip_address;
    std::string username;
    std::string password;

    k_api::TransportClientTcp* transport;
    k_api::RouterClient* router;

    k_api::TransportClientUdp* transport_real_time;
    k_api::RouterClient* router_real_time;
    
    k_api::Session::CreateSessionInfo create_session_info;
    
    k_api::SessionManager* session_manager;
    k_api::SessionManager* session_manager_real_time;

    // Create services
    k_api::Base::BaseClient* base;
    k_api::BaseCyclic::BaseCyclicClient* base_cyclic;
    k_api::ActuatorConfig::ActuatorConfigClient* actuator_config;

    k_api::BaseCyclic::Feedback base_feedback;
    k_api::BaseCyclic::Command  base_command;

    k_api::Base::ServoingModeInformation servoingMode;

    float velocity;         // Default velocity of the actuator (degrees per seconds)
    float time_duration; // Duration of the example (seconds)

    float pos_p[6];
    float vel_p[6];
    float cur_p[6];
    float tau_p[6];

    float pos_d[6];
    float vel_d[6];
    float cur_d[6];
    float tau_d[6];

    int actuator_count;
    int cnt;
    

    public:

    
    std::function<void(k_api::Base::ActionNotification)> 
    check_for_end_or_abort(bool& finished);

    int64_t GetTickUs();
    void move_to_home_position();

    void set_current_mode();
    bool current_control();

    void get_position();
    void get_velocity();
    void get_current();
    void get_torque();

    void state_update();
    void scheme();
    // void calculation();

    void connet();
    void disconnet();
};

