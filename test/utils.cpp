#include "utils.h"
#include "Robotics.h"

std::function<void(k_api::Base::ActionNotification)>
Robotarm::check_for_end_or_abort(bool& finished)
{
    return [&finished](k_api::Base::ActionNotification notification)
    {
        std::cout << "EVENT : " << k_api::Base::ActionEvent_Name(notification.action_event()) << std::endl;

        // The action is finished when we receive a END or ABORT event
        switch(notification.action_event())
        {
        case k_api::Base::ActionEvent::ACTION_ABORT:
        case k_api::Base::ActionEvent::ACTION_END:
            finished = true;
            break;
        default:
            break;
        }
    };
}

int64_t Robotarm::GetTickUs()
{
#if defined(_MSC_VER)
    LARGE_INTEGER start, frequency;

    QueryPerformanceFrequency(&frequency);
    QueryPerformanceCounter(&start);

    return (start.QuadPart * 1000000)/frequency.QuadPart;
#else
    struct timespec start;
    clock_gettime(CLOCK_MONOTONIC, &start);
    
    return (start.tv_sec * 1000000LLU) + (start.tv_nsec / 1000);
#endif
}

void Robotarm::move_to_home_position()
{
    // Make sure the arm is in Single Level Servoing before executing an Action
    auto servoingMode = k_api::Base::ServoingModeInformation();
    servoingMode.set_servoing_mode(k_api::Base::ServoingMode::SINGLE_LEVEL_SERVOING);
    base->SetServoingMode(servoingMode);
    std::this_thread::sleep_for(std::chrono::milliseconds(500));

    // Move arm to ready position
    std::cout << "Moving the arm to a safe position" << std::endl;
    auto action_type = k_api::Base::RequestedActionType();
    action_type.set_action_type(k_api::Base::REACH_JOINT_ANGLES);
    auto action_list = base->ReadAllActions(action_type);
    auto action_handle = k_api::Base::ActionHandle();
    action_handle.set_identifier(0);
    for (auto action : action_list.action_list()) 
    {
        if (action.name() == "Home") 
        {
            action_handle = action.handle();
        }
    }

    if (action_handle.identifier() == 0) 
    {
        std::cout << "Can't reach safe position, exiting" << std::endl;
    } 
    else 
    {
        bool action_finished = false; 
        // Notify of any action topic event
        auto options = k_api::Common::NotificationOptions();
        auto notification_handle = base->OnNotificationActionTopic(
            check_for_end_or_abort(action_finished),
            options
        );

        base->ExecuteActionFromReference(action_handle);

        while(!action_finished)
        { 
            std::this_thread::sleep_for(ACTION_WAITING_TIME);
        }

        base->Unsubscribe(notification_handle);
    }
}

void Robotarm::set_position_mode()
{
    auto control_mode_message = k_api::ActuatorConfig::ControlModeInformation();
    control_mode_message.set_control_mode(k_api::ActuatorConfig::ControlMode::POSITION);

    //1~5: 1번 ~ 5번 조인트, 7: 6번 조인트
    actuator_config->SetControlMode(control_mode_message, 1);
    actuator_config->SetControlMode(control_mode_message, 2);
    actuator_config->SetControlMode(control_mode_message, 3);
    actuator_config->SetControlMode(control_mode_message, 4);
    actuator_config->SetControlMode(control_mode_message, 5);
    actuator_config->SetControlMode(control_mode_message, 7);

    for(int i = 0; i < actuator_count; i++)
    {
        pos_d[i] = base_feedback.actuators(i).position();
        base_command.mutable_actuators(i)->set_position(pos_d[i]);
    }
}

void Robotarm::set_current_mode()
{
    servoingMode.set_servoing_mode(k_api::Base::ServoingMode::LOW_LEVEL_SERVOING);
    base->SetServoingMode(servoingMode);

    // Initialize each actuator to its current position
    for(int i = 0; i < actuator_count; i++)
    {
        pos_d[i] = base_feedback.actuators(i).position();
        base_command.add_actuators()->set_position(pos_d[i]);
    }

    auto control_mode_message = k_api::ActuatorConfig::ControlModeInformation();
    control_mode_message.set_control_mode(k_api::ActuatorConfig::ControlMode::CURRENT);

    //1~5: 1번 ~ 5번 조인트, 7: 6번 조인트
    actuator_config->SetControlMode(control_mode_message, 1);
    actuator_config->SetControlMode(control_mode_message, 2);
    actuator_config->SetControlMode(control_mode_message, 3);
    actuator_config->SetControlMode(control_mode_message, 4);
    actuator_config->SetControlMode(control_mode_message, 5);
    actuator_config->SetControlMode(control_mode_message, 7);
}

void Robotarm::scheme()
{
    VectorXf g(6);
    VectorXf q(6);
    VectorXf qdot(6);

    qdot.setZero();
    
    
    // std::cout << pos_p << std::endl;
    for (int i=0; i<6; i++)
    {
        q(i) = pos_p[i]*M_1_PI/180;
    }
    // std::cout << q << std::endl;
    updateFKList(q, qdot);
    g = systemGravity();

    std::cout << g << std::endl << std::endl;

    for (int i=0; i<6; i++)
    {
        cur_d[i] = g(i) / Kt[i];
    }
}

bool Robotarm::current_control()
{
    bool return_status = true;
    int timer_count = 0;
    int64_t now = 0;
    int64_t last = 0;

    int timeout = 0;

    std::cout << "Initializing the arm for low-level control" << std::endl;
    try
    {
        set_current_mode();
        setConstant();

        // Define the callback function used in Refresh_callback
        // auto lambda_fct_callback = [](const Kinova::Api::Error &err, const k_api::BaseCyclic::Feedback data)
        // {
        //     // We are printing the data of the moving actuator just for the example purpose,
        //     // avoid this in a real-time loop
        //     std::string serialized_data;
        //     google::protobuf::util::MessageToJsonString(data.actuators(0), &serialized_data);
        //     std::cout << serialized_data << std::endl << std::endl;
        // };

        // Real-time loop
        
        
        while(timer_count < (time_duration * 1000))
        // while(1)
        {
            now = GetTickUs();
            if(now - last > 1000)
            {
                get_position();
                scheme();
                // 0~5: 1번 ~ 6번 joint
                for (int i=0; i<actuator_count; i++) 
                {
                    // std::cout << cur_d[i] << " ";
                    base_command.mutable_actuators(i)->set_position(pos_p[i]);
                    base_command.mutable_actuators(i)->set_current_motor(0);

                    cur_p[i] = base_feedback.actuators(i).current_motor();
                    std::cout << cur_p[i] << std::endl;

                }
                std::cout << std::endl;
                try
                {
                    // base_cyclic->Refresh_callback(base_command, lambda_fct_callback, 0);
                    base_feedback = base_cyclic->Refresh(base_command, 0);
                }
                catch(...)
                {
                    timeout++;
                }
                
                timer_count++;
                last = GetTickUs();
            }
        }
    }
    catch (k_api::KDetailedException& ex)
    {
        std::cout << "Kortex error: " << ex.what() << std::endl;
        return_status = false;
    }
    catch (std::runtime_error& ex2)
    {
        std::cout << "Runtime error: " << ex2.what() << std::endl;
        return_status = false;
    }
    
    set_position_mode();
    // for (int i=0; i<3000000; i++);
    // for(int i = 0; i < actuator_count; i++)
    // {
    //     cur_p[i] = base_feedback.actuators(i).current_motor();
    //     std::cout << cur_p[i] << " ";
    // }
    // Set back the servoing mode to Single Level Servoing
    servoingMode.set_servoing_mode(k_api::Base::ServoingMode::SINGLE_LEVEL_SERVOING);
    base->SetServoingMode(servoingMode);

    // Wait for a bit
    std::this_thread::sleep_for(std::chrono::milliseconds(2000));

    return return_status;
}

void Robotarm::get_position()
{
    for (int i=0; i<actuator_count; i++){
        pos_p[i]=base_feedback.actuators(i).position();
    }
}

void Robotarm::get_velocity()
{
    for (int i=0; i<actuator_count; i++){
        vel_p[i]=base_feedback.actuators(i).velocity();
    }
}

void Robotarm::get_current()
{
    for (int i=0; i<actuator_count; i++){
        cur_p[i]=base_feedback.actuators(i).current_motor();
    }
}

void Robotarm::get_torque()
{
    for (int i=0; i<actuator_count; i++){
        tau_p[i]=base_feedback.actuators(i).torque();
    }
}

void Robotarm::connet() //init members and connect kinova
{
    ip_address = "192.168.1.10";
    username = "admin";
    password = "admin";

    transport = new k_api::TransportClientTcp();
    router = new k_api::RouterClient(transport, [](k_api::KError err){ cout << "_________ callback error _________" << err.toString(); });

    transport_real_time = new k_api::TransportClientUdp();
    router_real_time = new k_api::RouterClient(transport_real_time, [](k_api::KError err){ cout << "_________ callback error _________" << err.toString(); });

    session_manager = new k_api::SessionManager(router);
    session_manager_real_time = new k_api::SessionManager(router_real_time);

    // Create services
    base = new k_api::Base::BaseClient(router);
    base_cyclic = new k_api::BaseCyclic::BaseCyclicClient(router_real_time);
    actuator_config = new k_api::ActuatorConfig::ActuatorConfigClient(router);

    velocity = 20.0f;
    time_duration = DURATION;

    // Create API objects
    std::cout << "Creating transport objects" << std::endl;
    transport->connect(ip_address, PORT);

    std::cout << "Creating transport real time objects" << std::endl;
    transport_real_time->connect(ip_address, PORT_REAL_TIME);

    // Set session data connection information
    create_session_info.set_username(username);
    create_session_info.set_password(password);
    create_session_info.set_session_inactivity_timeout(60000);   // (milliseconds)
    create_session_info.set_connection_inactivity_timeout(2000); // (milliseconds)

    // Session manager service wrapper
    std::cout << "Creating sessions for communication" << std::endl;
    session_manager->CreateSession(create_session_info);
    session_manager_real_time->CreateSession(create_session_info);
    std::cout << "Sessions created" << std::endl;

    base_feedback = base_cyclic->RefreshFeedback();

    actuator_count = base->GetActuatorCount().count();
    base->ClearFaults();

    // for (int i=0; i<actuator_count; i++) 
    // {
    //     base_command.mutable_actuators(i)->clear_flags();
    // }
    // base_cyclic->Refresh(base_command, 0);
    
    cnt=0;
}

void Robotarm::disconnet()
{
    // Close API session
    session_manager->CloseSession();
    session_manager_real_time->CloseSession();

    // Deactivate the router and cleanly disconnect from the transport object
    router->SetActivationStatus(false);
    transport->disconnect();
    router_real_time->SetActivationStatus(false);
    transport_real_time->disconnect();

    // Destroy the API
    delete base;
    delete base_cyclic;
    delete actuator_config;
    delete session_manager;
    delete session_manager_real_time;
    delete router;
    delete router_real_time;
    delete transport;
    delete transport_real_time;
}
