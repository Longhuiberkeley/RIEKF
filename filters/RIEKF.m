%%
% Author: Long Hui
% modified from github repo: Invariant-ekf

% this assumes ALL IMU, GPS, and vehicle are in the SAME frame

%%
% Use RIEKF to be the main method for IMU prediction and vehicle constraint
% correction

% then use LIEKF for GPS update




classdef RIEKF < handle
    %% Left-Invariant filter class, predicts next state, corrects prediction
    properties
        mu;         %Pose Mean
        bias;       %Bias of gyro and accelerometer = [wb,ab]';
        g;          %gravity 
        Sigma;      %Pose Sigma

        cov_g;      %gyro noise
        cov_a;      %acc noise
        cov_gb;     %gyro bias
        cov_ab;     %acc bias

        V_gps;          %observation noise of position, GPS
        V_odometer_nonholonomic; %observation noise of odometer and noise of vehicle dynamics slippage

        Q;          %all covariance in process model, might be WRONG
    end



    methods

        % init
        function obj = RIEKF(R0, p0, v0, cov_g_, cov_a_, cov_gb_, cov_ab_, V_gps, V_odometer_, g_)
            % Set the initial state and noise profile
            if nargin == 0
                R0 = eye(3);
                p0 = zeros(3,1);
                v0 = zeros(3,1);
                cov_g_ = eye(3) * 20;
                cov_a_ = eye(3) * 20;
                cov_gb_ = eye(3);
                cov_ab_ = eye(3);
                V_gps = [4.6778    1.9437    0.0858;
                      1.9437   11.5621    5.8445;
                      0.0858    5.8445   22.4051]*1000;
                V_odometer_ = diag([1, 1, 1] );
                g_ = [0; 0; -9.81];
            end

            obj.mu = blkdiag(R0, eye(2));
            obj.mu(1:3,4) = v0;
            obj.mu(1:3,5) = p0;
            obj.g = g_; 

            obj.Sigma = eye(15);
            obj.bias = zeros(6,1);

            obj.cov_g = cov_g_; 
            obj.cov_a = cov_a_;
            obj.cov_gb = cov_gb_;
            obj.cov_ab = cov_ab_;
            obj.V_gps = V_gps;
            obj.V_odometer_nonholonomic = V_odometer_;

            obj.Q = blkdiag([
                obj.cov_g,zeros(3),zeros(3),zeros(3),zeros(3);
                zeros(3),obj.cov_a,zeros(3),zeros(3),zeros(3);
                zeros(3),zeros(3),eye(3),zeros(3),zeros(3);
                zeros(3),zeros(3),zeros(3),obj.cov_gb,zeros(3),;
                zeros(3),zeros(3),zeros(3),zeros(3),obj.cov_ab]);
        end

        function prediction(obj, w, a, dt) 

            % Store state in convenient variables
            [R, v, p] = obj.getState();
            wb = obj.bias(1:3);
            ab = obj.bias(4:6);



            R = R * expm( skew( ( w - wb ) * dt ) );
            v = v + R * ( a - ab ) * dt + obj.g * dt;
            p = p + v * dt + 0.5 * R * ( a - ab ) * dt ^ 2 + 0.5 * obj.g * dt ^ 2;
            obj.mu = [
                R, v, p;
                zeros(2,3), eye(2)
                ];

                
            ja = obj.jacobian_A();
            phi = expm( ja * dt );

            jb = obj.jacobian_B();
            cov_w = blkdiag(obj.cov_g, obj.cov_a, zeros(3,3), obj.cov_gb, obj.cov_ab);

            Q_ =  jb* cov_w * jb'; % cpv_w or obj.Q
            Q_disc = phi * Q_ * phi' *dt;
            obj.Sigma = phi * obj.Sigma * phi' + Q_disc;

        end

        %% GPS update methods, as a LIEFK method
        function correction(obj,GPS)
            Y = [GPS;0;1];
            b = [0;0;0;0;1];
            H = [zeros(3),zeros(3), eye(3), zeros(3), zeros(3)];

            [R, ~, ~] = obj.getState();

            N = R' * obj.V_gps * R;

            % P_r = Ad * P_l * Ad'; 
            ad = obj.Adjoint_bias();

            left_sigma = eye(15)/ad * obj.Sigma / ad'; 
            S = H * left_sigma * H' + N;
            K = left_sigma * H' / S;     % Kalman gain
            K_X = K(1:9,:);
            K_B = K(10:15,:);
            PI = [eye(3), zeros(3,2)];
            nu = eye(5)/obj.mu * Y - b;        % Innovation
            delta_X = K_X * PI * nu;
            delta_B = K_B * PI * nu;
            xi = makeTwist(delta_X);     % basically the obj.wedge operator for our 9-dimensional state space

            obj.mu = obj.mu * expm(xi);
            obj.bias = obj.bias + delta_B;
            left_sigma = (eye(15)- K * H) * left_sigma * (eye(15)- K * H)' + K * N * K';

            % map it back to RIEKF
            obj.Sigma = ad * left_sigma * ad'; 

        end


        %% odometry + non-holonomic are methods in the RIGHT-Invariant EKF

        % odometry function corrects our state with odometer and non-holonomic
        % constraint



        function odometry(obj, odo)
%             odo_positive = -odo;
            odo_positive = odo;
            Y_odo_measurement = [odo_positive; 0; 0; 1; 0];
            b_odo = [0;0;0;1;0];


            z = zeros(3);
            [R, ~, ~] = obj.getState();
            H = [z, -eye(3), z, z, z];
            
            N = R * obj.V_odometer_nonholonomic * R' ;

            S = H * obj.Sigma * H' + N;
            K = obj.Sigma * H' / S;     % Kalman gain
            K_X = K(1:9,:);
            K_B = K(10:15,:);


            PI = [eye(3), zeros(3,2)];

            nu = obj.mu * Y_odo_measurement - b_odo;        % Innovation
            delta_X = K_X * PI * nu;
            delta_B = K_B * PI * nu;
            xi = makeTwist(delta_X);     % basically the obj.wedge operator for our 9-dimensional state space

            obj.mu = expm(xi) * obj.mu ;
            obj.bias = obj.bias + delta_B;

            obj.Sigma= (eye(15)- K * H) * obj.Sigma * (eye(15)- K * H)' + K * N * K';
        end


        %% wheel_encoders, assume we have encoder data from both wheels


        %% nonholonomic function corrects our state with non-holonomic
        % constraint only
 
        function nonholonomic(obj)
            

            z = zeros(3);
            [R, v, ~] = obj.getState();

%             reduction_matrix = [eye(2), zeros(2,1)]; % this is WRONG
%             H_vehicle = [z, R', z, -skew(p), z]; % this is WRONG

            reduction_matrix = [zeros(2,1), eye(2)]; %correct
            H_vehicle = [z, R', z, z, z];
            H = reduction_matrix * H_vehicle; 
            
            N_ = obj.V_odometer_nonholonomic(2:3, 2:3);

            S = H * obj.Sigma * H' + N_;
            K = obj.Sigma * H' / S;     % Kalman gain

            pseudo_measurement = [0;0] ;
            vehicle_frame_velocity = R' * v;
            vehicle_slippage = vehicle_frame_velocity(2:3); 
            
            error = pseudo_measurement -vehicle_slippage; 

            delta = K* error ;

            delta_X = delta(1:9);
            delta_B = delta(10:15);
            xi = makeTwist(delta_X);

            obj.mu = expm(xi) * obj.mu ;
            obj.bias = obj.bias + delta_B;

            obj.Sigma = (eye(15)- K * H) * obj.Sigma * (eye(15)- K * H)' + K * N_ * K';

        end

        %% ================================================================
        %% helper functions

        function [R, v, p] = getState(obj)
            R = obj.mu(1:3, 1:3);
            v = obj.mu(1:3, 4);
            p = obj.mu(1:3, 5);
        end


        % Adjoint of SE_2 (3)
        % Ad_X is a 9*9 matrix
        function  Ad_X = Adjoint_noBias(obj)
            [R, v, p] = obj.getState();
            z = zeros(3);
            Ad_X = [R, z, z;
                skew(v)*R, R, z;
                skew(p)*R, z, R];
        end



        % Adjoint of SE_4 (3) in R^5*
        % Ad_X is a 15* 15 matrix
        function  Ad_X = Adjoint_bias(obj)
            [R, v, p] = obj.getState();
            wb = obj.bias(1:3);
            ab = obj.bias(4:6);

            z = zeros(3);
            Ad_X = [R, z, z, z, z;
                skew(v)*R, R, z, z, z;
                skew(p)*R, z, R, z, z;
                skew(wb)*R, z, z, R, z;
                skew(ab)*R, z, z, z, R];
        end


        function y = jacobian_A(obj)
            [R, v, p] = getState(obj); 
            z =zeros(3,3);
            y = [z, z, z, -R, z; 
                skew(obj.g), z, z, -skew(v)*R, -R;
                zeros(3), eye(3), zeros(3), -skew(p)*R, z;
                z, z, z, z, z;
                z, z, z, z, z];
        end

        function y = jacobian_B(obj)
            % this is a 15*15 matrix
            ad = obj.Adjoint_noBias(); 
            y = [ad, zeros(9,6);
                zeros(6,9) eye(6)];
        end

        function m = G(obj)
            % same as jacobian_B
            % this is a 15*15 matrix 
            [R, v, p] = getState(obj);
            z =zeros(3,3);
            m = [R, z, z, z;
                skew(v)*R, R, z, z;
                skew(p)*R, z, z, z; 
                z, z, eye(3), z; 
                z, z, z, eye(3)];

        end

    end
end