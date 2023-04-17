%%
% Author: Long Hui
% modified from github repo: Invariant-ekf

% this assumes ALL IMU, GPS, and vehicle are in the SAME frame

%%
% Use LIEKF to be the main method for IMU prediction and GPS correction,
% then use




classdef LIEKF < handle
    %% Left-Invariant filter class, predicts next state, corrects prediction
    properties
        mu;         %Pose Mean
        bias;       %Bias of gyro and accelerometer = [wb,ab]';
        g;          %gravity 
        Sigma;      %Pose Sigma
        A;          %Process model - Not currently being used
        cov_g;      %gyro noise
        cov_a;      %acc noise
        cov_gb;     %gyro bias
        cov_ab;     %acc bias
        V;          %observation noise of position, GPS

        V_odometer_nonholonomic; %observation noise of odometer and noise of vehicle dynamics slippage

        Q;          %all covariance in process model
        sigma_cart
    end



    methods

        % init
        function obj = LIEKF(R0, p0, v0, cov_g_, cov_a_, cov_gb_, cov_ab_, V_, V_odometer_, g_)
            % Set the initial state
            if nargin == 0
                R0 = eye(3);
                p0 = zeros(3,1);
                v0 = zeros(3,1);
                cov_g_ = eye(3) * 20;
                cov_a_ = eye(3) * 20;
                cov_gb_ = eye(3);
                cov_ab_ = eye(3);
                V_ = [4.6778    1.9437    0.0858;
                      1.9437   11.5621    5.8445;
                      0.0858    5.8445   22.4051]*1000;
                V_odometer_ = diag([0.5, 0.5, 0.5]);
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
            obj.V = V_;
            obj.V_odometer_nonholonomic = V_odometer_;

            obj.Q = blkdiag([
                obj.cov_g,zeros(3),zeros(3),zeros(3),zeros(3);
                zeros(3),obj.cov_a,zeros(3),zeros(3),zeros(3);
                zeros(3),zeros(3),eye(3),zeros(3),zeros(3);
                zeros(3),zeros(3),zeros(3),obj.cov_gb,zeros(3),;
                zeros(3),zeros(3),zeros(3),zeros(3),obj.cov_ab]);
            obj.A = @(wt,at) [
                -skew(wt), zeros(3),  zeros(3), -eye(3), zeros(3);
                -skew(at), -skew(wt), zeros(3), zeros(3), -eye(3);
                zeros(3), eye(3), -skew(wt), zeros(3),zeros(3);
                zeros(3), zeros(3), zeros(3), zeros(3), zeros(3);
                zeros(3), zeros(3), zeros(3), zeros(3), zeros(3)
                ];
        end

        function prediction(obj, w, a, dt)  %TBC bias
            %skew = @(u) skew(u);

            % Predicts position from gyro/accelerometer data;
            % LONG: I think we can just eye(3) as an approximation
            if norm(w) > 1e-8
                gamma0 = @(phi) eye(3) + sin(norm(phi))/norm(phi) * skew(phi) + (1-cos(norm(phi)))/(norm(phi)^2) * skew(phi)^2;

                gamma1 = @(phi) eye(3) + (1-cos(norm(phi)))/(norm(phi)^2) * skew(phi) + (norm(phi,2) - sin(norm(phi)))/(norm(phi)^3) * skew(phi)^2;

                gamma2 = @(phi) 0.5*eye(3) + (norm(phi,2) - sin(norm(phi,2)))/(norm(phi,2)^3) * skew(phi)+ (norm(phi,2)^2 + 2*cos(norm(phi,2)) - 2)/(2*(norm(phi,2)^4)) * skew(phi)^2;
            else
                gamma0 = @(phi) eye(3);
                gamma1 = @(phi) eye(3);
                gamma2 = @(phi) 1/2*eye(3);
            end

            % Bias stuff?
            wb = obj.bias(1:3); %TBC
            ab = obj.bias(4:6);
            wt = w - wb;    %true noisy value
            at = a - ab;

            % Store state in convenient variables
            [R, v, p] = obj.getState();

            % Integrate the angular rates
            % This is equivalent to
            %   R_k = R*expm(skew(w*dt))
            % only using the gamma function to
            % construct the rotation matrix expm(skew(w*dt))
            Rk = R*gamma0(w*dt);

            % An accelerometer can't tell the difference between
            % the pull of a gravitational field in one direction
            % and an acceleration in the opposite direction.
            % Let g = [0;0;-9.81], then the measured accel a_m is
            %   a_m = R^T*(a_e - g)
            % where a_e is the actual acceleration in the earth frame
            % and R is the rotation of the body frame with respect to earth.
            % Solving for a_e (to integrate), we get
            %   a_e = R*a_m + g
            %
            % With the earth frame acceleration, integrate once for velocity
            % and a second time for position.  Formulas with the gamma function
            % from slides (I don't know the derivation but they appear to work)
            
            vk = v + R*gamma1(wt*dt)*at*dt + obj.g *dt;
            pk = p + v*dt + R*gamma2(wt*dt)*at*(dt^2) + 0.5* obj.g *(dt^2);
            phi = expm(obj.A(wt,at)*dt);

            % Set the mean to the predicted value
            obj.mu = [
                Rk, vk, pk;
                zeros(2,3), eye(2)
                ];
            %             obj.bias = zeros(6,1);
            obj.Sigma = phi*obj.Sigma*(phi') + phi*obj.Q*(phi')*dt;
        end

        %% GPS update methods

        %GPS 3x1 is this in R^3 ECEF/NED/ENU??
        function correction(obj,GPS)
            Y = [GPS;0;1];
            H = [zeros(3),zeros(3), eye(3), zeros(3), zeros(3)];

            [R, ~, ~] = obj.getState();

            N = R' * obj.V * R;

            S = H * obj.Sigma * H' + N;
            K = obj.Sigma * H' / S;     % Kalman gain
            K_X = K(1:9,:);
            K_B = K(10:15,:);
            PI = [eye(3), zeros(3,2)];
            nu = eye(5)/obj.mu * Y;        % Innovation
            delta_X = K_X * PI * nu;
            delta_B = K_B * PI * nu;
            xi = makeTwist(delta_X);     % basically the obj.wedge operator for our 9-dimensional state space

            obj.mu = obj.mu * expm(xi);
            obj.bias = obj.bias + delta_B;
            obj.Sigma = (eye(15)- K * H) * obj.Sigma * (eye(15)- K * H)' + K * N * K';
        end


        %% odometry and non-holonomic are methods in the RIGHT-Invariant EKF

        % odometry function corrects our state with odometer and non-holonomic
        % constraint

        % nonholonomic function corrects our state with non-holonomic
        % constraint only


        % Procedure:
        % -- Need adjoint map;
        % -- Need xi_r;
        % -- Need covariance_r, Sigma_r
        % -- Need Sigma_r turn back into LIEKF covariance, Simga

        function odometry(obj, odo)


            odo_positive = -odo;

            Y_odo_measurement = [odo_positive; 0; 0; 1; 0];
%             b_odo = [0;0;0;1;0];


            z = zeros(3);
            [R, ~, ~] = obj.getState();
            %             H = [z, R', z, z, z];
            H = [z, -eye(3), z, z, z];
            


            N = R * obj.V_odometer_nonholonomic * R' ;

            Ad_X15 = obj.Adjoint_bias;

            Sigma_right = Ad_X15 * obj.Sigma * Ad_X15';

            S = H * Sigma_right * H' + N;
            K = Sigma_right * H' / S;     % Kalman gain
            K_X = K(1:9,:);
            K_B = K(10:15,:);


            PI = [eye(3), zeros(3,2)];

            nu = obj.mu * Y_odo_measurement;        % Innovation
            delta_X = K_X * PI * nu;
            delta_B = K_B * PI * nu;
            xi = makeTwist(delta_X);     % basically the obj.wedge operator for our 9-dimensional state space

            obj.mu = expm(xi) * obj.mu ;
            obj.bias = obj.bias + delta_B;

            Sigma_right= (eye(15)- K * H) * Sigma_right * (eye(15)- K * H)' + K * N * K';
            obj.Sigma = eye(15)/Ad_X15 * Sigma_right * Ad_X15; % Put it back to LIEKF

        end


        %% wheel_encoders, assume we have encoder data from both wheels

        %% this is odometer with 2 entries method, RIGHT invariant EKF
        % enc = [enc_left; enc_right]; 
        function wheel_encoders(obj, enc)
        end 


        function nonholonomic(obj)

%             Y_nonholo = [0; 0; 1; 0];
%             b_extract_v = [0;0;0;1;0];


            z = zeros(3);
            [R, ~, ~] = obj.getState();
            %             H = [z, R', z, z, z];
            H = [z, -eye(3), z, z, z];
            H = [z, R', z, z, z]; 

            H_ = H(2:3, :); 


            N = R * obj.V_odometer_nonholonomic * R' ;
            N_ = N(2:3,2:3);

            Ad_X15 = obj.Adjoint_bias;

            Sigma_right = Ad_X15 * obj.Sigma * Ad_X15';

            S = H_ * Sigma_right * H_' + N_;
            K = Sigma_right * H_' / S;     % Kalman gain
            K_X = K(1:9,:);
            K_B = K(10:15,:);


            nu = obj.mu(2:3,4); % we might have problem doing this, we might lose log-linearity property

            delta_X = K_X * nu;
            delta_B = K_B * nu;
            xi = makeTwist(delta_X);

            obj.mu = expm(xi) * obj.mu ;
            obj.bias = obj.bias + delta_B;

            Sigma_right= (eye(15)- K * H_) * Sigma_right * (eye(15)- K * H_)' + K * N_ * K';
            obj.Sigma = eye(15)/Ad_X15 * Sigma_right * Ad_X15; % Put it back to LIEKF

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

        function invAd_X = Inverse_Adjoint_noBias(obj)
            invAd_X = inv(Adjoint_noBias(obj));
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

        function invAd_X = Inverse_Adjoint_bias(obj)
            invAd_X = inv(Adjoint_bias(obj));
        end



    end
end