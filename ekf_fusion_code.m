%% main

% Read data
csvFile = 'sensor_data/sensor_data.csv';
data = readtable(csvFile);
numColumns = size(data, 2);
% Data Filter
t = transpose(data{(2:end), 1});
for  i = 2 : 10
    meas_f = 100; % sample freq
    Fc = 1; % cutoff frequency
    [Blp, Alp] = cheby2(4, 100, Fc / meas_f); 
    Xi = filtfilt(Blp, Alp, data{(2: end), i});
    if i < 5
        gyro(i - 1, :) = Xi;
    end
    if i > 4 && i < 8
        acc(i - 4, :) = Xi;
    end
    if i > 7
        mag(i - 7, :) = Xi;
    end  
end
gyro = deg2rad(gyro);
% Gyro to angle
angle = zeros(3, length(t));
angle(1, :) = intt(gyro(1, :), t);
angle(2, :) = intt(gyro(2, :), t);
angle(3, :) = intt(gyro(3, :), t);

% Filtering
% init
Q = [1e-1 0 0
    0 1e-1 0
    0 0 1e-1];
R = [1e-1 0 0
    0 1e-1 0
    0 0 1e-1];
P_new = [1e1 0 0
        0 1e1 0
        0 0 1e1];
predict_angle = zeros(3, length(t));
o_acc = zeros(3, length(t));
o_gyro = zeros(3, length(t));
% loop for estimation
for i = 1 : length(t)
    if i == 1
        pre_angle = predict_angle(:, i);
    end
    if i > 1
        pre_angle = predict_angle(:, i-1);
    end 
    % p increase
    [predict_angle(1, i), predict_angle(2, i), predict_angle(3, i), P_new] = ekf_prediction(pre_angle, gyro(:, i), Q, P_new); 
    % p decrease

    if mod(i, 2) == 0
        [predict_angle(1, i), predict_angle(2, i), predict_angle(3, i), P_new, o_acc(:,i), o_gyro(:, i)] = ekf_update(predict_angle(:, i), acc(:, i), R, P_new, mag(:,i));
        
    end
    
end
%% plot
% gyro to attitude
figure (1); 
subplot(2, 2, 1); 
plot(t, angle(1, :), 'b-');   
hold on;    
plot(t, angle(2, :), 'g-');   
plot(t, angle(3, :), 'r-'); 
hold off;    
title('Gyro readings');    
xlabel('Time(s)');    
ylabel('Angle(rad)');    
legend('roll', 'pitch', 'yaw');    
grid on;
% predict attitude
subplot(2, 2, 2); 
plot(t, predict_angle(1, :), 'b-'); 
hold on;    
plot(t,  predict_angle(2, :), 'g-'); 
plot(t,  predict_angle(3, :), 'r-'); 
hold off;    
title('Prediction');    
xlabel('Time(s)');    
ylabel('Angle(rad)');      
grid on;
% o_acc
subplot(2, 2, 3); 
plot(t, o_acc(1, :), 'b-'); 
hold on;    
plot(t,  o_acc(2, :), 'g-'); 
%plot(t,  o_acc(3, :), 'r-'); 
hold off;    
title('acc measurement');    
xlabel('Time(s)');    
ylabel('g/s');      
grid on;
% o_gyro
subplot(2, 2, 4); 
plot(t, o_gyro(1, :), 'b-'); 
hold on;    
plot(t,  o_gyro(2, :), 'g-'); 
%plot(t,  o_gyro(3, :), 'r-'); 
hold off;    
title('gyro prediction of acc measurement');    
xlabel('Time(s)');    
ylabel('g/s');      
grid on;

%% functions

% x = [roll, pitch, yaw]^T
% w = [dot(roll), dot(pitch), dot(yaw)]^T
% Q is 3x3 diagonal noise matrix, P is last estimated cov
function [roll, pitch, yaw, P_new] = ekf_prediction(x, w, Q, P)
    roll = x(1);
    pitch = x(2);
    yaw = x(3);
    dt = 0.01;
    % Jacobian of the state transition model 
    C = [1+dt*w(2)*cos(roll)*tan(pitch)-dt*w(3)*sin(roll)*tan(pitch) dt*w(2)*sin(roll)*sec(pitch)*sec(pitch)+dt*w(3)*cos(roll)*sec(pitch)*sec(pitch) 0
        0 1 0
        0 0 1];
     % C = [1 0 0
     %    0 1 0
     %    0 0 1];
    % Gyro model
    rot_b2e = [1 sin(roll)*tan(pitch) cos(roll)*tan(pitch)
            0 cos(roll) -sin(roll)
            0 sin(roll)/cos(pitch) cos(roll)/cos(pitch)];
    w = rot_b2e * w;
    x = x + dt * w;

    % Estimate covariance
    P_new = C * P * C' + Q;

    % Output
    roll = x(1);
    pitch = x(2);
    yaw = x(3);
end

% x = [roll, pitch, yaw]^T
% e = [roll, pitch, yaw]^T from sensor
% R is 3x3 diagonal noise matrix, P is last cov
function [roll, pitch, yaw, P_updated, o, o_expected] = ekf_update(x, o, R, P, m)
    % Jacobian of observation model
    roll = x(1);
    pitch = x(2);
    yaw = x(3);
    m1 = m(1);
    m2 = m(2);
    m3 = m(3);
    g = 9.81;
    % mag_dr = (m1*m2*sin(roll)*cos(pitch) - m1*m3 - m2^2*sin(pitch)*cos(pitch)*cos(roll) - m2*m3*cos(pitch)^2*cos(roll)) / ...
    % ((m2*cos(pitch) - m3*sin(roll))^2 + (m1*cos(roll) + m2*sin(pitch)*sin(roll) + m3*sin(roll)*cos(pitch))^2);
    % 
    % mag_dp = (-m1*m2*sin(pitch)*cos(roll) - m2^2*sin(roll) + m2*m3*sin(roll)^2*cos(pitch) - m3^2*sin(pitch)*sin(roll)^2) / ...
    % ((m2*cos(pitch) - m3*sin(roll))^2 + (m1*cos(roll) + m2*sin(pitch)*sin(roll) + m3*sin(roll)*cos(pitch))^2);

    % H = [0 -1*cos(pitch)*g 0
    %     cos(pitch)*cos(roll)*g -1*sin(pitch)*sin(roll)*g 0
    %     -cos(pitch)*sin(roll)*g -sin(pitch)*cos(roll)*g 0];
    H = [0 -1*cos(pitch)*g 0
        cos(pitch)*cos(roll)*g -1*sin(pitch)*sin(roll)*g 0
        0 0 1];
    
    % Kalman gain
    K = (P * transpose(H)) / (H *P * transpose(H) + R);
    
    % Updated estimation
    % acc_z to mag_yaw
    mx = m(1) * cos(roll) + m(2) * sin(roll) * sin(pitch) + m(3) * sin(roll) * cos(pitch);
    my = m(2) * cos(pitch) - m(3) * sin(roll);
    yaw_mag = atan2(my, mx) + deg2rad(11);  
    o = o + [0; 0; -o(3)+yaw_mag]
    o_expected = h(x);
    x = x + K*(o - o_expected);
    
    % Update covariance
    P_updated = (eye(3) - K * H) * P;
    
    % Output
    roll = x(1);
    pitch = x(2);
    yaw = x(3);
end

% Input current state, output expected observation (of acc and mag)
function o_expected = h(x)
   
    g = 9.81;
    roll = x(1);
    pitch = x(2);
    yaw = x(3);
    g_stary = [0
               0
               g];
    % rotation
    % acc
    % rot_a2m = [cos(pitch)*cos(yaw) cos(pitch)*sin(yaw) -1*sin(pitch)
    %                sin(roll)*sin(pitch)*cos(yaw)-cos(roll)*sin(yaw) sin(roll)*sin(pitch)*sin(yaw)+cos(roll)*cos(yaw) sin(roll)*cos(pitch)
    %                cos(roll)*sin(pitch)*cos(yaw)+sin(roll)*sin(yaw) cos(roll)*sin(yaw)*sin(pitch)-sin(roll)*cos(yaw) cos(roll)*cos(pitch)];
    rot_a2m = [-1*sin(pitch)*g
                cos(pitch)*sin(roll)*g
                cos(pitch)*cos(roll)*g];
    
    
    o_expected = rot_a2m + [0; 0; -cos(pitch)*cos(roll)*g+yaw];

 

end
