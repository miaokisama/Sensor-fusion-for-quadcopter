clear;
%% Read data
csvFile1 = 'sensor_data/gyro_data.csv';
csvFile2 = 'sensor_data/acc_data.csv';
csvFile3 = 'sensor_data/mag_data.csv';
data1 = readtable(csvFile1);
data2 = readtable(csvFile2);
data3 = readtable(csvFile3);

% Data Filter
% Table
t = transpose(data3{(2:end), 2});
l = length(t);
for  i = 3 : 5
    meas_f = 100; % sample freq
    Fc = 1; % cutoff frequency
    [Blp, Alp] = cheby2(4, 100, Fc / meas_f); 
    % X1 =  data1{(2: l+1), i};
    % X2 =  data2{(2: l+1), i};
    % X3 =  data3{(2: l+1), i};
    X1 = filtfilt(Blp, Alp, data1{(2: l+1), i});
    X2 = filtfilt(Blp, Alp, data2{(2: l+1), i});
    X3 = filtfilt(Blp, Alp, data3{(2: l+1), i});
    if i == 3 
        gyro(1, :) = X1;
        acc(1, :) = X2;
        mag(1, :) = X3;
    end
    if i == 4
        gyro(2, :) = X1;
        acc(2, :) = X2;
        mag(2, :) = X3;
    end
    if i == 5
        gyro(3, :) = X1;
        acc(3, :) = X2;
        mag(3, :) = X3;
    end  
end
gyro = deg2rad(gyro);
%% Cf loop
% Gyro to angle
angle = zeros(3, length(t));
angle(1, :) = intt(gyro(1, :), t);
angle(2, :) = intt(gyro(2, :), t);
angle(3, :) = intt(gyro(3, :), t);

% Acc/mag to angle
acmag = zeros(3, length(t));
filered_angle = zeros(3, length(t));
compensatory = zeros(3, length(t));
alpha = 0.005;
for i = 1 : length(t)
    acmag(1, i) = atan(acc(2, i) / acc(3, i));
    acmag(2, i) = asin(-1 * acc(1, i) / 9.81);
    pitch = acmag(2, i);
    roll = acmag(1, i);
    % Mag
    mx = mag(1, i) * cos(roll) + mag(2, i) * sin(roll) * sin(pitch) + mag(3, i) * sin(roll) * cos(pitch);
    my = mag(2, i) * cos(pitch) - mag(3, i) * sin(roll);
    acmag(3, i) = atan2(my, mx) + deg2rad(11);  
    
    % CF loop
    gyro_val = gyro(:, i);
    accmag_val = acmag(:, i);
    if i == 1
        pre_angle = filered_angle(:, i);
        dt = t(1);
    end
    if i > 1
        pre_angle = filered_angle(:, i - 1);
        dt = t(i) - t(i - 1);
    end
    [filered_angle(:, i), compensatory(:, i)] = cf_fusion(pre_angle, gyro_val, accmag_val, alpha, dt); 

end

%% Plot
figure (2); 
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
% Acmag
subplot(2, 2, 2); 
plot(t, acmag(1, :), 'b-'); 
hold on;    
plot(t,  acmag(2, :), 'g-'); 
plot(t,  acmag(3, :), 'r-'); 
hold off;    
title('AcMag observation');    
xlabel('Time(s)');    
ylabel('Angle(rad)');    
grid on;
% Filtered  
subplot(2, 2, 3); 
plot(t, filered_angle(1, :), 'b-'); 
hold on;    
plot(t,  filered_angle(2, :), 'g-'); 
plot(t,  filered_angle(3, :), 'r-'); 
hold off;    
title('Filtered');    
xlabel('Time(s)');    
ylabel('Angle(rad)');      
grid on;

subplot(2, 2, 4); 
plot(t, compensatory(1,:), 'b-'); 
hold on;    
plot(t,  compensatory(2,:), 'g-'); 
plot(t,  compensatory(3,:), 'r-'); 
hold off;    
title('Compensatory');    
xlabel('Time(s)');    
ylabel('Angle(rad)');      
grid on;

%% Functions
% Integrate
function x_i = intt(x, t)
    x_i = zeros(1,length(t));
    for i = 2:length(t)-1
        x_i(i) = ((x(i-1)+x(i))*0.5+ (x(i+1)+x(i))*0.5)*0.5*((t(i+1)+t(i))/2-(t(i-1)+t(i))/2)+ x_i(i-1);
    end
    x_i(1) = x_i(2);
    x_i(length(t)) = x_i(length(t)-1);
end

% Cf fusion
function  [state, compensatory] = cf_fusion(x, sensor1_val, sensor2_val, alpha, dt)
    % init
    x_filtered = zeros(3, 1);
    rot_b2e = zeros(3, 3, 1);
    roll = x(1);
    pitch = x(2);
    yaw = x(3);

    rot_b2e = [1 sin(roll)*tan(pitch) cos(roll)*tan(pitch)
        0 cos(roll) -sin(roll)
        0 sin(roll)/cos(pitch) cos(roll)/cos(pitch)];
    w = rot_b2e * sensor1_val;
 
    % fuse
    x_filtered = x + w * dt + alpha * (sensor2_val - x);  
    state = x_filtered;
    compensatory = alpha * (sensor2_val - x);  
end


