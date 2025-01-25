
%--------------------------------------------------------------------------------------------------
% Main program to read IMU data, run the fuction that estimates quaternion/Euler angles, and biases 
% and plots IMU data, external acceleration and estimated Euler angles
%--------------------------------------------------------------------------------------------------

R2D = 180/pi;                           % Conversion from rad to deg

% Question 2

% Load IMU data from the file
    load('IMU_Data.mat');

    % acc: 3d external acceleration Accelerometer measurements (3x333 matrix), unit: m/s^2
    % ba: Accelerometer bias (1x3 vector), unit: m/s^2
    % bg: Gyroscope bias (1x3 vector), unit: rad/s
    % euler: Euler angles (3x333 matrix), unit: radians %verify
    % N: Number of measurements (scalar)
    % tt: Time vector (333x1 matrix), unit: seconds
    % wa: %noises Accelerometer measurements in the body frame (3x333 matrix), unit: m/s^2
    % wg: % noisesGyroscope measurements in the body frame (3x333 matrix), unit: rad/s
    % wm: Magnetometer measurements in the body frame (3x333 matrix), unit: Gauss
    % ya: Accelerometer measurements in the navigation frame (3x333 matrix), unit: m/s^2
    % yg: Gyroscope measurements in the navigation frame (3x333 matrix), unit: rad/s
    % ym: Magnetometer measurements in the navigation frame (3x1903 matrix), unit: Gauss

% Question 3
figure;

   % Subplot 1: Accelerometer Measurements X
    subplot(3, 3, 1);
    plot(tt(1:333), ya(1, 1:333), 'b');
    xlabel('');
    ylabel('ya_x (m/s^2)');
    title('3D acceleration');
    grid on;

   % Subplot 2: Gyroscope Measurements X
    subplot(3, 3, 2);
    plot(tt(1:333), yg(1, 1:333), 'b');
    xlabel('');
    ylabel('yg_x (rad/s)');
    title('3D angular velocity');
    grid on;

   % Subplot 3: Magnetometer Measurements X
    subplot(3, 3, 3);
    plot(tt(1:333), ym(1, 1:333), 'b');
    xlabel('');
    ylabel('ym_x (Gauss)');
    title('3D magnetic field');
    grid on;

   % Subplot 4: Accelerometer Measurements Y
    subplot(3, 3, 4);
    plot(tt(1:333), ya(2, 1:333), 'b');
    xlabel('');
    ylabel('ya_yn (m/s^2)');
    title('');
    grid on;

   % Subplot 5: Gyroscope Measurements Y
    subplot(3, 3, 5);
    plot(tt(1:333), yg(2, 1:333), 'b');
    xlabel('');
    ylabel('yg_y (rad/s)');
    title('');
    grid on;
    
   % Subplot 6: Magnetometer Measurements Y
    subplot(3, 3, 6);
    plot(tt(1:333), ym(2, 1:333), 'b');
    xlabel('');
    ylabel('ym_y (Gauss)');
    title('');
    grid on;
    
   % Subplot 7: Accelerometer Measurements Z
    subplot(3, 3, 7);
    plot(tt(1:333), ya(3, 1:333), 'b');
    xlabel('Time (s)');
    ylabel('ya_z (m/s^2)');
    title('');
    grid on;

   % Subplot 8: Gyroscope Measurements Z
    subplot(3, 3, 8);
    plot(tt(1:333), yg(3, 1:333), 'b');
    xlabel('Time (s)');
    ylabel('yg_z (rad/s)');
    title('');
    grid on;

   % Subplot 9: Magnetometer Measurements Z
    subplot(3, 3, 9);
    plot(tt(1:333), ym(3, 1:333), 'b');
    xlabel('Time (s)');
    ylabel('ym_z(Gauss)');
    title('');
    grid on;

% Figure for External Acceleration
    figure;
    % Subplot 1: External Acceleration X
    subplot(3, 1, 1);
    plot(tt(1:333), acc(1, 1:333), 'b');
    xlabel('');
    ylabel('a (m/s^2)');
    title('3D external acceleration');
    grid on;

    % Subplot 2: External Acceleration Y
    subplot(3, 1, 2);
    plot(tt(1:333), acc(2,1:333), 'b');
    xlabel('');
    ylabel('a_{b,2} (m/s^2)');
    title('');
    grid on;

    % Subplot 3: External Acceleration Z
    subplot(3, 1, 3);
    plot(tt(1:333), acc(3, 1:333), 'b');
    xlabel('Time (s)');
    ylabel('a_{b,3} (m/s^2)');
    title('');
    grid on;


% Question 4
    Ra = 0.0056 * eye(3); 
    Rg = 0.003 * eye(3);   
    Rm = 0.001 * eye(3); 

% Question 5

[q4, eulercom4, bahat, bghat] = Compute_Attitude_student(yg,ya,ym,tt,Rg,Ra,Rm);

%Inputs: yg,ya,ym,tt,Rg,Ra,Rm
%Outputs: 
        % q4 - Quaternion representing attitude,
        % eulercom4 - Euler angles for attitude, 
        % bahat - Estimated accelerometer bias, 
        % bghat - bghat
        
% Mean of each quaternion element
m0=mean(q4(1,:))
m1=mean(q4(2,:))
m2=mean(q4(3,:))
m3=mean(q4(4,:))

% Question 26

figure(3)
subplot(3,1,1);
plot(tt, eulercom4(1,:),tt, euler(1,:))
title("Estimated Euler angles");
ylabel('pitch (phi)');
grid on

subplot(3,1,2);
plot(tt, eulercom4(2,:),tt, euler(2,:))
ylabel('roll (theta)');
grid on

subplot(3,1,3);
plot(tt, eulercom4(3,:),tt, euler(3,:))
xlabel('Time(sec)');
ylabel('yaw');
grid on







