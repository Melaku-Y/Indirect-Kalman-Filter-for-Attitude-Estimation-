% -------------------------------------------------------------------
% Attitude estimation using a quaternion-based indirect Kalman filter
% -------------------------------------------------------------------

function [q4, eulercom4, bahat, bghat] = Compute_Attitude_student(yg,ya,ym,tt,Rg,Ra,Rm)

D2R = pi/180;                               % Conversion from deg to rad

% Number of data: N
N = max(size(ya));

% Question 6
    alpha = 50;  % Dip angle in degrees
    alpha_rad = alpha * D2R;  % Convert dip angle to radians
  
    g= [0; 0; 9.81];  % Gravity vector in NED frame (z-axis points upwards)
    m= [cos(alpha_rad); 0; -sin(alpha_rad)];  % Magnetic field vector in NED frame
        %cos(alpha_rad) is the horizontal component (North direction).
        %0 for the East component.
        %-sin(alpha_rad) is the downward (negative z) component in the NED frame.

% Question 7
    Qba=0.000001*eye(3); %process noise covariance for the accelerometer bias. The small value 0.000001 suggests low process noise, indicating that the accelerometer bias is expected to change slowly over time.
    Qbg=0.000001*eye(3); %process noise covariance for the gyroscope bias. Similarly, the small value implies that the gyroscope bias is also expected to be relatively stable with minimal noise over time.
    Q=blkdiag(0.25*Rg,Qbg,Qba); %process noise covariance matrix

% --------------------------------------------------------
% Kalman Filter
% --------------------------------------------------------

% q4: quaternion
q4 = zeros(4,N);

% eulercom4: euler angles
eulercom4 = zeros(3,N);

% Estimated bias for gyroscope (bghat) and accelerometer (bahat) 
bghat = zeros(3,1);
bahat = zeros(3,1);

% inital orientation estimation using the TRIAD method
yabar = ya(:,1) / norm(ya(:,1));
ymbar = ym(:,1) / norm(ym(:,1));

foo1 = cross(yabar,ymbar) / norm( cross(yabar,ymbar) );
C = [ -cross(yabar,foo1)  , foo1 , yabar ] ;
q4(:,1) = dcm2quaternion_student(C);

% Kalman filter state creation
x = zeros(9,1);

% Question 8
    P=blkdiag(0.01*eye(3), 0.000001*eye(3), 0.000001*eye(3)); %initial state error covariance matrix

% Question 9
    wx = yg(1,1);  % Angular velocity around the x-axis
    wy = yg(2,1);  % Angular velocity around the y-axis
    wz = yg(3,1);  % Angular velocity around the z-axis

        omega_init = [0, -wx, -wy, -wz; 
                      wx, 0, wz, -wy; 
                      wy, -wz, 0, wx; 
                      wz, wy, -wx, 0];

% variable used in the adaptive algorithm      
r2count = 100;

% H1 and H2 for measurement update
Ha = zeros(3,9);
Ha(1:3,7:9) = eye(3);
Hm =  zeros(3,9);

% parameter for adaptive algorithm
M1 = 3;
M2 = 3;
gamma = 0.1;

R = zeros(3,3*N);

% Kalman filter loop
for i = 2:N
    % Question 10
        %Ts=tt(2)-tt(1); %sampling time
        Ts = tt(i)-tt(i-1);
    
    % Question 11
    % p1=yg(1,i)-bghat(1);
        % p2=yg(2,i)-bghat(2);
        % p3=yg(3,i)-bghat(3);
        % yg_skew =[0, -p3, p2;
        %           p2, 0 , -p1;
        %          -p2, p1, 0];
        % 
        % A_kt=[-yg_skew(1,:) ,    -0.5, 0, 0,   0 0 0;
        %       -yg_skew(2,:) ,    0, -0.5, 0,   0 0 0;
        %       -yg_skew(3,:) ,    0,   0, -0.5, 0 0 0;
        %       0 0 0              0    0    0   0 0 0;
        %       0 0 0              0    0    0   0 0 0;
        %       0 0 0              0    0    0   0 0 0;
        %       0 0 0              0    0    0   0 0 0;
        %       0 0 0              0    0    0   0 0 0;
        %       0 0 0              0    0    0   0 0 0];
          
A_kt = [-vec2product_student(yg(:,i-1) - bghat) -0.5*eye(3) zeros(3);zeros(6,9)];
    
    % Question 12
    % Compute the transition matrix  from Taylor expansion 2
        phi_k = eye(9) + A_kt * Ts + 0.5 * A_kt^2 * Ts^2;  % Taylor expansion of order 2

    x = phi_k * x;
    
    % Question 13
        Qd = Q*Ts + 0.5* A_kt * Q* Ts^2 + 0.5* Q * A_kt' * Ts^2;
        %Qd = Discretized process noise covariance matrix
        %Q  = Continuous-time process noise covariance matrix
        %Ts = Sampling period (time step between updates)
        %A_kt = System matrix (state transition matrix)
        %A_kt = Transpose of the system matrix A_kt
    
    % Question 14
        Pk= phi_k*P*phi_k'+Qd;
        %Pk is The updated state error covariance matrix at time step
        %Phi_k is the state transition matrix at time step k,
        %P is the updated state error covariance matrix at time step k,
        %Qd is the discretized process noise covariance matrix.

    % Question 15
       
        wx = yg(1, i-1) - bghat(1);
        wy = yg(2, i-1) - bghat(2);
        wz = yg(3, i-1) - bghat(3);
        omega_old = [0, -wx, -wy, -wz; 
                      wx, 0, wz, -wy; 
                      wy, -wz, 0, wx; 
                      wz, wy, -wx, 0];
        wk = yg(:,i)-bghat;
    % Question 16
    q4(:,i)= (eye(4) + (3/4)*omega_old*Ts -(1/4)*omega_init*Ts - (1/6)*norm(wk,2)^2*(Ts^2)*eye(4)-(1/24)*omega_init*omega_old*Ts^2-(1/48)*norm(wk)^2*omega_old*Ts^3)*q4(:,i-1);
    omega_init = omega_old;
    Cq = quaternion2dcm_student(q4(:,i));
     
    % ----------------------------------------------------
    % two step measurement update
    % ----------------------------------------------------
    
    % Question 17
    
    Ha=[2*vec2product_student(Cq*g), zeros(3), eye(3)];
    za=(ya(:,i)-bahat)- Cq*g;

    
    % adaptive algorithm
    fooR1 = (za - Ha*x) * (za - Ha*x)';
    R(:,3*(i-1)+1:3*i) = fooR1;
    uk = fooR1;
    for j = i-1:min([i-(M1-1),1])
        uk = uk + R(:,3*(j-1)+1:3*j);
    end
    uk = uk / M1;
    fooR2 = Ha*P*Ha' + Ra;
    [u,s,v] = svd(uk);
    u1 = u(:,1);
    u2 = u(:,2);
    u3 = u(:,3);
    lambda = [ s(1) , s(2) , s(3)];
    mu =  [ u1' * fooR2 * u1 , u2' * fooR2 * u2 , u3' * fooR2 * u3];
    if ( max(lambda - mu) > gamma )
      r2count = 0;
      Qa_b = max(lambda(1)-mu(1),0)*u1*u1' + max(lambda(2) -mu(2),0)*u2*u2' + max(lambda(3)-mu(3),0)*u3*u3';
    else
      r2count = r2count + 1;
      if ( r2count < M2 )
        Qa_b = max(lambda(1)-mu(1),0)*u1*u1' + max(lambda(2) -mu(2),0)*u2*u2' + max(lambda(3)-mu(3),0)*u3*u3';
      else
        Qa_b = zeros(3,3);
      end
    end
    
    % Question 18
    
     Ka = P*Ha' * inv(Ha*P*Ha' + Ra + Qa_b);
    % Question 19
   
    x=x+Ka*(za-Ha*x);
    % Question 20
    P = (eye(9)-Ka*Ha)* P * (eye(9)-Ka*Ha)' + Ka*(Ra+Qa_b)*Ka';
    
    % Question 21
    qe_h=[1;x(1:3)];
    q4(:,i) = quaternionmul_student(q4(:,i),qe_h);
    q4(:,i) = q4(:,i)/norm(q4(:,i));
    Cq = quaternion2dcm_student(q4(:,i));
    x(1:3) = 0;
    %%%%%%%%%%%%%
    % Question 22
    
     Hm = [2*vec2product_student(Cq*m), zeros(3), zeros(3)];
    zm = ym(:,i)-Cq*m;

    % Question 23
    
    Pm = [P(1:3,1:3) zeros(3,6);
            zeros(6,3) zeros(6,6)];
    r3 = Cq*[0;0;1];
    
    Km = [r3*r3' zeros(3,6);
          zeros(6,3) zeros(6)]*Pm*Hm'*inv(Hm*Pm*Hm'+Rm);
    x = x + Km*(zm-Hm*x);
    P = (eye(9)-Km*Hm)*P*(eye(9) - Km*Hm)' + Km*Rm*Km';
 
    bghat = bghat + x(4:6);
    x(4:6) = zeros(3,1);
    
    bahat = bahat + x(7:9);
    x(7:9) = zeros(3,1);
    
    
    % Question 24
    qe = [1;x(1:3)];
    q4(:,i) = quaternionmul_student(q4(:,i),qe);
    q4(:,i) = q4(:,i)/norm(q4(:,i),2);
    Cq = quaternion2dcm_student(q4(:,i));
    x(1:3) = 0;
    
    % Question 25
    eulercom4(:,i) = quaternion2euler_student(q4(:,i))/D2R;
     % quaternion2euler
     
end
