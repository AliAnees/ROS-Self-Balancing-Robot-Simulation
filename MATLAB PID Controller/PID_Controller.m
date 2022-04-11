%% Linearization - State Space Template File!

clear all;
close all;
clc;

% Defines the symbolic variables for taking differential and partials

syms q1 q2 q1_dot q2_dot q1_double_dot q2_double_dot m_body m_wheel l_body r_wheel g Tau J
q_old = [q1 q2 q1_dot q2_dot m_body m_wheel l_body r_wheel g Tau];

% Constants (Change as required)
m_body_val = 3.12193; % mass of rod (body of robot),in kilograms
m_wheel_val = 0.11752; % mass of wheel, in kilograms
l_body_val = 0.09675;
r_wheel_val = 0.05; % radius of wheel, in meters
g_val = 9.81; % gravity , units kg/ms^2 
Tau_val = 0.1907; % N*m

EoM_1 = (m_body*(l_body*cos(q1)*(r_wheel*q2_double_dot + (l_body*cos(q1)*q1_double_dot)/2 - (l_body*sin(q1)*q1_dot*q1_dot)/2) + (l_body^2*sin(q1)^2*q1_double_dot)/2 - l_body*sin(q1)*(r_wheel*q2_dot + (l_body*cos(q1)*q1_dot)/2)*q1_dot + l_body^2*cos(q1)*sin(q1)*q1_dot*q1_dot))/2 - (m_body*((l_body^2*q1_dot^2*cos(q1)*sin(q1))/2 - l_body*q1_dot*sin(q1)*(q2_dot*r_wheel + (l_body*q1_dot*cos(q1))/2)))/2 + (l_body^2*m_body*q1_double_dot)/3 - (g*l_body*m_body*sin(q1))/2 == -Tau;
EoM_2 = (3*m_wheel*r_wheel^2*q1_double_dot)/2 + m_body*r_wheel*(r_wheel*q2_double_dot + (l_body*cos(q1)*q1_double_dot)/2 - (l_body*sin(q1)*q1_dot*q1_dot)/2) == Tau;
EoMs = [EoM_1 EoM_2];


%%%%%%%%%%%%%%%%%%% System Linearization %%%%%%%%%%%%%%%%%%%%%%%

% Solve for second derivative of the system
q1_dd(q_old) = isolate(EoM_1, q1_double_dot); %rearanges the equation to solve for q1_double_dot and q2_double_dot
q2_dd(q_old) = isolate(EoM_2, q2_double_dot);


a1(q_old) = rhs(diff(q1_dd,q1)); %all partial derivitives that are input into the function
a2(q_old) = rhs(diff(q1_dd,q1_dot));
a3(q_old) = rhs(diff(q1_dd,q2));
a4(q_old) = rhs(diff(q1_dd,q2_dot));

a5(q_old) = rhs(diff(q2_dd,q1));
a6(q_old) = rhs(diff(q2_dd,q1_dot));
a7(q_old)= rhs(diff(q2_dd,q2));
a8(q_old) = rhs(diff(q2_dd,q2_dot));

b1 (q_old) = rhs(diff(q1_dd,Tau));
b2 (q_old) = rhs(diff(q2_dd,Tau));

q1_val = 0;      q2_val = 0; %initial conditions
q1_dot_val = 0;  q2_dot_val = 0;

 %dQ/dt = AQ + Bu
 % state space form we must solve for when stable (all input parameters are
 % 0)

J(q_old) =  [0 1 0 0;    %creates jacobian based on partial derivatives, is a function of q1,q1_dot,q2,q2_dot   
            a1 a2 a3 a4;           
             0 0 0 1;
            a5 a6 a7 a8];
 
 q_in = [q1_val, q1_dot_val, q2_val, q2_dot_val, m_body_val, m_wheel_val, l_body_val, r_wheel_val, g_val, Tau_val];%creates initial values to be input into A and B matrices of stat-space model

 A_lin = J;
 B_lin = [0;
          b1;
          0;
          b2];
%%    
A_lin_val = vpa(subs(A_lin,q_old,q_in)); %inputs the given values into the matrices to creat a linearized matrix 
B_lin_val = vpa(subs(B_lin,q_old,q_in));
     
% To display the matrices: dynamic matrix and control matrix.
A_lin
B_lin
    
%%%%%%%%%%%%%%%%%%% Simulation Parameters %%%%%%%%%%%%%%%%%%%%

% number of measurements
T = 0.001;
tf = 25;     % simulation finish time
% number of states

% Discretization and Simulation of the Model
% Discretized A and B
A_lin_discrete = [1 T 0 0;    %creates dsicretized jacobian     
                  T*a1 (1+T)*a2 T*a3 T*a4;           
                  0 0 1 T;
                  T*a5 T*a6 T*a7 (1+T)*a8];
B_lin_discrete = [0 ;
                T*b1 ;
                  0 ;
                 T*b2];              
% time intialization (Change as required)

A_lin_val_discrete = vpa(subs(A_lin_discrete,q_old,q_in)); %inputs vaslues to solve for the matix (output is symfun)
B_lin_val_discrete = vpa(subs(B_lin_discrete,q_old,q_in)); 

% Manually creat the matrices based on results of line 81 and 85, from
% symfun to 4x4 matrix (for A) and 4x1 (for B)

% discretized values for the A and B maticies as shown in lines 81 and 85
% this was used to calculate and copy the values into A_lin_val_discrete
% and B_lin_val_discrete on lines 113 and 118
A1 = (T*a1(0,0,0,0,3.12196,0.11752,0.09675,0.05,9.81,0.1907));
A2 = ((1+T)*a2(0,0,0,0,3.12196,0.11752,0.09675,0.05,9.81,0.1907));
A3 = T*a3(0,0,0,0,3.12196,0.11752,0.09675,0.05,9.81,0.1907);
A4 = T*a4(0,0,0,0,3.12196,0.11752,0.09675,0.05,9.81,0.1907);
A5 = T*a5(0,0,0,0,3.12196,0.11752,0.09675,0.05,9.81,0.1907);
A6 = T*a6(0,0,0,0,3.12196,0.11752,0.09675,0.05,9.81,0.1907);
A7 = T*a7(0,0,0,0,3.12196,0.11752,0.09675,0.05,9.81,0.1907);
A8 = (1+T)*a8(0,0,0,0,3.12196,0.11752,0.09675,0.05,9.81,0.1907);

B1 = T*b1(0,0,0,0,3.12196,0.11752,0.09675,0.05,9.81,0.1907);
B2 = T*b2(0,0,0,0,3.12196,0.11752,0.09675,0.05,9.81,0.1907);


A_lin_val_discrete = [1.0, 0.001, 0, 0;
                     0.086910299003322259136212624584718, 0, 0, 0;
                     0, 0, 1.0, 0.001;
                     0, 0, 0, 0];

B_lin_val_discrete =[0 ;
                      -0.0587;
                     0 ;
                     0.1281];
A = A_lin_val_discrete;         %simplifying the A and B matricies to make it easier for the PID controlling
B = B_lin_val_discrete;
                
%%%%%%%%%%%%%%%%%%%% PID setup %%%%%%%%%%%%%%%%%%%%%%%%


t = T:T:tf;         %time vector

% PID Gain Initialization (Kp, Kd, Ki?)
Kp = -13;
Kd = -5; 
Ki = -1;
% Initialization Parameters for Error calculations

e = zeros(1, length(t));    %initializing the error for the angular pos. of the body
sum_error = e(1);           %initializing the summation of the area under the error curve
xd = zeros(1, length(t));   %initializing the desired angular body position
xd(tf/T/10:end) = xd(tf/T/10:end) + 1;  % superimpose +1 starting at 1/10th of the simulation

%% Simulating System dynamics 
C = eye(4);                     % measurement matrix
n = 4;                  % number of states
m = 4;                  % number of measurements
x = zeros(n, length(t));     % initialize state equation
z = zeros(m, length(t));     % initialize measurement z to 0
u = zeros(1, tf/T);          % initializing input

for k = 2:length(t)-1               % Generates the true state trajectories and measurements
   u(k) = xd(k)-x(1, k);            % setting the input to be the error (Desired - Actual)
   x(:, k+1) = A_lin_val_discrete*x(:, k) + B_lin_val_discrete*(-u(k));     % state equation
   z(:, k+1) = C*x(:, k+1);            % measurement equation
end


%% Simulation No Controller
figure; 
plot(T:T:tf, xd(:)); xlabel('Time (sec)'); ylabel('Angular Position of the Body (rad)'); title('Closed-Loop Without Controller');     % plot desired position
hold on
plot(T:T:tf, x(1, :));     % plot actual position
legend('Desired Position', 'Actual Position');
hold off
figure; plot(T:T:tf, x(2, :)); xlabel('Time (sec)'); ylabel('Angular Velocity of the Body (rad/sec)'); title('Closed-Loop Without Controller'); % plot velocity
figure; plot(T:T:tf, x(3, :)); xlabel('Time (sec)'); ylabel('Angular Position of the Wheel (rad)'); title('Closed-Loop Without Controller'); % plot acceleration
figure; plot(T:T:tf, x(4, :)); xlabel('Time (sec)'); ylabel('Angular Velocity of the Wheel (rad/sec)'); title('Closed-Loop Without Controller'); % plot acceleration
figure; plot(T:T:tf, u(1, :)); xlabel('Time (sec)'); ylabel('Angular Position (rad)'); title('Input Without PID Controller'); % plot input

%% Simulation PID Controller

x = zeros(n, length(t));     % initialize state equation
u = zeros(1, tf/T);          % intitiate input
for k = 2:length(t)-1
    e(k) = xd(k)-x(1, k);                   % calculating error between the current desired position and the current body angle position
    sum_error = sum_error + e(k)*T;         % calculating culmutive error using trapazoid method
    u(k) = Kp*e(k) + Ki*sum_error + Kd*(e(k)-e(k-1));      % calculate controller input using the PID gains and error
    x(:,k+1) = A*x(:, k) + B*u(k);                      % state equation
    z(:,k+1) = C*x(:, k+1);                             % measurement equation
end
RMSE = sqrt(sum(e.^2)/length(e));                       %calculate Root mean Square Error to check how good PID is working, good is below 0.5

disp(['RSME Results:    ', num2str(RMSE)]);             %display root mean square error
figure; 
plot(T:T:tf, xd(:)); xlabel('Time (sec)'); ylabel('Angular Position of the Body (rad)'); title('Closed-Loop With PID Controller');     % plot desired position
hold on
plot(T:T:tf, x(1, :));     % plot position
legend('Desired Position', 'Actual Position');
hold off
figure; plot(T:T:tf, x(2, :)); xlabel('Time (sec)'); ylabel('Angular Velocity of the Body (rad/sec)'); title('Closed-Loop With PID Controller'); % plot velocity
figure; plot(T:T:tf, x(3, :)); xlabel('Time (sec)'); ylabel('Angular Position of the Wheel(rad)'); title('Closed-Loop With PID Controller');     % plot acceleration
figure; plot(T:T:tf, x(4, :)); xlabel('Time (sec)'); ylabel('Angular Velocity of the Wheel (rad/s)'); title('Closed-Loop With PID Controller');  % plot acceleration
figure; plot(T:T:tf, u(1, :)); xlabel('Time (sec)'); ylabel('Angular Position (rad)'); title('PID Controller Input');                          % plot input

