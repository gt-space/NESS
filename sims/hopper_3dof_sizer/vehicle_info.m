mt = 114; 

g = 9.80665; 
C_G = 1;  % center of gravity change over time
m0 = 48;
mp =  mt - m0;
tb= 65;
mdot  = 1.3;

D_m    = 0.2;
r_m = D_m/2;
L_m    = 2.3;


% Closed thin-walled cylinder (side + two caps)
I_xx = 0.5 * mt * r_m^2;
I_yy = (1/12) * mt * (3*r_m^2 + L_m^2);
I_zz = I_yy;

% Inertia Matrix
I = [ I_xx, 0,     0;
      0,    I_yy,  0;
      0,    0,     I_zz ];

% % A matrix (12x12)
% A =
% 
% % B matrix (12x6)
% B = 
% 
% % C matrix (identity - outputs = all states)
% C = eye(7);
% 
% % D matrix (no direct feedthrough)
% D = 0;
% 
% 
% p = [ -2.0, -2.2, -2.4, -2.6, -2.8, -3.0, ...
%       -3.5-2.0j, -3.5+2.0j, ...
%       -4.0-2.5j, -4.0+2.5j, ...
%       -5.0-3.0j, -5.0+3.0j ]';
% 
% K = place(A,B,p);
% 
% 
% disp('K (place) ='); disp(K);
% eig(A-B*K)
% 
% sys = ss(A, B, C, D);
% 
% TF = ss2tf(A, B, C, D)
% 
% 
% 
% 
% %%Trajectory:
% dt = 100
% x = linspace(0, 15, dt)'; 
% z = -x.^2 + 15*x;           
% t= [0:dt-1]';
% input_xz = [x, z];  
% 
% X = [x, t];
% Z = [z, t];
% plot(x, z)
% grid on
% % traj = [X; Z];
% % maxAltitude_m      = 500;
% % dryMass_kg         = 48;
% % wetMass_kg         = 72;
% % diameter_m         = 0.2;
% % TWR_takeoff        = 1.5;
% % maxVelocity_mps    = 26;
% % L_over_D           = 11.3;
% % nominalThrust_N    = 1334.5;
% % totalImpulse_Ns    = 86740.3;
% % g0                 = 9.80665;
% % 
% % % Derived numeric results
% % propellantMass_kg      = wetMass_kg - dryMass_kg;
% % Isp_s                  = totalImpulse_Ns / (propellantMass_kg * g0);
% % TWR_from_thrust        = nominalThrust_N / (wetMass_kg * g0);
% % impulse_from_Tburn_Ns  = nominalThrust_N * burnTime_s;