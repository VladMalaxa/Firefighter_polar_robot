% --- Robot Parameters ---
I1 = 1.5; 
I2 = 0.5; 
l_c1 = 0.5;
m1 = 15;
m2 = 8;
g = 9.81;
alpha = -pi/2;

% --- Damping ---
d1 = 2.0; 
d2 = 5.0; 
D = [d1 0; 0 d2];

% --- Control / Simulation Setup ---
S = eye(2);

% --- Control Gains (Tuned for Stability) ---

% Rotation Joint (q1)
Kp1 = 500; 
Ki1 = 150;   % Small value to start
Kd1 = 250;

% Prismatic Joint (q2)
Kp2 = 400; 
Ki2 = 120;   % Higher for joint 2 to fight gravity
Kd2 = 150;



% Define initial q and dq for your integrators if needed
q_init = [0; 0];
dq_init = [0; 0]; % Initial joint velocities

% position mark
x = 0;
y = 0;
