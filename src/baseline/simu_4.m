clc;
clear;
close all;

h= -0.1;

B=[10 -10 0 0;
    -10 20 -10 0;
    0 -10 20 -10;
    0 0 -10 10];
L = [2 -1 0 -1;
     -1 2 -1 0;
     0 -1 2 -1;
     -1 0 -1 2];
D = diag([1 2 4 8]);

A= h*B*D*L*D;



B = [ 1 1 1 1]';
C = [1 0 0 0];
D_d = [0];

sys = ss(A,B,C,D_d);

% Time vector
dt = 0.001;
t = 0:dt:10;

% Approximate impulse at t = 1
u = zeros(size(t));

% Pulse width = dt, area = 1
u(t == 1) = 1/dt;

% Simulate
[y,t,x] = lsim(sys,u,t);

x_n = x ./ [8 4 2 1];
% Plot
plot(t,x_n,'LineWidth',1.5);
xlabel('Time (s)');
ylabel('States');
grid on;