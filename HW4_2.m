close all;
odefun = 'myode2';

T0 = -1; % initial time
TF = 1; % final time

TSPAN = [T0, TF];

U0 = [-0.001]; % initial condition

% Butcher Array from assignment 3 has order 4
A = [0, 0, 0, 0; 1/2, 0, 0, 0;...
      0, 1/2, 0, 0; 0, 0, 1 0];
b = [1/6; 2/6; 2/6; 1/6];
c = [0; 1/2; 1/2; 1];

for k = 4:10
    [t, U] = RKw17sc(odefun, TSPAN, U0, 10^-k, A, b, c, 4);
    hold on;
    plot(t,U(1,:));
end