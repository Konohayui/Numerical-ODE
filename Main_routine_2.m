close all;
% Main routine (driver) for Euler method demonstration 
odefun = 'myode1'; %Specify name of user supplied                    
                   % function M-file with rhs of ode 
t0 = 0; 
tfinal = 1;             % Specify initial and final times 

U0 = [0;1];             
% Specify column vector of initial values 

NSTEP = [500];       % Specify number of steps.            
                     % Stepsize Delta_t = (tfinal-t0)/NSTEP.
           
TSPAN = [t0,tfinal];    

% Butcher Array
A = [0, 0, 0, 0;1/2, 0, 0, 0;...
    0, 1/2, 0, 0;0, 0, 1 0];
b = [1/6;2/6;2/3;1/6];
c = [0;1/2;1/2;1];

[t,U] = eulerw17_3(odefun,TSPAN,U0,NSTEP,A,b,c);

% plot numerical solution; 
figure; 
hold off; 
plot(U(1,:),U(2,:));
xlabel('U_1');
ylabel('U_2');
axis([-1.5 1 -1.5 1.5]);