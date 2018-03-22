close all;
odefun = 'myode1'; %Specify name of user supplied                    
                   % function M-file with rhs of ode 
t0 = 0; 
tfinal = 17.1;          % Specify initial and final times 

U0 = [0.994;0;0;-2.00158510637908252240537862224];          
                        % Specify column vector of initial values 
           
TSPAN = [t0,tfinal];    

% Butcher Array from assignment 3 has order 4
A1 = [0, 0, 0, 0; 1/2, 0, 0, 0;...
      0, 1/2, 0, 0; 0, 0, 1 0];
b1 = [1/6; 2/6; 2/6; 1/6];
c1 = [0; 1/2; 1/2; 1];

% Butcher Array from assignment 2 has order 2
A2 = [0, 0; 1, 0];
b2 = [1/2; 1/2];
c2 = [0; 1];

% Butcher Array for Euler method
A3 = [0];
b3 = [1];
c3 = [0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[t, U1, count1] = RKw17sc(odefun, TSPAN, U0, 10^-5, A1, b1, c1, 4);
figure;
hold off;
plot(U1(1,:),U1(2,:));

[t, U2, count2] = RKw17sc(odefun, TSPAN, U0, 10^-5, A2, b2, c2, 2);
figure;
hold off;
plot(U2(1,:),U2(2,:));

[t, U3, count3] = RKw17sc(odefun, TSPAN, U0, 10^-5, A3, b3, c3, 1);
figure;
hold off;
plot(U3(1,:),U3(2,:));

count = [count1; count2; count3];

[t,y] = ode45v4(odefun, t0, tfinal, U0, 0);