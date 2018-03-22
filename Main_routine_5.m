close all;
global count
count = 0; % initialize 
odefun = 'myode1'; %Specify name of user supplied                    
                   % function M-file with rhs of ode 
t0 = 0; 
tfinal = 17.1;          % Specify initial and final times 

U0 = [0.994;0;0;-2.00158510637908252240537862224];          
                        % Specify column vector of initial values 

[t,y] = ode45v4(odefun, t0, tfinal, U0, 10^-5, 0);

figure;
plot(y(:,1), y(:,2));
title('Orbit');