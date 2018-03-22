close all;
% Main routine (driver) for Euler method demonstration 
odefun = 'myode1'; %Specify name of user supplied                    
                   % function M-file with rhs of ode 
t0 = 0; 
tfinal = 2;             % Specify initial and final times 
U0 = [0;0];             % Specify column vector of initial values 
% U0 = [1;2];
NSTEP = [50,100,200];  
% NSTEP = [2000, 4000, 8000, 16000, 32000];
                        % number of step for part e 
                        % Specify number of steps.            
                        % Stepsize Delta_t = (tfinal-t0)/NSTEP.
L = length(NSTEP);           
TSPAN = [t0,tfinal];    
maxerr = zeros(L,1);    % store maximum error for each different step size
delta = zeros(L,3);     % contain differences of U(2) for different step size

delta(:,1) = NSTEP;
for s = 1:L
[t,U] = eulerw17(odefun,TSPAN,U0,NSTEP(s));

delta(s,2) = U(1,end); % save each U(2) with different step size 

% plot numerical solution; 
figure; 
hold off; 
plot(t,U(1,:));

% If the IVP is the demonstration  IVP 
% u'' + t^2u' + y^2 = t^6 + 3t^4 + 6t
% U(0) = U'(0) = 0
% plot the exact solution and find the maximum error 

if odefun == 'myode1' & U0 == [0;0]
    uexact = t.^3;     
    % Compute maximum error     
    maxerr(s) = max(abs(uexact-U(1,:)));     
    % Plot exact solution     
    hold on;     
    plot(t,uexact,'r');     
    xlabel(['red = exact sol., blue = numerical sol.']);     
    title(['NSTEP = ' num2str(NSTEP(s)) ',  max. error = ' ...             
        num2str(maxerr(s))])
end
end

for l = 2:L
delta(l,3) = abs(delta(l , 2) - delta(l - 1, 2)); % compute differences
end