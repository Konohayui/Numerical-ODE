function [t,U] = eulerw17_3(odefun, TSPAN ,U0, NSTEP, A, b, c, d)
T0 = TSPAN(1); 
TFINAL = TSPAN(2); 
dt = (TFINAL - T0)/NSTEP; 
m = length(U0); 
U = zeros(m,NSTEP+1); 
U(:,1) = U0; 
t = T0:dt:TFINAL; 

if d == 0
for k = 1:NSTEP     
    U(:,k+1) = eulerstep(odefun,t(k),U(:,k),dt); 
end
end

if d == 1
for k = 1:NSTEP     
    U(:,k+1) = RKexplicitstep(odefun,t(k),U(:,k),dt,A,b,c);
end
end

if d == 2
[T,U] = RKw17sc(odefun, TSPAN, U0, 10^-2, A, b, c, 4);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function U = eulerstep(odefun,t,U0,dt) 
U = U0 + dt*feval(odefun,t,U0); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U = RKexplicitstep(odefun,t,U0,dt,A,b,c) 
% One step of a general explicit Runge Kutta method 
% A is assumed to be strictly lower triangular 
% b and c must be column vectors 
% 
r = length(b); 
m = length(U0); 
K = zeros(m,r);
%matrix whose j-th column is K_j 
K(:,1) = feval(odefun,t,U0); 
for j = 2:r     
    Y = U0 + dt*K(:,1:j-1)*(A(j,1:j-1).');      
    % Y = U0 + dt*sum_{l=1}^{j-1} a_{jl}K_l     
    K(:,j) = feval(odefun,t + c(j)*dt, Y); 
end
U = U0 + dt*K*b; 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T,U] = RKw17sc(odefun, TSPAN , U0, TOL, A, b, c, qorder)
%RK method with stepsize control using Richardson extrapolation
% qorder = order of method used

dtmin = 1.e-5; % minimum allowed step size
sf = 0.8; %safety factor   

q1 = qorder+1; 
T0 = TSPAN(1); TFINAL = TSPAN(2);
dt = (TFINAL-T0)/10;
m = length(U0);
U = U0;
T(1) = T0;
kstep = 1; 
while T(kstep) < TFINAL
    t = T(kstep);
    flag = 1; 
    while flag == 1
        %one big step with dt
        Ubig = RKexplicitstep(odefun,t,U(:,kstep),dt,A,b,c); 
        
         % two small steps with dt/2
        dt2 = dt/2;
        V = RKexplicitstep(odefun,t,U(:,kstep),dt2,A,b,c);
        Usmall = RKexplicitstep(odefun,t+dt2,V,dt2,A,b,c);
        
        %Richardson extrapolation
        fac = 2^qorder; 
        Unew = (fac*Usmall-Ubig)/(fac-1);  
   
        
        %Estimate of one-step error for Usmall method 
        locerr =norm(Unew-Usmall,1);
         
        % conditions for accepting current dt 
        if (locerr <= TOL ) | (dt < 1.001*dtmin) 
            flag = 0;
            kstep = kstep+1;
            U(:,kstep) = Unew; %local extrapolation 
            T(kstep) = t + dt;
        end
        dt = max(sf*((TOL/locerr)^(1/q1))*dt,dtmin); %new value for dt 
        dt = min(dt,TFINAL-t);
    end
   
  
end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function yprime = myode1(t,y)  
yprime = zeros(size(y));

mu = 0.012277471;
muh = 1 - mu;

D1 = ((y(1) + mu).^2 + y(2)).^(3/2);
D2 = (((y(1) - muh)).^2 + y(2)).^(3/2);

yprime(1) = y(3);
yprime(2) = y(4);
yprime(3) = y(1) + 2*y(4) - muh*(y(1) + mu)./D1 - mu*(y(1) - muh)./D2;
yprime(4) = y(2) - 2*y(3) - muh*y(2)./D1 - mu*y(2)./D2;
end