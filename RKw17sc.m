function [T,U, count] = RKw17sc(odefun, TSPAN ,U0,TOL,A,b,c,qorder);
%RK method with stepsize control using Richardson extrapolation
% qorder = order of method used

count = 0; %count number of evaluations 

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
    while flag == 1; 
        %one big step with dt
        Ubig = RKexplicitstep(odefun,t,U(:,kstep),dt,A,b,c); 
        count = count + 1;
        
        % two small steps with dt/2
        dt2 = dt/2;
        V = RKexplicitstep(odefun,t,U(:,kstep),dt2,A,b,c);
        Usmall = RKexplicitstep(odefun,t+dt2,V,dt2,A,b,c);
        count = count + 1;
        
        %Richardson extrapolation
        fac = 2^qorder; 
        Unew = (fac*Usmall-Ubig)/(fac-1);  
        count = count + 1;
        
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