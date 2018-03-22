function [t,U] = eulerw17_3(odefun, TSPAN ,U0, NSTEP,A,b,c)
T0 = TSPAN(1); 
TFINAL = TSPAN(2); 
dt = (TFINAL - T0)/NSTEP; 
m = length(U0); 
U = zeros(m,NSTEP+1); 
U(:,1) = U0; 
t = T0:dt:TFINAL; 

if size(A,1) = 0 
for k = 1:NSTEP     
    U(:,k+1) = RKexplicitstep(odefun,t(k),U(:,k),dt); 
end
end 
if size(A,1) ~= 0
for k = 1:NSTEP     
    U(:,k+1) = eulerstep(odefun,t(k),U(:,k),dt); 
end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function U = eulerstep(odefun,t,U0,dt) 
U = U0 + dt*feval(odefun,t,U0); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U = RKexplicitstep(odefun,t,U0,dt, A, b, c) 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function yprime = myode1(t,y)  
yprime = zeros(size(y));
yprime(1) = y(2);
yprime(2) = -y(1);
end