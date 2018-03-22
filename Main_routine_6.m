global fcount;
fcount = 0;

C0 = [0.06; 3.3*10^(-7); 5.01*10^(-11); 0.03; 2.4*10^(-8)];

odefun = 'homework6ode';

t0 = 0;
tf = 1;
options = odeset('RelTol',10^(-10),'AbsTol',10^(-7));
[t,y] = ode15s(odefun, [t0, tf], C0, options);

figure
for i = 1:size(y,2)
subplot(2,3,i);
plot(t,y(:,i));
end
