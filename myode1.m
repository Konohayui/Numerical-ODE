function yprime = myode1(t,y)  
global count
yprime = zeros(size(y));

mu = 0.012277471;
muh = 1 - mu;

D1 = ((y(1) + mu).^2 + y(2).^2).^(3/2);
D2 = (((y(1) - muh)).^2 + y(2).^2).^(3/2);

yprime(1) = y(3);
yprime(2) = y(4);
yprime(3) = y(1) + 2*y(4) - muh*(y(1) + mu)./D1 - mu*(y(1) - muh)./D2;
yprime(4) = y(2) - 2*y(3) - muh*y(2)./D1 - mu*y(2)./D2;
count = count + 1;
end