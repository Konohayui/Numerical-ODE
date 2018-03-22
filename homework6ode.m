function cprime = homework6ode(t,c)
global fcount;
cprime = zeros(size(c));
k = [1.34; 1.6*10^9; 8*10^3; 4.0*10^7; 1.0];

cprime(1) = -k(1)*c(1)*c(2) - k(3)*c(1)*c(3);
cprime(2) = -k(1)*c(1)*c(2) - k(2)*c(2)*c(3) + k(5)*c(5);
cprime(3) =  k(1)*c(1)*c(2) - k(2)*c(2)*c(3) + k(3)*c(1)*c(3)...
    - 2*k(4)*c(3).^2;
cprime(4) = k(2)*c(2)*c(3) + k(4)*c(3).^2;
cprime(5) = k(3)*c(1)*c(3) - k(5)*c(5);

fcount = fcount + 1;
end