function uprime = myode2(t,u)
uprime = zeros(size(u));
if t == 0 && u(1) == 0
    uprime(1) = 0;
else
    uprime(1) = 2*t*u(1)/(t^2 + u^2);
end
end 