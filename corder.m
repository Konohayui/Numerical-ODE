function [q, L]= corder(alpha, beta, k)
% set k be a boundary where q<=k
c0 = sum(alpha);
r = length(alpha) - 1;
j = [0:r];

if c0 == 0
    L = zeros(k + 1,1);
    L(1) = c0;
    for l = 1:k
        ck = j.^l*alpha - l*j.^(l - 1)*beta;
        L(l+1) = sum(ck);
    end
    q = sum(L == 0) - 1;
end

end
