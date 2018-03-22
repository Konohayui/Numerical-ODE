function oder = order(A, b, c)
  oder = zeros(4,4);
  for n = 1:4
  o1 = (c.').^(n-1)*b;
  oder(1,n) = o1;
  end
  o23 = b.'*A*c;
  o34 = b.'*A*c.^2;
  o24 = (b.*c).'*A*c;
  o4 = b.'*A^2*c;
  oder(2,3) = o23;
  oder(2,4) = o24;
  oder(3,4) = o34;
  oder(4,4) = o4;
end