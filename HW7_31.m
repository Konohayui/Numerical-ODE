theta = 0:0.0001:2*pi;
z = exp(1i*theta);

rho = z.^3 + z.^2/4 - z/2 - 3/4; 
phi = 1/8*(19*z.^2 + 5);         
R = rho./phi;                    

plot(real(R), imag(R));
axis square

A = [200 398 198; -500 -696 -296; 500 694 294];
e = eig(A);
			
% obtain angles for each eigenvalue
% and all approximation
ang = angle(R);                   
an1 = angle(e(1));
an2 = angle(e(2));
an3 = angle(e(3));           
			
% find the angle for each point 
% that have the same angle with eigenvalues
ang_ind1 = find(abs(ang - an1) < 0.001);
approx1 = R(ang_ind1);
ang_ind2 = find(abs(ang - an2) < 0.001);
approx2 = R(ang_ind2);
ang_ind3 = find(abs(ang - an3) < 0.001);
approx3 = R(ang_ind3);
			
% get h_hat for each approximation
h1 = abs(approx1)./abs(e(1));
h2 = abs(approx2)./abs(e(2));
h3 = abs(approx3)./abs(e(3));
			
% obtain the minimum h_hat such
% that h_hat = h*lambda
minh = min([h1 h2 h3]);