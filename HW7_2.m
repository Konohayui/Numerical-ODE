[x,y] = meshgrid(-6:0.01:6, -6:0.01:6);
z = x + 1i*y;

R = 1 + z + z.^2/2 + z.^3/6 + z.^4/24;% Taylor polynomial
zlevel = abs(R);

figure
contour(x, y, zlevel,[1 1]);
axis square