visc    = 10;                           % Viscosity [Pas] 
U       = 5;                            % Runner speed [m/s]
L       = 10;                           % Bearing length [m]
hmin    = 0.0001;                       % Trailing edge film thickness [m]
k       = linspace(0,1e-3,100);         % Slope parameter 
N       = ???;

x(i) = i.*L./N;
h = hmin(1+k-(k./L));
a = h.^3/(12.*visc);
F = U/2;

c(i) = (a(i-1)+a(i))/(2.*deltax.^2);
d(i) = -(a(i-1)+2.*a(i)+a(i+1))/(2.*deltax.^2);
e(i) = (a(i)+a(i+i))/(2.*deltax.^2);
z(i) = (F(i+1)-F(i-1))/(2.*deltax);
g(i) = -(F(i)/deltax);
k(i) = F(i-1)/deltax;
