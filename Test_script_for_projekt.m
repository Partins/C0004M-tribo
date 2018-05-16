visc    = 10;                           % Viscosity [Pas] 
U       = 5;                            % Runner speed [m/s]
L       = 10;                           % Bearing length [m]
hmin    = 0.0001;                       % Trailing edge film thickness [m]
k       = linspace(0,1e-3,100);         % Slope parameter 
%k        = 50;
N       = 10;
u(1)    = 0;

h = hmin.*(1+k-(k./L));
a = (h.^3)./(12.*visc);
F = (U/2).*h;
deltax = L/N;
    
    for i = 2:N-1
        x(i) = i.*L./N;

        c(i) = (a(i-1)+a(i))./(2.*deltax.^2);
        d(i) = -(a(i-1)+2.*a(i)+a(i+1))/(2.*deltax.^2);
        e(i) = (a(i)+a(i+1))/(2.*deltax.^2);
        z(i) = (F(i+1)-F(i-1))./(2.*deltax);
        g(i) = -(F(i)/deltax);
        K(i) = F(i-1)/deltax;
    end
A = spdiags([c' d' e'],[-1 0 1],N-1,N-1);
B = spdiags([g' K'],[-1 0],N-1,N-1);
f = [z(1)-c(1)*u(1) z(2:end-1) z(N-1)-e(N-1)]';