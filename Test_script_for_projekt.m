%Setup

visc    = 10;                           % Viscosity [Pas] 
U       = 5;                            % Runner speed [m/s]
L       = 10;                           % Bearing length [m]
hmin    = 0.0001;                       % Trailing edge film thickness [m]
k       = linspace(0,1e-3,100);         % Slope parameter 
N       = 10;                           % Precision
p(1)    = 0;
p(N)    = 0;
h = hmin.*(1+k-(k./L));

a = (h.^3)./(12.*visc);
F = (U/2).*h;
deltax = L/N;
    
    for i = 2:N-1
        c1(i) = (a(i-1)+a(i))./(2.*deltax.^2);
        d(i) = -(a(i-1)+2.*a(i)+a(i+1))/(2.*deltax.^2);
        e1(i) = (a(i)+a(i+1))/(2.*deltax.^2);
        z(i) = (F(i+1)-F(i-1))./(2.*deltax);
        g(i) = -(F(i)/deltax);
        K(i) = F(i-1)/deltax;
    end

e = [0 e1(1:end-1)];
c = [c1(2:end) 0];
A = spdiags([c' d' e'],[-1 0 1],N-1,N-1); %Se ekv.31
B = spdiags([g' K'],[-1 0],N-1,N-1); %Tror inte vi beh?ver den h?r
f = [z(1)-c(1)*p(1) z(2:end-1) z(N-1)-e(N-1)*p(N)]'; %Se ekv.32

u = A\f;
p = [p(1) u' p(N)];
plot(0:N,p);