%% Increasing x, while keeping konstant k.
% A function for finding the pressure over the length of the bearing with a
% konstant k.

%Setup
N       = 1000;                         % Number of unknowns
mu      = 0.01;                         % Viscosity [Pas] 
U       = 1;                            % Runner speed [m/s]
l       = 0.1;                          % Bearing length [m]
hmin    = 10^-6;                        % Trailing edge film thickness [m]
k       = 4;                         % Slope parameter [No dimension].
p0      = 100000;                       % left limit [Pa]
pL      = 100000;                       % right limit [Pa]
xi      = l/(N);                        % step length
x       = [xi:xi:l]                     % x axis

h = hmin*(1+k-(k/l).*x);

a = (h.^3)/(12*mu);
F = (U/2).*h;

%% Diagonal vectors

c1 = (a(1:N-1)+a(2:N))./(2.*x(1:N-1).^2); %lower diagonal
d1 = -(a(1:N-1)+2*a(2:N)+a(1:N-1))./(2.*x(1:N-1).^2); %main diagonal
e1 = (a(1:N-1)+a(2:N))./(2.*x(1:N-1).^2); %top diagonal

%% Right hand side
z = (F(2:N)-F(1:N-1))./(2.*x(1:N-1));

%% Setup vectors and sparse matrix A

e = [0 e1(1:end-1)];
c = [c1(2:end) 0];
A = spdiags([c' d1' e'],[-1 0 1],N-1,N-1);

f = [z(1)-c1(1)*p0 z(2:end-1) z(N-1)-e(N-1)*pL]';

%% Solving the linear system of equations.
u = A\f;

%% Plotting p over x.
x = [0 x]
p = [p0 u' pL];
plot(x,p);

%% Finding max values

[maxp,xval] = max(p);
xloc = xval/N;
