%% Increasing x, while keeping konstant k.
% A function for finding the pressure over the length of the bearing with a
% konstant k.

%Setup
N       = 1000;                         % Number of unknowns
mu      = 0.01;                         % Viscosity [Pas] 
U       = 1;                            % Runner speed [m/s]
l       = 0.1;                          % Bearing length [m]
hmin    = 10^-6;                        % Trailing edge film thickness [m]
k       = 1.18;                         % Slope parameter [No dimension].
p0      = 100000;                       % left limit [Pa]
pL      = 100000;                       % right limit [Pa]
xi      = l/(N);                        % step length
x       = [xi:xi:l]                     % x axis

h = hmin*(1+k-(k/l).*x);

a = (h.^3)/(12*mu); %L�ngd N
F = (U/2).*h;

%% Diagonaler matris A

c1 = (a(1:N-1)+a(2:N))./(2.*x(1:N-1).^2); %underdiag
d1 = -(a(1:N-1)+2*a(2:N)+a(1:N-1))./(2.*x(1:N-1).^2); %huvudiag
e1 = (a(1:N-1)+a(2:N))./(2.*x(1:N-1).^2); %�verdiag

%% H�gerled
z = (F(2:N)-F(1:N-1))./(2.*x(1:N-1));

%% Upps�ttning matris A

e = [0 e1(1:end-1)];
c = [c1(2:end) 0];
A = spdiags([c' d1' e'],[-1 0 1],N-1,N-1); %Se ekv.31, gles matris A

f = [z(1)-c1(1)*p0 z(2:end-1) z(N-1)-e(N-1)*pL]'; %Se ekv.32, h�gerled 

%% L�sning av ekvationssystem
u = A\f;

%% Plot av tryck med varierande x, k konstant.
x = [0 x]
p = [p0 u' pL];
plot(x,p);

[maxp,xval] = max(p);
xloc = xval/N;
