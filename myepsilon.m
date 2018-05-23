%% This function calculates the realtive error between the analytical and numerical values for the LCC.

function epsilon = myepsilon(N) %N = Number of elements N used for the approximation.


%% Setup

% Variables vareing over N.
l       = 0.1;                          % Bearing length [m]
xi      = l/(N);                        % step length
x       = [0:xi:l];                     % x axis
k       = linspace(0,4,N);              % Slope parameter where 0<=k<=4

% Constants
mu      = 0.01;                         % Viscosity [Pas] 
U       = 1;                            % Runner speed [m/s]
hmin    = 10^-6;                        % Trailing edge film thickness [m]
p0      = 0;                            % left limit [Pa]
pL      = 0;                            % right limit [Pa]

%% LOOP FOR VALUES OF FORCE(N) ANALYTICALLY AND NUMERICALLY
% In the following loop the function will approximate the pressure p over x
% given a k where k = k(i) 
% The integral of p(x) is then approximated using the trapezoidal rule
% using the function trap.m

for i = 1:N
    
%% Creating a linear system of equations given k(i) and solving for p.

% Setup.
h = hmin.*(1+k(i)-(k(i)./l).*x);
F = (U/2).*h;
a = (h.^3)/(12*mu);

% Diagonal vectors

c1 = (a(1:N-1)+a(2:N))./(2*xi^2); %lower diagonal
d1 = -(a(1:N-1)+2*a(2:N)+a(3:N+1))./(2*xi^2); %main diagonal
e1 = (a(2:N)+a(3:N+1))./(2*xi^2); %upper diagonal

z = (F(3:N+1)-F(1:N-1))./(2.*xi);

% Setup vectors and sparse matrix A

e = [0 e1(1:end-1)];
c = [c1(2:end) 0];
A = spdiags([c' d1' e'],[-1 0 1],N-1,N-1);

f = [z(1)-c1(1)*p0 z(2:end-1) z(N-1)-e(N-1)*pL]';

% Solving the linear system of equations.
p = A\f;

% pressure including boundary conditions.
p = [p0 p' pL];

%% Numerical approximation for Integrating p(x)

if i == 1
    fa = trapz(x,p);
else
    fa = [fa trapz(x,p)];
end

end

%Compute the realtive error.
fx = log(1+k)./k.^2 - 2./(2+k)./k;
fx(1) = p0;
f1 = (6*mu*U*l^2)/(hmin^2);
fa = fa/f1; %Remove dimension.

epsilon = (fa-fx)/fx;

end