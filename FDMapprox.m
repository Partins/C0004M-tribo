%% This scrip and its related functions realp(),trap() and force() use
% well known numerical methods for approximating the Load Carrying Capacaty
% over a 1D padtilt-bearing.
cc;

%%Setup
N       = 100;                          % Number of unknowns
mu      = 0.01;                         % Viscosity [Pas] 
U       = 1;                            % Runner speed [m/s]
l       = 0.1;                          % Bearing length [m]
hmin    = 10^-6;                        % Trailing edge film thickness [m]
p0      = 0;                            % left limit [Pa]
pL      = 0;                            % right limit [Pa]
xi      = l/(N);                        % step length
x       = [0:xi:l];                     % x axis
k       = linspace(0,4,N);              % Slope parameter from 0<=k<=4


%% LOOP FOR VALUES OF FORCE(N) ANALYTICALLY AND NUMERICALLY
% In the following loop the script will approximate the pressure p over x
% given a k where k = k(i) 
% The integral of p(x) is then approximated using the trapezoidal rule
% using function trap.m

for i = 1:N
    
% Creating a linear system of equations given k(i) and solving for p.
h = hmin.*(1+k(i)-(k(i)./l).*x);
F = (U/2).*h;
a = (h.^3)/(12*mu);

% Diagonal vectors

c1 = (a(1:N-1)+a(2:N))./(2*xi^2); %lower diagonal
d1 = -(a(1:N-1)+2*a(2:N)+a(3:N+1))./(2*xi^2); %main diagonal
e1 = (a(2:N)+a(3:N+1))./(2*xi^2); %top diagonal

z = (F(3:N+1)-F(1:N-1))./(2.*xi);

% Setup vectors and sparse matrix A

e = [0 e1(1:end-1)];
c = [c1(2:end) 0];
A = spdiags([c' d1' e'],[-1 0 1],N-1,N-1);

f = [z(1)-c1(1)*p0 z(2:end-1) z(N-1)-e(N-1)*pL]';

%  Solving the linear system of equations.
p = A\f;

%pressure including boundary conditions.
p = [p0 p' pL];

%% Controll p(x) analytically

%Call function returning analyticall value for p(x)
    if i ==1
        px = zeros(length(x),1);
    else
        px = realp(mu,U,l,hmin,x,k(i));
    end
    
    
%% Numerical approximation: Integrating p(x) over x for the load carrying capacity given a specifik k. using trapezoid rule.
    if i == 1
        fa = trap(p,0,l,N);
    else
        fa = [fa trap(p,0,l,N)];
    end
%% The analytical integratoin for N(k)
    if i==1   
        fx = 0;
    else
        fx = [fx force(mu,U,l,hmin,k(i))];
    end

end
fdx = fx/((6*mu*U*l^2)/(hmin^2)); %Dimensionslös
fda = fa/((6*mu*U*l^2)/(hmin^2)); %Dimensionslös
plot(k,fda,k,fdx,'.')
legend('LCC numerical','LCC analytical')

%%BASELINE BDUMTSS

%Dimensionslös k från Matlabscriptet.

kb = linspace(0.01,4,10001);
LCC = log(1+k)./k.^2 - 2./(2+k)./k;
[maxN,indmax] = max(LCC);
disp([LCC(indmax),k(indmax)]);
figure
plot(k,fda,k,LCC,'.');
legend('LCC numerical2','LCC analytical2')