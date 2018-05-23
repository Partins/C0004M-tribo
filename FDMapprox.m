%% A FDM based solver for the Load carrying capacity of a 1D tilted pad bearing. 
clear all; clc;

%% Setup

% Variables varying over N.
N       = 10;                          % Number of unknowns, "resolution"
l       = 0.1;                          % Bearing length [m]
xi      = l/(N);                        % step length
x       = [0:xi:l];                     % x axis
k       = linspace(0,4,N);              % Slope parameter where 0<=k<=4

% Setup for Convergence study

NL = 10;
NU = 50;

% Constants
mu      = 0.01;                         % Viscosity [Pas] 
U       = 1;                            % Runner speed [m/s]
hmin    = 10^-6;                        % Trailing edge film thickness [m]
p0      = 0;                            % left limit [Pa]
pL      = 0;                            % right limit [Pa]

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

z = (F(3:N+1)-F(1:N-1))./(2.*xi); %right hand side

% Setup vectors and sparse matrix A

e = [0 e1(1:end-1)];
c = [c1(2:end) 0];
A = spdiags([c' d1' e'],[-1 0 1],N-1,N-1);

f = [z(1)-c1(1)*p0 z(2:end-1) z(N-1)-e(N-1)*pL]'; %right hand side

% Solving the linear system of equations.
p = A\f;

% pressure including boundary conditions.
p = [p0 p' pL];

%% Controll p(x) analytically

%Call function returning analyticall value for p(x)
    if i == 1
        px = zeros(length(x),1);
    else
        px = analyticp(mu,U,l,hmin,x,k(i));
    end
%     
%     
%% Numerical approximation for Integrating p(x)

% Numerical approximation over x for the load carrying capacity given a 
% specifik k. using the trapezoid rule.

    if i == 1
        fa = trap(p,0,l,N);
    else
        fa = [fa trap(p,0,l,N)]; 
    end

end
%% The analytical integration for N(k) without dimension.
fx = log(1+k)./k.^2 - 2./(2+k)./k;
fx(1) = p0;

%% Plotting the approximated pressure and the analyticall pressure and plotting the approximated Load Carrying Capacity against the analytical result
%Removing dimensions from the approximation fa.

f1 = (6*mu*U*l^2)/(hmin^2);
fa = fa/f1;

%Plotting the approximated pressure p against the analyticall pressure px
%over x for k = K(N) = 4.
figure('Name','Pressure over 0<=x<=l','NumberTitle','off');
plot(x,p,x,px,'.')
legend('p(x) numerical','p(x) analytical')
xlabel('$x$ [m]','interpreter','latex','fontsize',16);
ylabel('$p(x) [Pa], k=4$','interpreter','latex','fontsize',16);


%Plotting the results for the LCC

figure('Name','LCC','NumberTitle','off');
plot(k,fa,k,fx,'.')
legend('LCC numerical','LCC analytical')
xlabel('$k$','interpreter','latex','fontsize',16);
ylabel('$\overline{N}$','interpreter','latex','fontsize',16);

%% Convergence study for the relative error

epsilon = ones(1,NU-NL);
for t = 1:NU-NL+1
   epsilon(t) = myepsilon(t+NL);
end
figure('Name','Relative error','NumberTitle','off');
plot(NL:NU,abs(epsilon))
xlabel('$N$','interpreter','latex','fontsize',18);
ylabel('$\epsilon$','interpreter','latex','fontsize',24);