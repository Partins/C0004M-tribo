%% Part 2 The analytical function for discribing pressure, p(x) =f1*f2

%SETUP

N = 10;
L = 0.02;
%x = 1/N;    %Dimensionlös! x/L

f1 = @(mu,U,L,hmin) 6*mu*U*L/(hmin^2);
f2 = @(k,x) 1/k*(1/(1 + k +x*k) -((1+k)/(2+k))*(1/(1+k+k*x)^2)-(1/(2+k)));

px = f1(mu,U,L,hmin)*f2(k,x);

