% The discrete representation of the analytical pressure and load carrying
% capacity for numerical evaluation. Code is taken from the presentation
% material. 

visc    = 10;                            % Viscosity [Pas] 
U       = 5;                            % Runner speed [m/s]
L       = 10;                            % Bearing length [m]
hmin    = 0.0001;                            % Trailing edge film thickness [m]
k       = linspace(0,1e-3,100);         % Slope parameter 
x       = 1;                            % Position

tryckUppg = ((6.*visc.*U.*L)./(hmin.^2)); % Tryckuppgyggnadens variation med alla designparametrar utom k.

p = @(x) tryckUppg.*((1./k).*(1./(1+k-k.*x./L))-((1+k)./(2+k)).*(1./(1+k-k.*x./k).^2)-1./(2+k));
N = log(1+k)./k.^2 - 2./(2+k)./k;

[maxN, indmax] = max(N);
disp([N(indmax),k(indmax)]);
plot(k,p(2))

