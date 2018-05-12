visc    = 10;                           % Viscosity [Pas] 
U       = 5;                            % Runner speed [m/s]
L       = 10;                           % Bearing length [m]
hmin    = 0.0001;                       % Trailing edge film thickness [m]
k       = linspace(0,1e-3,100);         % Slope parameter 
x1       = 1;                            % Position
x2      = 0:0.1:10;
y0      = 0;

tryckUppg = ((6.*visc.*U.*L)./(hmin.^2));
p = @(x) tryckUppg.*((1./k).*(1./(1+k-k.*x./L))-((1+k)./(2+k)).*(1./(1+k-k.*x./k).^2)-1./(2+k));
y1 = p(x1);

func = @(x,h) (6*visc*U)/(hmin(1+k-(k/L)*x));
y2 = rungekutta4(func,x2,y0);

plot(x1,y1,x2,y2);