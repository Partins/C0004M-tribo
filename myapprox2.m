%% Varible slope parameter k, konstant xs (length of "xi")

%Setup
N       = 10;                           % Antal obekanta/Antal delar p� intervallet.
mu      = 0.01;                         % Viscosity [Pas] 
U       = 1;                            % Runner speed [m/s]
l       = 0.02;                         % Bearing length [m]
hmin    = 10^-6;                        % Trailing edge film thickness [m]
k       = linspace(0,1e-3,N);           % Slope parameter
p0      = 0;                            % left limit
pL      = 0;                            % right limit
xs      = l/N                         % steplength
h       = hmin*(1+k-(k/l));             % vektor h, with length N+1


a = (h.^3)/(12*mu); %L�ngd N
F = (U/2).*h;
% X = [0 x:x:l-x]

%% Diagonaler matris A

c1 = (a(1:N-1)+a(2:N))/(2*xs^2); %underdiag
d1 = -(a(1:N-1)+2*a(2:N)+a(1:N-1))/(2*xs^2); %huvudiag
e1 = (a(1:N-1)+a(2:N))/(2.*xs.^2); %�verdiag

%% H�gerled
z = (F(2:N)-F(1:N-1))/(2*xs);

%% Upps�ttning matris A

e = [0 e1(1:end-1)];
c = [c1(2:end) 0];
A = spdiags([c' d1' e'],[-1 0 1],N-1,N-1); %Se ekv.31, gles matris A

f = [z(1)-c1(1)*p0 z(2:end-1) z(N-1)-e(N-1)*pL]'; %Se ekv.32, h�gerled 

%% L�sning av ekvationssystem
u = A\f;

%% Plot av tryck med varierande k.

p = [p0 u' pL];
figure
plot([0 k],p);