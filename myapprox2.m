%% Varible slope parameter k.

%Setup
N       = 512;                           % Antal obekanta/Antal delar på intervallet.
mu      = 0.1;                          % Viscosity [Pas] 
U       = 1;                            % Runner speed [m/s]
l       = 0.02;                         % Bearing length [m]
hmin    = 10^-6;                        % Trailing edge film thickness [m]
k       = linspace(0,4,N+1);            % Slope parameter 0 till 1
p0      = 100000;                       % left limit [Pa]
pL      = 100000;                       % right limit [Pa]
xs      = l/N                           % steplength
h       = hmin.*(1+k-(k./l).*xs);       % vektor h, with length N+1


a = (h.^3)./(12.*mu); %Längd N+1
F = -(U/2).*h;

%% Diagonaler matris A

c1 = (a(1:N-1)+a(2:N))/(2*xs^2); %underdiag
d1 = -(a(1:N-1)+2*a(2:N)+a(3:N+1))/(2*xs^2); %huvudiag
e1 = (a(2:N)+a(3:N+1))/(2.*xs.^2); %överdiag

%% Högerled
z = (F(3:N+1)-F(1:N-1))/(2*xs); %kontrollera 

%% Uppsättning matris A
c = [c1(2:end) 0]; %övre
e = [0 e1(1:end-1)]; %undre
A = spdiags([c' d1' e'],[-1 0 1],N-1,N-1); %Se ekv.31, gles matris A

f = [z(1)-c1(1)*p0 z(2:end-1) z(N-1)-e(N-1)*pL]'; %Se ekv.32, högerled 

%% Lösning av ekvationssystem
u = A\f;

%% Plot av tryck med varierande k.

p = [p0 u' pL];
figure
plot(k,p);
[maxval,maxindex] = max(p);
xmax = maxindex/N
