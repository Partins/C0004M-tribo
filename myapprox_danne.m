

%Setup
N       = 10;                          % Antal obekanta/Antal delar p� intervallet.
mu      = 0.01;                         % Viscosity [Pas] 
U       = 1;                            % Runner speed [m/s]
l       = 0.02;                         % Bearing length [m]
hmin    = 10^-6;                        % Trailing edge film thickness [m]
k       = 10^-14;                         % Slope parameter
p0      = 0;                            % V�nster Randvillkor
pL      = 0;                            % H�ger randvillkor
xi      = l/(N);                        % stegl�ngd
x       = [xi:xi:l]                     % x axel

h = hmin*(1+k-(k/l).*x);                % v�rden f�r h med l�ngd N


a = (h.^3)/(12*mu); %L�ngd N
F = (U/2).*h;
% X = [0 x:x:l-x]

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

p = [p0 u' pL];
plot([0 x],p);

%% Hitta f�r vilket k �r st�rst.

%% J�mf�r mot den analytiska l�sningen.

%pa = validationfunc(N,p0,pL,mu,U,L,hmin,k,x)

%% ALTERNATIV KOD

%Samma som ovanst�ende men med trycket varierande efter
%lutningskoefficienten.

%Innan �ndringar f�r x ist�llet f�r k.

%Setup
N       = 10;                          % Antal obekanta/Antal delar p� intervallet.
mu      = 0.01;                         % Viscosity [Pas] 
U       = 1;                            % Runner speed [m/s]
l       = 0.02;                         % Bearing length [m]
hmin    = 10^-6;                        % Trailing edge film thickness [m]
k       = linspace(0,1e-3,N+1);         % Slope parameter
p0      = 0;                            % V�nster Randvillkor
pL      = 0;                            % H�ger randvillkor

h = hmin*(1+k-(k/l));                % v�rden f�r h med l�ngd N
xs = l/N

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
plot(k,p);
