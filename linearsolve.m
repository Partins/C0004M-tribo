%Numerical approximation of reynolds equation by solving the linear system
%for p.
 %Initialisera Variabler
 
 N = 128;               %Antal bitar
 y0 = 0;   %Randvillkor
 yL = 0;   %Randvillkor
 l = 1;    % Lagrets "l�ngd"
 x = linspace(0,l,N); %intervall
 hmin = 4e-3; %minsta filmtjocklek
 k = 1; %Lutningskoefficient
 mu = 10; %Viskostitet
 i = 1.
 
 h = hmin.*(1 + k - k.*x./l)
 
 a = h.^3./12.*mu;
 
 ah(1) = (a(1) + a(1))./2;
 ah(N) = (a(N) + a(N))./2;
 
 for i = 2:length(a)-1
 
 ah(i) = (a(i) + a(i))./2; %ah = a halva
 
 end
 