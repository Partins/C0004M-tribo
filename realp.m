%% Evaluate p(x) analytically given a specifik k and a range of x.

function px = realp(mu,U,L,hmin,x,k)

X = x/L; %vektor

f1 = (6*mu*U*L)/(hmin^2); %konstant

px = f1.*(1./k.*(1./(1 + k - k.*X)-((1+k)/(2+k)).*(1./((1 + k - k.*X).^2))-1/(2+k)));

end