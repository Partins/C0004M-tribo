function fx = force(mu,U,L,hmin,k)

f1 = (6*mu*U*L^2)/(hmin^2);

fx = f1.*((1./(k.^2)).*log(1+k)-((1./k).*(2./(2+k))));

end