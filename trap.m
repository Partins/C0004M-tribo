%% Trapets integration of p(x)

function fa = trap(px,l,u,N)
h = (u-l)/N;
fa = h*0.5*sum([px(1) 2*px(2:end-1) px(end)]);

end