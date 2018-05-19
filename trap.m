%% Trapets integration of p(x)

function fa = trap(px,l,u,N)  %trycket i punkterna f�r vektor x,
                              %nedre v�rde av x, �vre v�rde f�r x, antal
                              %punkter att approximera �ver.
h = (u-l)/N;
fa = h*0.5*sum([px(1) 2*px(2:end-1) px(end)]);

end