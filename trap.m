%% Trapets integration of p(x)

function fa = trap(px,l,u,N)  %trycket i punkterna för vektor x,
                              %nedre värde av x, övre värde för x, antal
                              %punkter att approximera över.
h = (u-l)/N;
fa = h*0.5*sum([px(1) 2*px(2:end-1) px(end)]);

end