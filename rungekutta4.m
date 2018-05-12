function yrk = rungekutta4(func,x,y0)
    yrk = zeros(length(x));
    h = diff(x);
    yrk(1) = y0;
    
    for k = 1:length(x)-1
        s1 = func(x(k),yrk(k));
        s2 = func(x(k)+0.5*h(k),yrk(k)+0.5*h(k)*s1);
        s3 = func(x(k)+0.5*h(k),yrk(k)+0.5*h(k)*s2);
        s4 = func(x(k)+0.5*h(k),yrk(k)+0.5*h(k)*s3);
        yrk(k+1) = yrk(k) + (1/6)*h(k)*(s1+2*s2+2*s3+s4);
    end
end