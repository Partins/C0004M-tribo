function yh = heun(func,x,y0)
    yh = zeros(length(x));
    h = diff(x);
    yh(1) = y0;
    
    for k = 1:length(x)-1
        s1 = func(x(k),yh(k));
        s2 = func(x(k+1),yh(k)+h(k)*s1);
        yh(k+1) = yh(k) + h(k)*0.5*(s1+s2);
    end
end