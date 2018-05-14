function yr = ralston(func,x,y0)
    yr = zeros(length(x));
    h = diff(x);
    yr(1) = y0;
    
    for k = 1:length(x)-1
        s1 = func(x(k),yr(k));
        s2 = func(x(k)+0.75*h(k),yr(k)+0.75*h(k)*s1);
        yr(k+1) = yr(k) + h(k)*((1/3)*s1+(2/3)*s2);
    end
end