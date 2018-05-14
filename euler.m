function ye = euler(func,x,y0)
    ye = zeros(length(x));
    h = diff(x);
    ye(1) = y0;
    
    for k = 1:length(x)-1
        ye(k+1) = ye(k) + h(k).*func(x(k),ye(k));
    end
end