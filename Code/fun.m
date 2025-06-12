function y = fun(x)
    t1 = x(1:end-1);
    t2 = x(2:end);
    y = sum(100*(t2-t1.^2).^2+(t1-1).^2);
end
    