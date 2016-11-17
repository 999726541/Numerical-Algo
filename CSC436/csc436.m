v = linspace(a,b,1000);
u = f(v);
axis([0 1 -0.1 1.2]);
hold on
ep = zeros(6,1);
el = ep;
eh = ep;
ec = ep;
plot(v,u,'-')
for i = 2:7
    n = pow2(i);
    x = linspace(a,b, n + 1);
    y = f(x);
    p = polyval(polyfit(x,y,n),v);
    ep(i-1) = max(abs(p - u));
    l = interp1(x,y,v,'linear');
    el(i-1) = max(abs(l - u));
    c = interp1(x,y,v,'spline');
    ec(i-1) = max(abs(c - u));
    h = interp1(x,y,v,'pchip');
    eh(i-1) = max(abs(h - u));
    if i == 2
        fprintf('%4d %11.4e %11.4e %11.4e %11.4e\n',[n;ep(i-1);el(i-1);ec(i-1);eh(i-1)]);
    else
        fprintf('%4d %11.4e %11.4e %11.4e %11.4e %7.2f %5.2f %5.2f %5.2f\n',[n;ep(i-1);el(i-1);ec(i-1);eh(i-1);log(ep(i-2)/ep(i-1));log(el(i-2)/el(i-1));log(ec(i-2)/ec(i-1));log(eh(i-2)/eh(i-1))]);
    end
    if i == 5
        plot(v,p,'--')
    end
    if i == 7
        plot(v,p,':')
        hold off
    end
end

