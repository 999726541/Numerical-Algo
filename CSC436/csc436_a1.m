format compact

v = linspace(0,1,1000);
u = exp(v);


factors = [4,8,16,32,64,128];
err_p = [0,0,0,0,0,0];
err_l = [0,0,0,0,0,0];
err_h = [0,0,0,0,0,0];
err_c = [0,0,0,0,0,0];
for i = 2:7
    n = 2^i;
    x = linspace(0, 1, n+1);
    y = exp(x);
    p = polyval(polyfit(x,y,n),v);
    err_p(i-1) = max(abs(p - u));
    l = interp1(x,y,v,'linear');
    err_l(i-1) = max(abs(l - u));
    c = interp1(x,y,v,'spline');
    err_c(i-1) = max(abs(c - u));
    h = interp1(x,y,v,'pchip');
    err_h(i-1) = max(abs(h - u));
    if i ==2
        fprintf('%4d %11.4e %11.4e %11.4e %11.4e\n',[n;err_p(i-1);err_l(i-1);err_c(i-1);err_h(i-1)]);
    else
        fprintf('%4d %11.4e %11.4e %11.4e %11.4e %7.2f %5.2f %5.2f %5.2f\n',[n;err_p(i-1);err_l(i-1);err_c(i-1);err_h(i-1);log2(err_p(i-2)/err_p(i-1));log2(err_l(i-2)/err_l(i-1));log2(err_c(i-2)/err_c(i-1));log2(err_h(i-2)/err_h(i-1))]);
    end
    if i == 2
        plot(v,u,v,p,'--');hold on;
    end
    if i == 7
        plot(v,p,':')
        hold off
        figure
        axis([0 1 0.9 exp(1)+0.2])
        plot(log([2,3,4,5,6,7]),log(err_p),'-');hold on;
        plot(log([2,3,4,5,6,7]),log(err_l),'--');
        plot(log([2,3,4,5,6,7]),log(err_c),'-.');
        plot(log([2,3,4,5,6,7]),log(err_h),':');
    end
end


