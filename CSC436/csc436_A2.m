format compact

%Q2
xv = linspace(a,b,1000)'; 
e = zeros(5,1);
enodes = e;
ecube = e;
el = e;
u
fprintf( '       qurtic     cubic      linear        order\n')

for i = 3:7
    n = pow2(i);
    x = linspace(a,b,n+1)';
    for kk=1:length(x)
        x(i)=x(i)^2;
    end
    xmid = (x(1:n) + x(2:n+1))/2;
    y = u(x);
    yqrt = u([a;xmid;b])';%special input for qrtspline
    yv = qrtspline(x,yqrt,xv);
    yvcube = spline(x,y,xv);
    yvl = interp1(x,y,xv);
    y
    enodes(i-2) = max(abs(qrtspline(x,yqrt,x)-y)); e(i-2) = max(abs(yv-u(xv)));
    ecube(i-2) = max(abs(yvcube-u(xv)));
    el(i-2) = max(abs(yvl-u(xv)));
    if i == 3
        fprintf( '%4d %10.2e %10.2e %10.2e\n', ...
        [n, e(1),ecube(1),el(1)]);
    else
        fprintf( '%4d %10.2e %10.2e %10.2e %10.2f\n', ...
        [n, e(i-2), ecube(i-2),el(i-2),log2(e(i-3)/e(i-2))]);
    end
end
n = 2.^(3:7); 
loglog(n,enodes,'-',n,e,'--',n,ecube,'-.',n,el,':') 
legend('quartic error at nodes',...
'quartic error','cubic spline error','linear spline error')


% x Vector, y Vector of u(x), xv: Evalustion vactor
function yv = qrtspline(x,y,xv)
    n = length(x)-1;
    h = (x(n+1)-x(1))/n;
    % Construct sparse matrix
    L = [ones(1,n+4)' repmat(76,n+4,1),repmat(230,n+4,1)...
    repmat(76,n+4,1) ones(1,n+4)'];
    l1 = [1,-1,-1,1];
    l2 = [1,11,11,1];
    H =spdiags(L,[-2 -1 0 1 2],n+4,n+4);
    H(1,1:4) = l1;
    H(2,1:4) = l2;
    H(end,end-3:end) = l1;
    H(end-1,end-3:end) = l2;
    H = (1/384).*H;
    k = [3040 -5859 4809 -2835 999 -154]/378;
    u02nd = k * y(1:6)/h^2;
    un2nd = k(6:-1:1) * y(n-3:n+2)/h^2;
    u = [u02nd * h^2/192; y; un2nd * h^2/192];
    u(2) = u(2)/16;
    u(n+3) = u(n+3)/16;
    c = H\u;
    B = zeros(length(xv),n+4);
        for i = -1:n+2
           B(:,i+2) = funQ1( (xv-x(1))/h -i +3 );
        end
    yv = B*c;
    yv
end


function y = phi(x)
n = length(x);
if n >= 2
    y = zeros(0,1);
end
for i = 1:n
    if (x(i) >= 0) && (x(i)<=1)
        y(i) = x(i)^4;
    elseif (x(i)>1) && (x(i)<=2)
        y(i) = x(i)^4-5*(x(i)-1)^4;
    elseif (x(i)>2) && (x(i)<=3)
        y(i) = x(i)^4 - 5*(x(i)-1)^4 + 10 * (x(i)-2)^4;
    elseif (x(i)>3) && (x(i)<=4)
        y(i) = x(i)^4 - 5*(x(i)-1)^4 + 10 * (x(i)-2)^4 - 10*(x(i)-3)^4;
    elseif (x(i)>4) && (x(i)<=5)
        y(i) = x(i)^4 - 5*(x(i)-1)^4 + 10 * (x(i)-2)^4 - 10*(x(i)-3)^4 + 5*(x(i)-4)^4;
    else
       y(i) = 0;
    end
end
    y = y/24;
end
