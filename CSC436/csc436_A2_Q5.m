
for m = 7:2:13
    f = @(x)x.^(m/2);
    ER=2/(m+2);
    fprintf('Function of x^(m/2) @ m= %4d\n',m)
    fprintf('panel#   Simpson    twopntgauss    lobattoGuass\n')
    for i = 0:5
        n = 2^i;
        fprintf('%4d %10.2e %10.2e %10.2e\n',[n,ER-simpson(f,0,1,n),ER-twopntgauss(f,0,1,n),...
ER-gausslobatto(f,0,1,n)]);

    end

end


    f = @(x)max(x-1/3,0);
    ER=2/9;
    fprintf('Function of Max(x-1/3,0)\n')
    fprintf('panel#   Simpson    twopntgauss    lobattoGuass\n')
    for i = 0:5
        n = 2^i;
        fprintf('%4d %10.2e %10.2e %10.2e\n',[n,ER-simpson(f,0,1,n),ER-twopntgauss(f,0,1,n),...
ER-gausslobatto(f,0,1,n)]);

    end



    f = @(x)ff(x);
    ER=2/3;
    fprintf('Function of decision function\n')
    fprintf('panel#   Simpson    twopntgauss    lobattoGuass\n')
    for i = 0:5
        n = 2^i;
        fprintf('%4d %10.2e %10.2e %10.2e\n',[n,ER-simpson(f,0,1,n),ER-twopntgauss(f,0,1,n),...
ER-gausslobatto(f,0,1,n)]);

    end



function I = simpson(f,a,b,n)
x = linspace (a,b,2*n+1);
if n == 1
    I = (b-a)/6*(f(a)+4*f((a+b)/2)+f(b));
else
    I = (b-a)/6/n * (f(a) + f(b) + 4*sum(f(x(2:2:2*n+1)))+2*sum(f(x(3:2:2*n-1))));
end
end

function I = gausslobatto(f,a,b,n)
h = (b-a)/n;
x1 = (a+h/2-sqrt(5)/10*h):h:b;
x2 = (a+h/2+sqrt(5)/10*h):h:b;
if n == 1
    I = (b-a)/12 * (f(a) + 5* f(x1) +5* f(x2)+f(b));
else
    I =h/12 *(f(a) + 5 *sum( f(x1) + f(x2))+2*sum(f(a+h:h:b-h))+f(b));
end
end

function I = twopntgauss(f,a,b,n) 
h = (b-a)/n;
x1 = (a+h/2-sqrt(3)/6 *h) :h:b;
x2 = (a+h/2+sqrt(3)/6*h):h:b;
I = h/2 * (sum(f(x1)+f(x2)));
end

function y = ff(x)
n = length(x);
y = zeros(n,1);
for i = 1:n
    if x(i)<1/3 && x(i)>=0
        y(i) = 0;
    elseif x(i)>=1/3 && x(i)<=1
        y(i) = 1;
    end
end
end
