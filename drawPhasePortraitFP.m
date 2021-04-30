function drawPhasePortraitFP(a,b,g,N,mySS,evec1,evec2)
hold on;
f = @(t,y) [-1*(-g-a*y(1)+a*(2-2*y(1)-y(2)))/(2*N)-a*y(1)*(2-2*y(1)-y(2))+g*y(2); a*y(1)*(2-2*y(1)-y(2))+2*g*(1-y(1)-y(2))-g*y(2)-(b+0.5*a*(2-2*y(1)-y(2)))*y(2)-(b-g+a*y(1)-0.5*a*(2-2*y(1)-y(2))-a/2*y(2))/(2*N)];
f1 = @(x,y) -1*(-g-a*x+a*(2-2*x-y))/(2*N)-a*x*(2-2*x-y)+g*y;
f2 = @(x,y) a*x*(2-2*x-y)+2*g*(1-x-y)-g*y-(b+0.5*a*(2-2*x-y))*y-(b-g+a*x-0.5*a*(2-2*x-y)-a/2*y)/(2*N);
y1 = linspace(mySS(1)-0.15, mySS(1)+0.15,50);
y2 = linspace(mySS(2)-0.15, mySS(2)+0.15,50);
[x,y] = meshgrid(y1,y2);

u = zeros(size(x));
v = zeros(size(x));

t=0;
for i=1:numel(x)
    Yprime = f(t,[x(i);y(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end

% for i = 1:numel(x)
% Vmod = sqrt(u(i)^2 + v(i)^2);
% u(i) = u(i)/Vmod;
% v(i) = v(i)/Vmod;
% end

h1=quiver(x,y,u,v,'k');
set(h1,'AutoScale','on', 'AutoScaleFactor', 2)
% legend('Phase portrait')
figure(gcf);
xlabel('y_1');
ylabel('y_2');

% xeig = linspace(0,1,1000);
% yeig1 = mySS(2)+evec1(2)/evec1(1)*(xeig-mySS(1));
% yeig2 = mySS(2)+evec2(2)/evec2(1)*(xeig-mySS(1));
% p(1)=plot(xeig,yeig1,'r');
% p(2)=plot(xeig,yeig2,'r');
p(1)=fimplicit(f1);
p(2)=fimplicit(f2);
p(1).LineWidth = 3;
p(2).LineWidth = 3;
    