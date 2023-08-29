Nt = 150;
dt = 0.08;
iteration = 1e7;
epsilon= 1e-5;
everyPlot = 2000;
everyEL = 5000;

xv = pi/4;
ktwo = 0.2;
c = 2;

x_ = linspace(0,xv,Nt);
x = x_';
xdot = zeros(Nt,1);
xdotdot = zeros(Nt,1);
xdotdotdotdot = zeros(Nt,1);
elvalue = [];
G = @(x,xdot,xdotdot) 2*(c-3*c*x.^2/xv^2-ktwo^2/2).*xdotdot - 6*c*x/xv^2.*xdot.^2 + 3*c^2*x.^5/xv^4 - 4*c^2*x.^3/xv^2 + c^2*x;
F = @(x,xdot,xdotdot,xdotdotdotdot) xdotdotdotdot + G(x,xdot,xdotdot);

%construct the fourth order operator matrix
a1 = [3,-14,26,-24,11,-2,zeros(1,Nt-6)];
a2 = [2,-9,16,-14,6,-1,zeros(1,Nt-6)];
abend = [zeros(1,Nt-6),-1,6,-14,16,-9,2];
aend = [zeros(1,Nt-6),-2,11,-24,26,-14,3];
D = [a1;a2];
for ii = 1:Nt-4
    a = zeros(1,Nt);
    a(ii:ii+4) = [1,-4,6,-4,1];
    D = [D;a];
end
D = [D;abend;aend];

I = eye(Nt);
A = inv(I+epsilon*D/dt^4);


for iter = 1:iteration
    xdot = [0; (x(3:end)-x(1:end-2))/(2*dt); 0];

    xdotdot = [(-x(4)+4*x(3)-5*x(2))/dt^2; (x(3:end) - 2*x(2:end-1) + x(1:end-2))/dt^2; (2*xv-5*x(end-1)+4*x(end-2)-x(end-3))/dt^2];

    x = A*(x-epsilon*G(x,xdot,xdotdot));

    x(1) = 0;
    x(2)=3*x(1)/4 + x(3)/4;
    x(end) = xv;
    x(end-1)=3*x(end)/4 + x(end-2)/4;

    if mod(iter,everyPlot)==0
	  plot(x, xdot, 'x-')
	  drawnow
    end

    %if mod(iter,everyEL)==0
	%  xdotdotdotdot = D*x/dt^4;
    %  value = sum(F(x,xdot,xdotdot,xdotdotdotdot).^2);
    %  elvalue = [elvalue,sqrt(value)];
    %end
end
