Nt = 300;
dt = 0.08;
iteration = 1e7;
epsilon= 1e-5;
everyEL = 100000;

% c: depends on the ocean vehicle
c = 1;
h = 2;
% xv: angle of vanishing stability
xv = pi/4;
gamma = 2;
kone = 0.2;
ktwo = 0.2;

x_ = linspace(0,xv,Nt);
z_ = linspace(0,0.5*gamma*xv^2,Nt);
x = x_';
z = z_';
xdot = zeros(Nt,1);
xdotdot = zeros(Nt,1);
xdotdotdotdot = zeros(Nt,1);
zdot = zeros(Nt,1);
zdotdot = zeros(Nt,1);
zdotdotdotdot = zeros(Nt,1);
Gx = @(x,xdot,xdotdot,z,zdot,zdotdot) (3*h*gamma^2*x.^2-6*c*x.^2/xv^2+2*c-2*h*gamma*z-ktwo^2).*xdotdot -2*h*gamma*x.*zdotdot ...
      - (kone+ktwo)*h*gamma*x.*zdot - 2*h*gamma*zdot.*xdot + (h*gamma^2/2-c/xv^2)*6*x.*xdot.^2 - gamma*h^2*x.*(z-gamma*x.^2/2) ...
      + ((h*gamma^2/2-c/xv^2)*3*x.^2+c-h*gamma*z).*((h*gamma^2/2-c/xv^2)*x.^3+(c-h*gamma*z).*x);
Gz = @(x,xdot,xdotdot,z,zdot,zdotdot) (2*h-kone^2)*zdotdot - 2*h*gamma*x.*xdotdot + x.*xdot*h*gamma*(kone-ktwo) - gamma*h*xdot.^2 ...
      + h^2*(z-gamma*x.^2/2).*(1+gamma^2*x.^2) + c*h*gamma*x.^2.*(x.^2/xv^2-1);

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

    xdotdot = [(-x(4)+4*x(3)-5*x(2)+2*x(1))/dt^2; (x(3:end) - 2*x(2:end-1) + x(1:end-2))/dt^2; (2*x(end)-5*x(end-1)+4*x(end-2)-x(end-3))/dt^2];

    zdot = [0; (z(3:end)-z(1:end-2))/(2*dt); 0];

    zdotdot = [(-z(4)+4*z(3)-5*z(2)+2*z(1))/dt^2; (z(3:end) - 2*z(2:end-1) + z(1:end-2))/dt^2; (2*z(end)-5*z(end-1)+4*z(end-2)-z(end-3))/dt^2];

    x = A*(x-epsilon*Gx(x,xdot,xdotdot,z,zdot,zdotdot));
    z = A*(z-epsilon*Gz(x,xdot,xdotdot,z,zdot,zdotdot));

    x(1)=0; x(end)=xv;
    x(2)=3*x(1)/4 + x(3)/4;
    x(end-1)=3*x(end)/4 + x(end-2)/4;
    z(1)=0; z(end)=gamma*xv^2/2;
    z(2)=3*z(1)/4 + z(3)/4;
    z(end-1)=3*z(end)/4 + z(end-2)/4;

    if mod(iter,everyEL)==0
      figure(1)
	  plot(x, xdot, 'x-')
	  drawnow
    end

    if mod(iter,everyEL)==0
      figure(2)
	  plot(z,zdot,'x-')
      drawnow
    end
end

figure(1)
hold on
title('phase diagram of roll')
xlabel('configuration of roll')
ylabel('momentum of roll')
hold off

figure(2)
hold on
title('phase diagram of heave')
xlabel('configuration of heave')
ylabel('momentum of heave')
hold off

figure(3)
hold on
title('roll-heave coupling')
xlabel('configuration of roll')
ylabel('configuration of heave')
plot(x,z,'x-')
drawnow
hold off