
% c: depends on the ocean vehicle
c = 2;
h = 2;
% xv: angle of vanishing stability
xv = pi/4;
gamma = 2;
kone = 0.2;
ktwo = 0.2;

%initialising
dt = 0.08;
Nt = 150;
total_time = Nt*dt;
iterations = 1e7;
alpha = 5e-7;
everyPlot = 1000;

x = linspace(0,xv,Nt);
z = linspace(0,gamma*xv^2/2,Nt);

for iter = 1:iterations

  xdot = [0, (x(3:end)-x(1:end-2))/(2*dt), 0];

  xdotdot = [(-x(4)+4*x(3)-5*x(2)+2*x(1))/dt^2, (x(3:end) - 2*x(2:end-1) + x(1:end-2))/dt^2, (2*x(end)-5*x(end-1)+4*x(end-2)-x(end-3))/dt^2];

  %xdotdotdot = [(-3*x(5)/2+7*x(4)-12*x(3)+9*x(2)-5*x(1)/2)/dt^3,(-3*x(6)/2+7*x(5)-12*x(4)+9*x(3)-5*x(2)/2)/dt^3,...
      %(x(5:end)/2-x(4:end-1)+x(2:end-3)-x(1:end-4)/2)/dt^3,...
      %(5*x(end-1)/2-9*x(end-2)+12*x(end-3)-7*x(end-4)+3*x(end-5)/2)/dt^3,(5*x(end)/2-9*x(end-1)+12*x(end-2)-7*x(end-3)+3*x(end-4)/2)/dt^3];

  xdotdotdotdot = [(-2*x(6)+11*x(5)-24*x(4)+26*x(3)-14*x(2)+3*x(1))/dt^4,(-x(6)+6*x(5)-14*x(4)+16*x(3)-9*x(2)+2*x(1))/dt^4,...
      (x(5:end)-4*x(4:end-1)+6*x(3:end-2)-4*x(2:end-3)+x(1:end-4))/dt^4,(2*x(end)-9*x(end-1)+16*x(end-2)-14*x(end-3)+6*x(end-4)-x(end-5))/dt^4,...
      (3*x(end)-14*x(end-1)+26*x(end-2)-24*x(end-3)+11*x(end-4)-2*x(end-5))/dt^4];
 
  zdot = [0, (z(3:end)-z(1:end-2))/(2*dt), 0];

  zdotdot = [(-z(4)+4*z(3)-5*z(2)+2*z(1))/dt^2, (z(3:end) - 2*z(2:end-1) + z(1:end-2))/dt^2, (2*z(end)-5*z(end-1)+4*z(end-2)-z(end-3))/dt^2];

  %zdotdotdot = [(-3*z(5)/2+7*z(4)-12*z(3)+9*z(2)-5*z(1)/2)/dt^3,(-3*z(6)/2+7*z(5)-12*z(4)+9*z(3)-5*z(2)/2)/dt^3,...
      %(z(5:end)/2-z(4:end-1)+z(2:end-3)-z(1:end-4)/2)/dt^3,...
      %(5*z(end-1)/2-9*z(end-2)+12*z(end-3)-7*z(end-4)+3*z(end-5)/2)/dt^3,(5*z(end)/2-9*z(end-1)+12*z(end-2)-7*z(end-3)+3*z(end-4)/2)/dt^3];

  zdotdotdotdot = [(-2*z(6)+11*z(5)-24*z(4)+26*z(3)-14*z(2)+3*z(1))/dt^4,(-z(6)+6*z(5)-14*z(4)+16*z(3)-9*z(2)+2*z(1))/dt^4,...
      (z(5:end)-4*z(4:end-1)+6*z(3:end-2)-4*z(2:end-3)+z(1:end-4))/dt^4,(2*z(end)-9*z(end-1)+16*z(end-2)-14*z(end-3)+6*z(end-4)-z(end-5))/dt^4,...
      (3*z(end)-14*z(end-1)+26*z(end-2)-24*z(end-3)+11*z(end-4)-2*z(end-5))/dt^4];

  dEL_x = xdotdotdotdot + (3*h*gamma^2*x.^2-6*c*x.^2/xv^2+2*c-2*h*gamma*z-ktwo^2).*xdotdot -2*h*gamma*x.*zdotdot - (kone+ktwo)*h*gamma*x.*zdot ...
      - 2*h*gamma*zdot.*xdot + (h*gamma^2/2-c/xv^2)*6*x.*xdot.^2 - gamma*h^2*x.*(z-gamma*x.^2/2) ...
      + ((h*gamma^2/2-c/xv^2)*3*x.^2+c-h*gamma*z).*((h*gamma^2/2-c/xv^2)*x.^3+(c-h*gamma*z).*x);

  dEL_z = zdotdotdotdot + (2*h-kone^2)*zdotdot - 2*h*gamma*x.*xdotdot + x.*xdot*h*gamma*(kone-ktwo) - gamma*h*xdot.^2 ...
      + h^2*(z-gamma*x.^2/2).*(1+gamma^2*x.^2) + c*h*gamma*x.^2.*(x.^2/xv^2-1);


  x(1)=0; x(end)=xv;
  x(2)=3*x(1)/4 + x(3)/4;
  x(end-1)=3*x(end)/4 + x(end-2)/4;
  z(1)=0; z(end)=gamma*xv^2/2;
  z(2)=3*z(1)/4 + z(3)/4;
  z(end-1)=3*z(end)/4 + z(end-2)/4;
 
  x = x - alpha*dEL_x;
  z = z - alpha*dEL_z;

  if mod(iter,everyPlot)==0
   figure(1)
   plot(x,xdot,'x-')      
   drawnow
  end

  if mod(iter,everyPlot)==0
   figure(2)
   plot(z,zdot,'x-')      
   drawnow
  end

end
