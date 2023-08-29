
% c: depends on the ocean vehicle
c = 2;

% xv: angle of vanishing stability
xv = pi/4;

ktwo = 0.5;

%initialising
dt = 0.05;
Nt = 150;
total_time = dt * Nt;
iterations = 1e7;
alpha = dt^2/50;
everyPlot = 10000;

x = linspace(0,xv,Nt);


for iter = 1:iterations
  xdot = [0, (x(3:end)-x(1:end-2))/(2*dt), 0];

  xdotdot = [(-x(4)+4*x(3)-5*x(2))/dt^2, (x(3:end) - 2*x(2:end-1) + x(1:end-2))/dt^2, (2*xv-5*x(end-1)+4*x(end-2)-x(end-3))/dt^2];

  %xdotdotdot = [(-3*x(5)/2+7*x(4)-12*x(3)+9*x(2)-5*x(1)/2)/dt^3,((-3*x(6)/2+7*x(5)-12*x(4)+9*x(3)-5*x(2)/2))/dt^3,...
  %    (x(5:end)/2-x(4:end-1)+x(2:end-3)-x(1:end-4)/2)/dt^3,...
  %    (5*x(end-1)/2-9*x(end-2)+12*x(end-3)-7*x(end-4)+3*x(end-5)/2)/dt^3,(5*x(end)/2-9*x(end-1)+12*x(end-2)-7*x(end-3)+3*x(end-4)/2)/dt^3];

  xdotdotdotdot = [(-2*x(6)+11*x(5)-24*x(4)+26*x(3)-14*x(2))/dt^4,(-x(6)+6*x(5)-14*x(4)+16*x(3)-9*x(2))/dt^4,...
      (x(5:end)-4*x(4:end-1)+6*x(3:end-2)-4*x(2:end-3)+x(1:end-4))/dt^4,(2*xv-9*x(end-1)+16*x(end-2)-14*x(end-3)+6*x(end-4)-x(end-5))/dt^4,...
      (3*xv-14*x(end-1)+26*x(end-2)-24*x(end-3)+11*x(end-4)-2*x(end-5))/dt^4];
 

  dEL = xdotdotdotdot + 2*(c-3*c*x.^2/xv^2-ktwo^2/2).*xdotdot - 6*c*x/xv^2.*xdot.^2 + 3*c^2*x.^5/xv^4 - 4*c^2*x.^3/xv^2 + c^2*x;
  
    x(1)=0;
    x(2)=3*x(1)/4 + x(3)/4;
    x(end)= xv;
    x(end-1)=3*x(end)/4 + x(end-2)/4;
 
 
  x = x - alpha*dEL;
  
  if mod(iter,everyPlot)==0
   figure(1)
   plot(x,xdot,'x-')      
   drawnow
  end

  %if mod(iter,everyPlot)==0
  % figure(2)
  % plot(x)      
  % drawnow
  %end

  %if mod(iter,everyPlot)==0
  % figure(3)
  % plot(xdotdotdotdot)      
  % drawnow
  %end
end
