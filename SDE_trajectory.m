
steps = 2400;
N = 2000;
dt = 0.01;
total_time = steps*dt;
eps = 0.005;

% everyPlot = 3000;

u = zeros(N,1);
v = zeros(N,1);
x = zeros(N,1);
y = zeros(N,1);

us = zeros(N,steps);
vs = zeros(N,steps);
xs = zeros(N,steps);
ys = zeros(N,steps);

for step = 1:steps
  u_ = u + dt*bu(u,v,x,y) + sqrt(dt*eps)*randn(size(x));
  v_ = v + dt*bv(u,v,x,y) + sqrt(dt*eps)*randn(size(x));
  x_ = x + dt*bx(u,v,x,y) + sqrt(dt*eps)*randn(size(x));
  y_ = y + dt*by(u,v,x,y) + sqrt(dt*eps)*randn(size(x));
  x = x_; v = v_; u = u_; y = y_;

  x(find(abs(x)>2))=sign(x(find(abs(x)>2)))*2;
  y(find(abs(y)>2))=sign(y(find(abs(y)>2)))*2;

  us(:,step) = u;
  vs(:,step) = v;
  xs(:,step) = x;
  ys(:,step) = y;
  
end

figure(1)
plot(us(find(u>1.5),:)',vs(find(u>1.5),:)','DisplayName','Simulation')
drawnow

figure(2)
plot(xs(find(u>1.5),:)',ys(find(u>1.5),:)','DisplayName','Simulation')
drawnow

figure(4)
plot(xs(find(u>1.5),:)',us(find(u>1.5),:)','DisplayName','Simulation')

  