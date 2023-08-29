c = 2;
xv = pi/4;
ktwo = 0.2;

V = @(x) 0.5*c*x.^2 - 0.25*c*x.^4/xv^2;
Hamilton = @(x,y) V(x) + 0.5*y.^2;
bx = @(x,y) y;
by = @(x,y) c*x.^3/xv^2 - c*x - ktwo*y;
dbxdy = 1;
dbydx = @(x) 3*c*x.^2/xv^2 -c;
dbydy = -ktwo;

H = @(phix,phiy,thetax,thetay) bx(phix,phiy).*thetax + by(phix,phiy).*thetay + 0.5*(thetay.^2);
dH_dphix = @(phix,thetay) dbydx(phix).*thetay;
dH_dphiy = @(thetax,thetay) dbxdy*thetax + dbydy*thetay;
dH_dthetax = @(phix,phiy) bx(phix,phiy);
dH_dthetay = @(phix,phiy,thetay) by(phix,phiy) + thetay;

params.bx = bx;
params.by = by;
params.dH_dphix = dH_dphix;
params.dH_dphiy = dH_dphiy;
params.dH_dthetax = dH_dthetax;
params.dH_dthetay = dH_dthetay;

[x_,y_] = meshgrid(linspace(-1.5,1.5,100), linspace(-1.5,1.5,100));

Nt = 100;
T = 10;
t = linspace(0,T,Nt);
dt = t(2)-t(1);

params.Nt = Nt;
params.dt = dt;

phi0x = 0;
phi0y = 0;
phiTx = xv;
phiTy = 0;

params.phi0x = phi0x;
params.phi0y = phi0y;
params.phiTx = phiTx;
params.phiTy = phiTy;

phix = linspace(phi0x,phiTx,Nt);
phiy = linspace(phi0y,phiTy,Nt);
mux = zeros(1,params.Nt-1);
muy = zeros(1,params.Nt-1);
thetay = zeros(1,params.Nt-1);

params.lambda = 0.1;
params.betax = 0;
params.betay = 0;

for outer=1:100
  for iter=1:1e6
	gradx = gradient(thetay, params);          % compute gradient

  if iter == 25
     gradcompare = zeros(1,Nt-1);
     h = 1e-7;
     for i = 1:Nt-1
         vector = zeros(1,Nt-1);
         vector(i) = 1;
         costA = cost(thetay-h*vector,params);
         costB = cost(thetay+h*vector,params);
         gradcompare(i) = (costB-costA)/2/h;
     end
     %const(thetay,params)
     figure(4)
     plot(gradx,'x-'); hold on
     plot(gradcompare,'g-'); hold off
     drawnow
  end
  
	dtau = linesearch(thetay, -gradx, params); % get best step size
	thetay = thetay - dtau*gradx;              % gradient descent

	if mod(iter,25)==0
	  % plotting
	  for s=2:Nt
	    phix(s) = phix(s-1) + params.dt*params.dH_dthetax(phix(s-1),phiy(s-1));
	    phiy(s) = phiy(s-1) + params.dt*params.dH_dthetay(phix(s-1),phiy(s-1),thetay(s-1));
      end
	  figure(1)
      contourf(x_,y_,Hamilton(x_,y_), 13); hold on; 
	  plot(phix, phiy, 'k-x'); hold off
	  drawnow
	  action = 0.5*sum(thetay.^2*dt);
	  display(sprintf('action=%g, norm=%g, dtau=%g', action, norm(gradx),dtau))
	  drawnow
	end
	if norm(gradx)<5e-3
	  break
	end
  end
  params.betax = params.betax + params.lambda*(phix(end)-phiTx);
  params.betay = params.betay + params.lambda*(phiy(end)-phiTy);
  params.lambda = params.lambda*1.1;
  display(sprintf('OUTER LOOP NEXT STEP after %d steps: lambda=%g, betax=%g, betay=%g',...
				  iter, params.lambda, params.betax, params.betay));
end

function ret=cost(thetay, params)
  phix = zeros(1,params.Nt);
  phiy = zeros(1,params.Nt);
  phix(1) = params.phi0x;
  phiy(1) = params.phi0y;
  %integrate the phi from time 0 to time T
  for s=2:params.Nt
	phix(s) = phix(s-1) + params.dt*params.dH_dthetax(phix(s-1),phiy(s-1));
	phiy(s) = phiy(s-1) + params.dt*params.dH_dthetay(phix(s-1),phiy(s-1),thetay(s-1));
  end
  
  %cost function (theta + boundary condition and panelty)
  ret = 0.5*sum(thetay.^2)*params.dt  ...
      + params.lambda*(phix(end)-params.phiTx).^2 + params.betax*(phix(end)-params.phiTx) ...
      + params.lambda*(phiy(end)-params.phiTy).^2 + params.betay*(phiy(end)-params.phiTy);
end

function ret=gradient(thetay, params)
  phix = zeros(1,params.Nt);
  phiy = zeros(1,params.Nt);
  phix(1) = params.phi0x;
  phiy(1) = params.phi0y;
  for s=2:params.Nt  % forward equation
	phix(s) = phix(s-1) + params.dt*params.dH_dthetax(phix(s-1),phiy(s-1));
	phiy(s) = phiy(s-1) + params.dt*params.dH_dthetay(phix(s-1),phiy(s-1),thetay(s-1));
  end
  mux=zeros(1,params.Nt-1);
  muy=zeros(1,params.Nt-1);
  %why? mu = partial boundary condition partial phi
  mux(params.Nt-1) = -(2*params.lambda*(phix(end)-params.phiTx)+params.betax);
  muy(params.Nt-1) = -(2*params.lambda*(phiy(end)-params.phiTy)+params.betay);
  for s=params.Nt-1:-1:2  % backward/adjoint equation
	mux(s-1) = mux(s) + params.dt*params.dH_dphix(phix(s),muy(s));
	muy(s-1) = muy(s) + params.dt*params.dH_dphiy(mux(s),muy(s));
  end
  ret = (thetay - muy)*params.dt; % gradient computation
end

%deciding for eac h iteration what is the optimal step size (alpla)
function ret=linesearch(thetay, p, params)
  dtau = 1;
  for kk=1:100
   	thetay_ = thetay + dtau*p;
   	if cost(thetay_,params) < cost(thetay,params)
   	  break
   	else
   	  dtau = dtau*0.5;
   	end
  end
  ret = dtau;
end
