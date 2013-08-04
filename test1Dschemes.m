function test1Dschemes(n, CFL, nTurns, flipFlow, shape, varGrid)
% testSchemes(n,CFL,nTurns,flipFlow)
%
% n        - number of points (default 35)
% CFL      - Courant-Friedrichs-Lewy number to use (CFL=u*dt/dx)
% nTurns   - number of rotations through domain (default 1)
% flipFlow - If =1, change sign of flow each rotation
% shape    - string determining shape of test function
%            (see 'help testFunction')
% varGrid  - If <1, generate a random grid s.t. varGrid/n < dx < 1/n
%            and CFL value applies at dx=varGrid/n
%
% To run a simple scheme
% >> test1Dschemes
%
% To run specific configuration
% >> test1Dschemes(100,.05,3,1,'cosinebell')

% The following allows invoking test1Dschemes without all arguments
if ~exist('n','var'); n=35; end
if ~exist('CFL','var'); CFL=.1; end
if ~exist('nTurns','var'); nTurns=1; end
if ~exist('flipFlow','var'); flipFlow=0; end
if ~exist('shape','var'); shape='triangle'; end
if ~exist('varGrid','var'); varGrid=.3; end

dx=rand(1,n); dx=varGrid/n+(1-varGrid)*dx/sum(dx);
xg=cumsum([0 dx]); xc=(xg(1:n)+xg(2:n+1))/2;
u=ones(1,n+1);
nt=round( 1./( CFL*min(dx)./max(abs(u)) ) );
dt=1/nt;

Q=testFunction(xc,shape);
x0=(0:1000)/1000; q0=testFunction(x0,shape);

[~,Xpcm,Ppcm]=PCM(Q,dx,u,dt);
[~,Xplm,Pplm]=PLM(Q,dx,u,dt);
[~,Xppmh3,Pppmh3]=PPMh3(Q,dx,u,dt);
[~,XppmCW,PppmCW]=PPMcw(Q,dx,u,dt);
plot(x0,q0,'k:', ...
	Xpcm,Ppcm,'r', ...
	Xplm,Pplm,'m', ...
	Xppmh3,Pppmh3,'b', ...
	XppmCW, PppmCW,'k')
legend('Test function','PCM','PLM','PPMh3','PPMcw')

qPCM=Q; qPLM=Q; qPPMh3=Q; qPPMcw=Q;
for t=1:nt*nTurns
	if flipFlow && mod(t,nt)==0
		u=-u;
	end
	Fpcm=PCM( qPCM, dx, u, dt );
	qPCM = update( qPCM, Fpcm, dx , dt, 'PCM' );
	Fplm = PLM( qPLM, dx, u, dt );
	qPLM = update( qPLM, Fplm, dx , dt, 'PLM' );
	FppmH3 = PPMh3( qPPMh3, dx, u, dt );
	qPPMh3 = update( qPPMh3, FppmH3, dx , dt, 'PPMh3' );
	FppmCW = PPMcw( qPPMcw, dx, u, dt );
	qPPMcw = update( qPPMcw, FppmCW, dx , dt, 'PPMcw' );
	if mod(t,10)==0
		plot(xc,qPCM,'r.-',xc,qPLM,'m.-',xc,qPPMh3,'b.-',xc,qPPMcw,'k.-')
		drawnow
	end
end
plot(x0,q0,'k:',xc,qPCM,'r.-',xc,qPLM,'m.-',xc,qPPMh3,'b.-',xc,qPPMcw,'k.-')
legend('Test function','PCM','PLM','PPMh3','PPMcw')

function [Qnp1] = update(Qn, F, dx, dt, msg)
% Implements the time update with some sanity checking
Qmin=min(Qn(:)); Qmax=max(Qn(:)); Q2=dx.*Qn.*Qn; Q2=sum(Q2(:));
Qnp1 = Qn - dt*diff(F)./dx;
if min(Qnp1(:))-Qmin < -eps(Qmin)*0
	fprintf('Undershoot %e - %e = %e %s\n',min(Qnp1(:)),Qmin,min(Qnp1(:))-Qmin,msg)
% 	keyboard
end
if max(Qnp1(:))-Qmax > eps(Qmax)*0
	fprintf('Overshoot %e - %e = %e %s\n',max(Qnp1(:)),Qmax,max(Qnp1(:))-Qmax,msg)
% 	keyboard
end
Qn2=dx.*Qnp1.*Qnp1; Qn2=sum(Qn2(:));

