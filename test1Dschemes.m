function test1Dschemes(n,CFL,nTurns,flipFlow,shape)
% testSchemes(n,CFL,nTurns,flipFlow)
%
% n        - number of points (default 35)
% CFL      - Courant-Friedrichs-Lewy number to use (CFL=u*dt/dx)
% nTurns   - number of rotations through domain (default 1)
% flipFlow - If =1, change sign of flow each rotation
% shape    - string determining shape of test function
%            (see 'help testFunction')
%
% To run a simple scheme
% >> test1Dschemes
%
% To run specific configuration
% >> test1Dschemes(100,.05,3,1,'cosinebell')

% This is allow invocing testSchemes without arguments
if ~exist('n','var'); n=35; end
if ~exist('CFL','var'); CFL=.1; end
if ~exist('nTurns','var'); nTurns=1; end
if ~exist('flipFlow','var'); flipFlow=0; end
if ~exist('shape','var'); shape='triangle'; end

xg=(0:n)/n; dx=diff(xg); xc=(xg(1:n)+xg(2:n+1))/2;
nt=n/CFL; dt=1/nt;
u=CFL/dt*(1/n)+0*xg;

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
	[Fpcm]=PCM(qPCM,dx,u,dt);
	qPCM = qPCM - diff(Fpcm);
	[Fplm]=PLM(qPLM,dx,u,dt);
	qPLM = qPLM - diff(Fplm);
	[Fppmh3]=PPMh3(qPPMh3,dx,u,dt);
	qPPMh3 = qPPMh3 - diff(Fppmh3);
	[FppmCW]=PPMcw(qPPMcw,dx,u,dt);
	qPPMcw = qPPMcw - diff(FppmCW);
	if mod(t,10)==0
		plot(xc,qPCM,'r',xc,qPLM,'m',xc,qPPMh3,'b',xc,qPPMcw,'k')
		drawnow
	end
end
plot(x0,q0,'k:',xc,qPCM,'r',xc,qPLM,'m',xc,qPPMh3,'b',xc,qPPMcw,'k')
legend('Test function','PCM','PLM','PPMh3','PPMcw')

