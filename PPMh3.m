function [F,X,P] = PPMh3(q,dx,u,dt,uniform)
% [F,X,P] = PPMh3(q,dx,u,dt)
%
% Piecewise Parabolic Method (PPM) of Colella and Woodward, 1984.
% H.T. Huynh, Schemes and constraints for advection, in: Proceedings of the Fifteenth International Conference on Numerical Methods in Fluid Dynamics, Monterey, CA, USA, 24?28 June 1996, 1996.
%
% Use the netcdf convention for indexing (j,i). For 1-dimensional applications
% data should be sized (1,ni). For 2-dimensional applications data should be
% sized (nj,ni) and this function operates on the second index.
%
% Inputs:
%  q is the scalar field to be reconstructed/transported.
%  dx is the cell widths, same shape/size as q.
%  u is the flow, size(u)=size(q)+[0 1].
%  dt is the time-step (scalar value).
%  uniform defaults to 0. If =1, uses the uniform grid discretization.
%
% Outputs:
%  F is the flux = u q that q = q - dt*diff(F)./dx evolves the scalar field.
%    F has the shape/size of u.
%  X, P are position, values for visualization. Plot with plot(X,P).
%    X, P may have arbitrary lengths compared to q.

if ~exist('uniform','var'); uniform=0; end

sz = size(q);

% Edge values
qL = q(:,[end 1:end-1]); qR = q(:,[2:end 1]);
if uniform
	aL = (5*q + 2*qL - qR)/6;
	aR = (5*q + 2*qR - qL)/6;
else
	h0 = dx(:,[end 1:end-1]); h2 = dx(:,[2:end 1]);
	h01 = h0 + dx; h12 = dx + h2 ; h012 = dx + (h0+h2);
	aL = h12.*( dx.*qL + h0.*q )./( h01.*h012 ) ...
		+ ( ( 2*dx + h2 ).* h0.*q - h0.*dx.*qR )./( h12.*h012 );
	aR = ( ( h0 + 2*dx ).*h2.*q - dx.*h2.*qL )./( h01.*h012 ) ...
		+  h01.*( h2.*q + dx.*qR )./( h12.*h012 );
end

% Bound edge values
aL = max( min(qL,q), aL); aL = min( max(qL,q), aL);
aR = max( min(qR,q), aR); aR = min( max(qR,q), aR);

% Colella-Woodward Limiting
left = find( (aR-aL).*(6*q-3*(aR+aL)) < -(aR-aL).^2 );
aR(left) = 3*q(left) - 2*aL(left);
right = find( (aR-aL).*(6*q-3*(aR+aL)) > (aR-aL).^2 );
aL(right) = 3*q(right) - 2*aR(right);
extrema = find( (qR-q).*(q-qL) <=0 );
aL(extrema) = q(extrema); aR(extrema) = q(extrema);

% Curvature
a6 = 3*( 2*q - (aL+aR));

% Flux
Cp = u(:,2:end)./dx*dt; % CFL for use when u>0
Cm = -u(:,1:end-1)./dx*dt; % CFL for use when u<0
Ap = aR - Cp/2 .* ( (aR-aL) - a6 .* (1 - 2/3*Cp) ); % Average value leaving to right for u>0
Am = aL + Cm/2 .* ( (aR-aL) + a6 .* (1 - 2/3*Cm) ); % Average value leaving to left for u<0
% Combine fluxes from different signed flow
up = (u+abs(u))/2; um = (u-abs(u))/2;
F = up.*Ap(:,[end 1:end])+um.*Am(:,[1:end 1]);

if nargout > 1
	% Create plottable reconstruction
	d = dx(:,[1 1:end]); d(:,1) = 0; xg = cumsum(d,2);
	n=7; % Number of nodes per cell to plot
	X = zeros([sz(1) n sz(2)]); P = X*NaN;
	for j=1:n
		x=(j-1)/(n-1);
		X(:,j,:) = (1-x) * xg(:,1:end-1) + x * xg(:,2:end);
		P(:,j,:) = (1-x) * aL + x * aR + a6 .* x .* (1-x);
	end

	X = reshape(X,sz.*[1 n]);
	P = reshape(P,sz.*[1 n]);
end
