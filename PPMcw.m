function [F,X,P] = PPMcw(q,dx,u,dt,uniform)
% [F,X,P] = PPMcw(q,dx,u,dt)
%
% Piecewise Parabolic Method (PPM) of Colella and Woodward, 1984.
% http://dx.doi.org/10.1016/0021-991(84)90143-8
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

% Slope
qL = q(:,[end 1:end-1]); qR = q(:,[2:end 1]);
sMax = max( max(qR, qL), q) - q;
sMin = q - min( min(qR, qL), q);
if uniform
	slp2 = (qR - qL)/2;
else
	hL = dx(:,[end 1:end-1]); hR = dx(:,[2:end 1]);
	slp2 = dx ./ ( dx + (hL+hR) ) .* ( ...
		( 2*hL + dx )./( hR + dx ) .* ( qR - q ) + ...
		( dx + 2*hR )./( hL + dx ) .* ( q - qL ) );
end
slp = sign(slp2).*min( abs(slp2), 2*min( sMax, sMin ) );

% Edge values
sR = slp(:,[2:end 1]);
if uniform
% 	sL = slp(:,[end 1:end-1]);
% 	aL = (q+qL)/2-(slp-sL)/6;
	aR = (q+qR)/2-(sR-slp)/6;
else
	hL = dx(:,[end 1:end-1]); hR = dx(:,[2:end 1]);
	hRR = hR(:,[2:end 1]);
	f1 = 2*hR.*dx./( dx + hR ).*( ...
		( hL + dx )./( 2*dx + hR ) - ( hRR + hR )./( 2*hR + dx ) );
	f2 = dx .* ( hL + dx )./( 2*dx + hR );
	f3 = hR .* ( hR + hRR )./( dx + 2*hR );
	aR = q + dx./( dx + hR ).*( qR - q ) ...
		+ 1./( (hL + dx) + (hR + hRR) ) .* ( ...
		 ( f1 .* ( qR - q ) - f2 .* sR ) + f3 .* slp );
end
aL = aR(:,[end 1:end-1]);

% Colella-Woodwrd Limiting
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
Fp = aR - Cp/2 .* ( (aR-aL) - a6 .* (1 - 2/3*Cp) ); % Average value leaving to right for u>0
Fm = aL + Cm/2 .* ( (aR-aL) + a6 .* (1 - 2/3*Cm) ); % Average value leaving to left for u<0
% Combine fluxes from different signed flow
up = (u+abs(u))/2; um = (u-abs(u))/2;
F = up.*Fp(:,[end 1:end])+um.*Fm(:,[1:end 1]);

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
