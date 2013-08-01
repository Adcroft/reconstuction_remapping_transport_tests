function [F,X,P] = PLM(q,dx,u,dt)
% [F,X,P] = PLM(q,dx,u,dt)
%
% Piecewise Linear Method, due to van Leer, 1979.
% http://dx.doi.org/10.1016/0021-9991(79)90145-1
% This currently implements the uniform grid case.
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
%
% Outputs:
%  F is the flux * dt/dx so that q = q - diff(F) evolves the scalar field.
%    F has the shape/size of u.
%  X, P are position, values for visualization. Plot with plot(X,P).
%    X, P may have arbitrary lengths compared to q.

sz = size(q);

% Slope
qL = q(:,[end 1:end-1]); qR = q(:,[2:end 1]);
sMax = max( max(qR, qL), q) - q;
sMin = q - min( min(qR, qL), q);
slp2 = (qR - qL)/2;
slp = sign(slp2).*min( abs(slp2), 2*min( sMax, sMin ) );

% Edge values
aL = q - slp/2;
aR = q + slp/2;

% Flux
Cp = u(:,2:end)./dx*dt; % CFL for use when u>0
Cm = -u(:,1:end-1)./dx*dt; % CFL for use when u<0
Fp = (aR - slp/2 .* Cp)./dx*dt; % Flux out of right for u>0
Fm = (aL + slp/2 .* Cm)./dx*dt; % Flux out of left for u<0
% Combine fluxes from different signed flow
up = (u+abs(u))/2; um = (u-abs(u))/2;
F = up.*Fp(:,[end 1:end])+um.*Fm(:,[1:end 1]);

if nargout > 1
	% Create plottable reconstruction
	d = dx(:,[1 1:end]); d(:,1) = 0; xg = cumsum(d,2);
	X = zeros([sz(1) 2 sz(2)]);
	X(:,1,:) = xg(:,1:end-1);
	X(:,2,:) = xg(:,2:end);

	P = X*NaN;
	P(:,1,:) = aL;
	P(:,2,:) = aR;
	X = reshape(X,sz.*[1 2]);
	P = reshape(P,sz.*[1 2]);
end
