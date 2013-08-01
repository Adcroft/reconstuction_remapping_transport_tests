function [F,X,P] = PCM(q,dx,u,dt)
% [F,X,P] = PCM(q,dx,u,dt)
%
% Piecewise Constant Method
%   In one dimensions this is also known as the FTUS (forward in time, upstream)
% or first order upwind. It is the lowest order most diffusive explicit scheme
% possible.
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

% Flux
Fp = q./dx*dt; % Flux out of right for u>0
Fm = q./dx*dt; % Flux out of left for u<0
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
	P(:,1,:) = q;
	P(:,2,:) = q;
	X = reshape(X,sz.*[1 2]);
	P = reshape(P,sz.*[1 2]);
end
