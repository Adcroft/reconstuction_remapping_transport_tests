function [F,X,P] = PCM(q,dx,u,dt)
% [F,X,P] = PCM(q,dx,u,dt)

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