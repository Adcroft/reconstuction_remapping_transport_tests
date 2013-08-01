function [F,X,P] = PPMh3(q,dx,u,dt)
% [F,X,P] = PPMh3(q,dx,u,dt)

sz = size(q);

% Edge values
qL = q(:,[end 1:end-1]); qR = q(:,[2:end 1]);
aL = (5*q + 2*qL - qR)/6;
aR = (5*q + 2*qR - qL)/6;

% Bound edge values
aL = max( min(qL,q), aL); aL = min( max(qL,q), aL);
aR = max( min(qR,q), aR); aR = min( max(qR,q), aR);

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
Fp = (aR - Cp/2 .* ( (aR-aL) - a6 .* (1 - 2/3*Cp) ) )./dx*dt; % Flux out of right for u>0
Fm = (aL + Cm/2 .* ( (aR-aL) + a6 .* (1 - 2/3*Cm)) )./dx*dt; % Flux out of left for u<0
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