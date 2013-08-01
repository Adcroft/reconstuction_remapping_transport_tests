function [f] = testFunction(x,shape)
% f = testFunction(x,shape)
%
% Generates a test function as function of x
% shape can be any of 'cosinebell', 'wave', 'box', 'pulse'
% 'triangle'
%
% plot(0:.01:1, testFunction(0:.01:1,'pulse') )

switch shape
	case {'cosinebell'}
		f = cosbell(x, .5, .3);
	case {'wave'}
		f = cosbell(x, .5, .5);
	case {'box'}
		f = box(x, .5, .25);
	case {'pulse'}
		f = box(x, .5, .1);
	case {'triangle'}
		f = triangle(x, .5, .2);
	otherwise
		error('Unkown shape')
end

function [f] = box(x, xc, w)
% Chapeau hat function centered at xc with half width w
xp = w-abs(x-xc);
f = (1+sign(xp))/2;

function [f] = triangle(x, xc, w)
% Triagular function centered at xc with half width w
xp = (w-abs(x-xc))/w;
f = max(0, xp);

function [f] = cosbell(x, xc, w)
% Cosine bell centered at xc with width w
xp = max(-1, min(1, (x-xc)/w))*pi;
f = (1+cos(xp))/2;