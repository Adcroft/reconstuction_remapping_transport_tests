Tests of transport schemes that use reconstruction/remapping ideas, implemented
in Matlab.

To run a simple test

	>> test1Dschemes

To go further

	>> help test1Dschemes

Protocol for reconstruction functions is to return flux and plotting data, e.g.

>> [F,X,P] = PCM(q,dx,u,dt);

so that the evolution can be found by

>> q = q - diff(F);

and the reconstruction can be visualized by

>> plot(X,P)
