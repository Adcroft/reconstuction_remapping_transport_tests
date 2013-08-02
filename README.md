Tests of transport schemes that use reconstruction/remapping ideas, implemented
in Matlab.

To run a simple test

	>> test1Dschemes

To go further

	>> help test1Dschemes

Protocol for reconstruction functions is to use the netcdf convention for indexing (j,i).
For 1-dimensional applications data should be sized (1,ni).
For 2-dimensional applications data should be sized (nj,ni) and functions operate on the second index.
To operate in the j-direction, transpose all data before calling the function and transpose back.
The reconstruction functions should return both flux and plotting data, e.g.

>> [F,X,P] = PCM(q,dx,u,dt);

so that the evolution can be found by

>> q = q - diff(F);

and the reconstruction can be visualized by

>> plot(X,P)
