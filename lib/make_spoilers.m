function [gxSpoil, gySpoil, gzSpoil] = make_spoilers(Nx, Ny, Nz, fov, NcyclesSpoil, sys, CRT)
%% make_spoilers  Create x/y/z spoiler trapezoids for EPI sequence.
%   Nx, Ny, Nz:    matrix dimensions
%   fov:           field of view [x y z] (m)
%   NcyclesSpoil:  number of gradient spoiler cycles across the FOV
%   sys:           Pulseq system struct (mr.opts)
%   CRT:           common raster time for GE compatibility (s)

deltak = 1./fov;
tmp = 0.5; % scale factor <1 to avoid PNS
gxSpoil = trap4ge(mr.scaleGrad(mr.makeTrapezoid('x', sys, 'Area', Nx*deltak(1)*NcyclesSpoil/tmp), tmp, sys), CRT, sys);
gySpoil = trap4ge(mr.scaleGrad(mr.makeTrapezoid('y', sys, 'Area', Ny*deltak(2)*NcyclesSpoil/tmp), tmp, sys), CRT, sys);
gzSpoil = trap4ge(mr.scaleGrad(mr.makeTrapezoid('z', sys, 'Area', Nz*deltak(3)*NcyclesSpoil/tmp), tmp, sys), CRT, sys);
end
