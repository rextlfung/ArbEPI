function [gxPre, gyPre, gzPre] = make_prephasers(Nx, Ny, Nz, fov, sys, CRT)
%% make_prephasers  Create x/y/z prephaser trapezoids for EPI readout.
%   Nx, Ny, Nz: matrix dimensions
%   fov:        field of view [x y z] (m)
%   sys:        Pulseq system struct (mr.opts)
%   CRT:        common raster time for GE compatibility (s)
%
% Duration is set large enough to support the maximum-area direction without
% exceeding the PNS limit (scale factor < 1 avoids max slew rate).

deltak = 1./fov;
tmp = 0.5; % scale factor <1 to avoid PNS
gxPre = trap4ge(mr.scaleGrad(mr.makeTrapezoid('x', sys, 'Area', -Nx/2*deltak(1)/tmp), tmp, sys), CRT, sys);
gyPre = trap4ge(mr.scaleGrad(mr.makeTrapezoid('y', sys, 'Area', -Ny/2*deltak(2)/tmp), tmp, sys), CRT, sys);
gzPre = trap4ge(mr.scaleGrad(mr.makeTrapezoid('z', sys, 'Area', -Nz/2*deltak(3)/tmp), tmp, sys), CRT, sys);
end
