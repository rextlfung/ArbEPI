# ArbEPI: matlab/pulseq code for generating 3D-EPI sequences from arbitrary 2D sampling masks in the phase-encode-partition (ky-kz) plane
![An example of a Poisson-Disc 3D-EPI trajectory and corresponding PSF.](poisson3DEPI.png)

## Getting started
1. Set experimental parameters in ```params.m```. I recommend making a copy of this file for each experiment for reconstruction afte the scan.
2. Run ```ArbEPI.m``` to generate ```ArbEPI.seq``` and ```ArbEPI.pge```. Plug in your custom 2D sampling mask (```omega```) the section at line 44. TODO: move sampling mask selection to somewhere that's not them middle of the script.  
    a. Also generates ```samp_locs.mat```, which saves the sampling schedule for use during reconstruction.
3. Run ```EPIcal.m``` to generate ```EPIcal.seq``` and ```EPIcal.pge```. This sequence contains dummy excitations and blipless EPI trains for: 1) setting receiver gain and 2) odd/even and linear phase correction.  
    a. Also generates ```kxoe$Nx.mat```, which saves the odd/even locations for use during reconstruction.
4. (optional) Run ```GRE.m``` to generate ```GRE.seq``` and ```GRE.pge```. Acquires gradient echo images for SENSE reconstruction.

## On the scanner (GE)
1. Set up ```.pge``` files by following: https://github.com/jfnielsen/TOPPEpsdSourceCode/tree/UserGuide/v7#running-a-sequence-on-the-scanner (private repo).
2. Run a localizer sequence.
3. Run ```EPIcal.pge``` using auto-prescan.
4. Run ```ArbEPI.pge``` using manual-prescan, keeping the same receiver gain settings from the last step.
5. (optional) Run ```GRE.pge``` using auto-prescan.

## Methods overview
1. The core logic is contained in ```mask2epi()```, a function that accepts any 2D sampling mask of size (Ny, Nz) and partitions it into efficient 3D-EPI trajectories. The number of shots and echo train length must also be specified.
2. Some example 2D sampling mask generators are provided:  
    a. ```caipi_sample()``` creates uniformly spaced, CAIPI-shifted sampling masks.  
    b. ```pd_sample()``` creates Poisson-disc sampling masks based on the SigPy implementation [1].  
    c. ```rand_sample()``` creates random sampling masks according to provided 2D probability mass functions. As an example, ```gen_gaussian_pdf()``` is provided to generate a 2D Gaussian pmf.  
3. Tajectory generation details:  
    a. The number of samples about ky = 0 are approximately even so that the echo time, defined as time of crossing ky = 0, is approximately constant across shots. This is achieved by collecting samples in an in-to-out order when constructing each echo train.  
    b. ky is nondecreasing within each echo train, so that B0-induced distortion along the phase encoding direction stays coherent.  
4. Pulseq details:  
    a. All ky and kz encoding blips have the same duration for efficient representation.  
    b. All ky and kz encoding blips are contained within a block so that they can be written as the same waveform with scaled amplitudes for efficient representation. This is achieved by splicing the trapezoidal readout gradient into two segments, then appending a sign-flipped version of the first segment to the second segment, creating a new readout gradient as an arbitrary waveform.  
5. Unsolved problems:  
    a. Acqusition speed is bottlenecked by large jumps in k-space, so we need a way to constrain large jumps without compromising the randomness of sampling.  

## Why I think this works
More random sampling patterns  
--> More spatiotemporally incoherent aliases that are noise-like in the singular value domain  
--> Removable via low-rank regularization and/or other denoising methods.

## References
1. Ong, F., & Lustig, M. (2019, May). SigPy: a python package for high performance iterative reconstruction. In Proceedings of the ISMRM 27th Annual Meeting, Montreal, Quebec, Canada (Vol. 4819, No. 5).
