# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Generating sequences

All sequence generation runs in MATLAB. The normal workflow:

```matlab
% From the ArbEPI/ directory in MATLAB:
main          % generates all 4 sequences at once
```

Or step by step:
```matlab
omegas = gen_sampling_masks(R);   % Ny × Nz × Nframes logical mask
ArbEPI(omegas);                   % writes output/ArbEPI.seq/.pge, output/samp_locs.mat
EPIcal();                         % writes output/EPIcal.seq/.pge, output/kxoe<Nx>.mat
GRE();                            % writes output/GRE.seq/.pge
noise();                          % writes output/noise.seq/.pge
```

`EPIcal` and `noise` must run after `ArbEPI` — they load `output/samp_locs.mat`.

There is no automated test suite. Each sequence function runs `seq.checkTiming()` internally and prints pass/fail.

## Architecture

### Data flow

```
params.m  ──►  gen_sampling_masks(R)  ──►  omegas (Ny×Nz×Nframes logical)
                                                │
                                                ▼
                                          ArbEPI(omegas)
                                            │  mask2epi() called per frame
                                            │  schedules: Nframes×Nshots×ETL×2
                                            │  saved to output/samp_locs.mat
                                            ▼
                                 EPIcal() / noise()   ← load samp_locs.mat
```

### Key design decisions

**`mask2epi`** (`lib/mask2epi.m`) is the core algorithm. It partitions a 2D `(Ny, Nz)` sampling mask into `Nshots` EPI trajectories, each of length `ETL`. Ordering constraints: samples near ky = 0 are spread center-out across shots; ky is non-decreasing within each echo train.

**`make_readout_grads`** (`lib/make_readout_grads.m`) returns a struct `rg` with pre-built gradient objects. Blips (`rg.gyBlip`, `rg.gzBlip`) are stored at *unit amplitude* and scaled at assembly time via `mr.scaleGrad(rg.gyBlip, step_size)`. The readout trapezoid (`rg.gro`) is circularly shifted so blips fit within each Pulseq block boundary. `rg.gro1` and `rg.gro2` are the leading/trailing half-trapezoids played outside and inside the echo loop respectively.

**`EPIcal`** mirrors `ArbEPI`'s gradient design exactly (same `rg`) but sets all blip scale factors to 0, so it acquires unencoded lines at k-space center for EPI ghost correction.

**`write_to_ge`** (`lib/write_to_ge.m`) converts a Pulseq `mr.Sequence` to GE TOPPE `.pge` format. It also runs the acoustic frequency check and opens a plot of the first 10 blocks.

### Toolchain

- **Pulseq** (`mr.*`) — sequence object, gradient/RF construction
- **TOPPE** (`toppe.*`, `trap4ge`) — GE-compatible gradient rounding to CRT
- **pge2** (`pge2.*`, `seq2ceq`) — convert Pulseq → GE `.pge` format

`trap4ge(grad, CRT, sys)` rounds all gradient timing to the common raster time (`CRT = 20e-6 s`). Every gradient in the sequence must pass through it before being added to blocks.

### Configuration

`params.m` is a script (not a function) that sets all parameters in the caller's workspace. All top-level functions call it as:

```matlab
projRoot = fileparts(mfilename('fullpath'));
addpath(projRoot);
addpath(fullfile(projRoot, 'lib'));
params;
```

Do **not** use `run('params.m')` inside a function — MATLAB documents that `run()` may not populate the calling function's workspace. Use the direct `params;` call instead.

Avoid naming workspace variables after MATLAB built-in functions. For example, the EPI flip angle is `fa` (not `alpha`, which shadows the built-in `alpha()` transparency function and causes a runtime error).

### Output files

All generated files go to `output/` (gitignored). Key files consumed by reconstruction:

- `output/samp_locs.mat` — `schedules` (Nframes×Nshots×ETL×2 ky/kz indices) and `parts` (Ny×Nz partition map)
- `output/kxoe<Nx>.mat` — odd/even echo k-space trajectories from EPIcal, used for ghost correction
