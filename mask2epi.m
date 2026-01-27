%% mask2epi.m 
% Input arguments:
% mask = 2D Cartesian sampling mask (ky x kz)
% ETL = echo train legnth = number of samples per shot
% Nshots = number of shots/excitations per volume
%
% Output arguments:
% schedule = schedule of (ky, kz) locations to sample
% parts = partitionining of samples
%
function [schedule, parts] = mask2epi(mask, ETL, Nshots)
    assert(sum(mask, 'all') == Nshots*ETL, ...
           'Number of samples must == Nshots*ETL')

    schedule = zeros(Nshots, ETL, 2); % preallocate output array
    
    % first pass, partition samples
    [Ny, Nz] = size(mask);
    parts = zeros(Ny, Nz);
    echo_count = 0; part = 1;
    for iz = 1:Nz
        for iy = center_out(1:Ny)
            if mask(iy,iz)
                assert(part <= Nshots, ...
                           'Partitions exceeds number of shots.')
                parts(iy,iz) = part;
                echo_count = echo_count + 1;
                if echo_count == ETL
                    part = part + 1; % move to the next partition
                    echo_count = 0;  % reset echo count
                end
            end
        end
    end

    % second pass, infer optimal schedul subject to
    % 1. non-decreasing ky
    % 2. minimum distance traveled
    for shot = 1:Nshots
        samples = (parts == shot);
        [kz_indices, ky_indices] = find(samples.');
        kz_min = min(kz_indices);
        kz_max = max(kz_indices);

        echo_count = 1;
        for iy = unique(ky_indices).'
            if echo_count > 1
                last_kz = schedule(shot, echo_count - 1, 2);
                if abs(last_kz - kz_min) < abs(last_kz - kz_max)
                    kz_range = kz_min:kz_max;
                else
                    kz_range = kz_max:-1:kz_min;
                end
            else
                kz_range = kz_min:kz_max;
            end
            
            for iz = kz_range
                if samples(iy,iz)
                    schedule(shot, echo_count, 1) = iy;
                    schedule(shot, echo_count, 2) = iz;
                    echo_count = echo_count + 1;
                end
            end
        end
    end
end

% helper function for iterating through a sequence with center-out order
function out = center_out(in)
    N = length(in);
    out = zeros(1,N);
    tmp = flip(fftshift(in));
    j = 1;
    for i = 1:ceil(N/2)
        out(j) = tmp(i);
        j = j + 1;
        if j > N, break, end
        out(j) = tmp(end-i+1);
        j = j + 1;
    end
end