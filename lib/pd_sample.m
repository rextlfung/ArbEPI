function mask = pd_sample(img_shape, accel, varargin)
% PD_SAMPLE Generate variable-density Poisson-disc sampling pattern with EXACT acceleration.
%
%   mask = PD_SAMPLE(img_shape, accel, ...) generates a mask where the 
%   total number of samples is EXACTLY floor(total_pixels / accel).
%   It achieves this by slightly oversampling and then randomly pruning 
%   excess points.
%
%   Inputs:
%       img_shape (1x2 vector): Image dimensions [ny, nx].
%       accel (scalar): Target acceleration factor (must be > 1).
%
%   Optional Name-Value Arguments:
%       'calib'        : (1x2 vector) Calibration region size [calib_y, calib_x]. Default [0, 0].
%       'dtype'        : (string) Output data type ('logical', 'double', 'complex'). Default 'logical'.
%       'crop_corner'  : (logical) Toggle whether to crop sampling corners. Default true.
%       'seed'         : (scalar) Random seed. Default 0.
%       'max_attempts' : (scalar) Max attempts to generate point. Default 30.
%       'tol'          : (scalar) Tolerance for binary search loop. Default 0.1.
%
%   Example:
%       mask = poisson([256, 256], 4, 'calib', [24, 24]);
%       fprintf('Actual Accel: %.4f\n', prod([256, 256]) / sum(mask(:)));

    % --- 1. Parse Inputs ---
    p = inputParser;
    p.addRequired('img_shape', @(x) isnumeric(x) && numel(x)==2);
    p.addRequired('accel', @(x) isnumeric(x) && isscalar(x));
    p.addParameter('calib', [0, 0], @(x) isnumeric(x) && numel(x)==2);
    p.addParameter('dtype', 'logical', @(x) ischar(x) || isstring(x));
    p.addParameter('crop_corner', true, @islogical);
    p.addParameter('seed', 0, @isnumeric);
    p.addParameter('max_attempts', 30, @isnumeric);
    p.addParameter('tol', 0.1, @isnumeric);
    p.parse(img_shape, accel, varargin{:});

    calib = p.Results.calib;
    dtype = p.Results.dtype;
    crop_corner = p.Results.crop_corner;
    seed = p.Results.seed;
    max_attempts = p.Results.max_attempts;
    tol = p.Results.tol;

    if accel <= 1
        error('accel must be greater than 1, got %f', accel);
    end

    if seed ~= 0
        s_rng = rng; % Save random state
        rng(seed);
    end

    % --- 2. Calculate Exact Target Count ---
    ny = img_shape(1);
    nx = img_shape(2);
    total_pixels = nx * ny;
    target_samples = floor(total_pixels / accel);

    % --- 3. Setup Coordinates & Radius ---
    [Y, X] = ndgrid(0:ny-1, 0:nx-1);

    x_dist = max(abs(X - nx / 2) - calib(2) / 2, 0);
    x_dist = x_dist / max(x_dist(:)); 
    
    y_dist = max(abs(Y - ny / 2) - calib(1) / 2, 0);
    y_dist = y_dist / max(y_dist(:)); 
    
    r = sqrt(x_dist.^2 + y_dist.^2);

    % --- 4. Binary Search for Slope ---
    % TRICK: We target a slightly HIGHER density (lower accel) to ensure we 
    % generate enough points to prune down later.
    accel_search = accel * 0.95; 

    slope_max = max(nx, ny);
    slope_min = 0;
    
    mask = [];
    
    % Limit iterations to prevent infinite loops
    for i = 1:50 
        slope = (slope_max + slope_min) / 2;
        
        radius_x = max((1 + r * slope) * nx / max(nx, ny), 1);
        radius_y = max((1 + r * slope) * ny / max(nx, ny), 1);
        
        mask = poisson_disc_core(nx, ny, max_attempts, radius_x, radius_y, calib, 0); % Use 0 seed inside loop
        
        if crop_corner
            mask = mask .* (r <= 1);
        end
        
        num_samples = sum(mask(:));
        current_accel = total_pixels / num_samples;
        
        % Check if we are close enough OR if we have bracketed the target
        if abs(current_accel - accel_search) < tol
            break;
        end
        
        if current_accel < accel_search
            slope_min = slope; % Increase slope to reduce samples
        else
            slope_max = slope; % Decrease slope to increase samples
        end
    end

    % --- 5. ENFORCE EXACT ACCELERATION (Pruning/Filling) ---
    current_samples = sum(mask(:));
    
    % Identify Calibration Region (Protected from pruning)
    calib_mask = zeros(ny, nx);
    y_start = max(1, floor(ny/2 - calib(1)/2) + 1);
    y_end   = min(ny, floor(ny/2 + calib(1)/2));
    x_start = max(1, floor(nx/2 - calib(2)/2) + 1);
    x_end   = min(nx, floor(nx/2 + calib(2)/2));
    calib_mask(y_start:y_end, x_start:x_end) = 1;

    if current_samples > target_samples
        % Too many samples: Prune randomly (excluding calibration)
        num_to_remove = current_samples - target_samples;
        
        % Find candidates (points that are in mask but NOT in calib region)
        candidates = find(mask == 1 & calib_mask == 0);
        
        if ~isempty(candidates)
            % Random shuffle
            p_idx = randperm(length(candidates));
            % Select points to remove
            remove_idx = candidates(p_idx(1:min(num_to_remove, end)));
            mask(remove_idx) = 0;
        end
        
    elseif current_samples < target_samples
        % Too few samples: Add randomly (to empty spots)
        % This is rare due to the 0.95 factor above, but handled for robustness.
        num_to_add = target_samples - current_samples;
        
        % Find empty spots
        candidates = find(mask == 0);
        
        if ~isempty(candidates)
            p_idx = randperm(length(candidates));
            add_idx = candidates(p_idx(1:min(num_to_add, end)));
            mask(add_idx) = 1;
        end
    end

    % --- 6. Final Output ---
    switch lower(dtype)
        case 'complex', mask = complex(mask);
        case 'double',  mask = double(mask);
        case 'logical', mask = logical(mask);
    end

    if seed ~= 0
        rng(s_rng);
    end
end

% -------------------------------------------------------------------------
% LOCAL HELPER FUNCTION
% -------------------------------------------------------------------------
function mask = poisson_disc_core(nx, ny, max_attempts, radius_x, radius_y, calib, seed)
    mask = zeros(ny, nx);
    if seed ~= 0, rng(seed); end

    % Calibration Region
    calib_y = calib(1); calib_x = calib(2);
    y_start = max(1, floor(ny/2 - calib_y/2) + 1);
    y_end   = min(ny, floor(ny/2 + calib_y/2));
    x_start = max(1, floor(nx/2 - calib_x/2) + 1);
    x_end   = min(nx, floor(nx/2 + calib_x/2));
    mask(y_start:y_end, x_start:x_end) = 1;

    % Active List
    pxs = zeros(1, nx * ny); pys = zeros(1, nx * ny);
    pxs(1) = randi([0, nx-1]); pys(1) = randi([0, ny-1]);
    num_actives = 1;
    
    while num_actives > 0
        idx = randi([1, num_actives]);
        px = floor(pxs(idx)); py = floor(pys(idx));
        rx = radius_x(py+1, px+1); ry = radius_y(py+1, px+1);
        
        done = false; k = 0;
        while ~done && k < max_attempts
            v = sqrt(rand()*3 + 1); t = 2*pi*rand();
            qx = px + v*rx*cos(t); qy = py + v*ry*sin(t);
            
            if qx >= 0 && qx < nx && qy >= 0 && qy < ny
                sx = max(floor(qx - rx), 0); ex = min(ceil(qx + rx), nx - 1);
                sy = max(floor(qy - ry), 0); ey = min(ceil(qy + ry), ny - 1);
                
                valid = true;
                for x = sx:ex
                    for y = sy:ey
                        if mask(y+1, x+1) == 1
                            if ((qx-x)/radius_x(y+1, x+1))^2 + ((qy-y)/radius_y(y+1, x+1))^2 < 1
                                valid = false; break;
                            end
                        end
                    end
                    if ~valid, break; end
                end
                
                if valid
                    done = true; num_actives = num_actives + 1;
                    pxs(num_actives) = qx; pys(num_actives) = qy;
                    mask(floor(qy)+1, floor(qx)+1) = 1;
                end
            end
            k = k+1;
        end
        if ~done
            pxs(idx) = pxs(num_actives); pys(idx) = pys(num_actives);
            num_actives = num_actives - 1;
        end
    end
end