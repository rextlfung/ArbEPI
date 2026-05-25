function plot_kt_mask(data)
% PLOT_KT_MASK Interactive viewer for k-t sampling mask (ky x kz x time)
%   Scroll wheel, slider, or arrow keys to move through time

    % Ensure orientation matches expectation (x horizontal)
    data = permute(data, [2, 1, 3]);

    % Compute global intensity scaling
    dataMin = min(data(:));
    dataMax = max(data(:));

    % Create figure
    fig = figure('Name', '3D Time Viewer', 'NumberTitle', 'off', ...
                 'WindowScrollWheelFcn', @scrollCallback, ...
                 'WindowKeyPressFcn', @keyPressCallback);

    % Create axes — bottom margin for x-label/ticks, right margin for slider
    ax = axes('Parent', fig, ...
              'Position', [0.01 0.08 0.93 0.88]);
    ax.PositionConstraint = 'innerposition';

    % Frame number text (right side, above slider)
    frameText = uicontrol('Style', 'text', ...
                          'Parent', fig, ...
                          'Units', 'normalized', ...
                          'Position', [0.95 0.93 0.05 0.05], ...
                          'String', '', ...
                          'FontSize', 20, ...
                          'BackgroundColor', 'w', ...
                          'ForegroundColor', 'k', ...
                          'HorizontalAlignment', 'center');

    % Slider (right side)
    scrollBar = uicontrol('Style', 'slider', ...
                          'Parent', fig, ...
                          'Units', 'normalized', ...
                          'Position', [0.96 0.06 0.04 0.86], ...
                          'Min', 1, 'Max', size(data, 3), ...
                          'Value', 1, ...
                          'SliderStep', [1/(size(data,3)-1) , 10/(size(data,3)-1)], ...
                          'Callback', @sliderCallback);

    % Store handles
    handles.data = data;
    handles.ax = ax;
    handles.imgHandle = []; % placeholder for imagesc handle
    handles.currentTime = 1;
    handles.maxTime = size(data, 3);
    handles.frameText = frameText;
    handles.scrollBar = scrollBar;
    handles.dataMin = dataMin;
    handles.dataMax = dataMax;
    guidata(fig, handles);

    % Initial plot
    plotFrame(fig);
end

function scrollCallback(src, event)
    handles = guidata(src);
    if event.VerticalScrollCount > 0
        handles.currentTime = min(handles.currentTime + 1, handles.maxTime);
    else
        handles.currentTime = max(handles.currentTime - 1, 1);
    end
    set(handles.scrollBar, 'Value', handles.maxTime + 1 - handles.currentTime);
    guidata(src, handles);
    plotFrame(src);
end

function sliderCallback(src, ~)
    fig = ancestor(src, 'figure');
    handles = guidata(fig);
    handles.currentTime = round(handles.maxTime + 1 - get(src, 'Value'));
    guidata(fig, handles);
    plotFrame(fig);
end

function keyPressCallback(src, event)
    handles = guidata(src);
    switch event.Key
        case {'rightarrow','uparrow'}
            handles.currentTime = min(handles.currentTime + 1, handles.maxTime);
        case {'leftarrow','downarrow'}
            handles.currentTime = max(handles.currentTime - 1, 1);
        otherwise
            return;
    end
    set(handles.scrollBar, 'Value', handles.maxTime + 1 - handles.currentTime);
    guidata(src, handles);
    plotFrame(src);
end

function plotFrame(fig)
    handles = guidata(fig);

    % Get 2D frame at current time
    img = handles.data(:, :, handles.currentTime);

    % Normalize consistently across time
    img = (img - handles.dataMin) / (handles.dataMax - handles.dataMin);

    if isempty(handles.imgHandle) || ~isvalid(handles.imgHandle)
        Nz = size(handles.data, 1);
        Ny = size(handles.data, 2);
        handles.imgHandle = imagesc(handles.ax, [-Ny/2, Ny/2-1], [-Nz/2, Nz/2-1], img);
        axis(handles.ax, 'image');
        colormap(handles.ax, gray);
        set(handles.ax, 'FontSize', 20);
        xlabel(handles.ax, 'k_y', 'FontSize', 20);
        ylabel(handles.ax, 'k_z', 'FontSize', 20);
        title(handles.ax, 'k-t sampling mask', 'FontSize', 20);
    else
        set(handles.imgHandle, 'CData', img);
    end

    % Update frame text
    set(handles.frameText, 'String', ...
        sprintf('t = %d / %d', handles.currentTime, handles.maxTime));

    guidata(fig, handles);
    drawnow;
end
