function hIm = imdisp(I, varargin)
%IMDISP  Display one or more images nicely
%
% Examples:
%   imdisp
%   imdisp(I)
%   imdisp(I, map)
%   imdisp(I, lims)
%   imdisp(I, map, lims)
%   imdisp(..., param1, value1, param2, value2, ...)
%   h = imdisp(...)
%
% This function displays one or more images nicely. Images can be defined
% by arrays or filenames. Multiple images can be input in a cell array or
% stacked along the fourth dimension, and are displayed as a grid of
% subplots (an improvement over MONTAGE). The size of grid is calculated or
% user defined. The figure size is set so that images are magnified by an
% integer value.
%
% If the image grid size is user defined, images not fitting in the grid
% can be scrolled through using the following key presses:
%    Up - Back a row.
%    Down - Forward a row.
%    Left - Back a page (or column if there is only one row).
%    Right - Forward a page (or column if there is only one row).
%    Shift - 2 x speed.
%    Ctrl - 4 x speed.
%    Shift + Ctrl - 8 x speed.
%
% This allows fast scrolling through a movie or image stack, e.g. 
%    imdisp(imstack, 'Size', 1)
% The function can be used as a visual DIR, e.g. 
%    imdisp()
% to display all images in the current directory on a grid, or 
%    imdisp({}, 'Size', 1)
% to scroll through them one at a time.
%
% IN:
%   I - MxNxCxP array of images, or 1xP cell array. C is 1 for indexed
%       images or 3 for RGB images. P is the number of images. If I is a
%       cell array then each cell must contain an image. Images can equally
%       be defined by filenames. If I is an empty cell array then all the
%       images in the current directory are used. Default: {}.
%   map - Kx3 colormap to be used with indexed images. Default: gray(256).
%   lims - [LOW HIGH] display range for indexed images. Default: [min(I(:))
%          max(I(:))].
%   Optional parameters - name, value parameter pairs for the following:
%      'Size' - [H W] size of grid to display image on. If only H is given
%               then W = H. If either H or W is NaN then the number of rows
%               or columns is chosen such that all images fit. If both H
%               and W are NaN or the array is empty then the size of grid
%               is chosen to fit all images in as large as possible.
%               Default: [].
%      'Indices' - 1xL list of indices of images to display. Default: 1:P.
%      'Border' - [TB LR] borders to give each image top and bottom (TB)
%                 and left and right (LR), to space out images. Borders are
%                 normalized to the subplot size, i.e. TB = 0.01 gives a
%                 border 1% of the height of each subplot. If only TB is
%                 given, LR = TB. Default: 0.01.
%      'DisplayRange' - Same as lims input.
%      'Map' - Kx3 colormap or (additionally from above) name of MATLAB
%              colormap, for use with indexed images. Default: gray(256).
%
% OUT:
%   h - HxW array of handles to images.
%
%   See also IMAGE, IMAGESC, IMSHOW, MONTAGE.

% Parse inputs
[map layout gap indices lims] = parse_inputs(varargin);

if nargin == 0 || (iscell(I) && isempty(I))
    % Read in all the images in the directory
    I = get_im_names;
    if isempty(I)
        % No images found
        if nargout > 0
            hIm = [];
        end
        return
    end
end

% Check if input is filenames
if ischar(I)
    [x y c] = size(I);
    if (x > 1 && y > 1) || c > 1 
        I = num2cell(I, 2);
    else
        I = {I(:)'};
    end
end

% Get limits, etc.
if isnumeric(I) || islogical(I)
    [y x c n] = size(I);
    if isempty(lims)
        lims = min_max(I);
    elseif isequal(0, lims)
        lims = default_limits(I);
    elseif c == 3
        % Rescale
        if ~isfloat(I)
            I = single(I);
        end
        I = min(max((I - lims(1)) ./ (lims(2) - lims(1)), 0), 1);
    end
    if isfloat(I) && c == 3 && n > 1
        I = uint8(I * 256 - 0.5);
        lims = round(lims * 256 - 0.5);
    end
elseif iscell(I)
    n = numel(I);
    A = I{1};
    if ischar(A)
        % Read in the image (or images for multi-frame files)
        if n == 1
            I = imread_rgb_multi(A);
            if iscell(I)
                n = numel(I);
                A = I{1};
                [y x c] = size(A);
            else
                [y x c n] = size(I);
                A = I;
            end
        else
            A = imread_rgb(A);
            I{1} = A;
            [y x c] = size(A);
        end
    else
        [y x c] = size(A);
    end
    % Assume all images are the same size and type as the first
    if isempty(lims) || isequal(0, lims)
        lims = default_limits(A);
    end
else
    error('I not of recognized type.');
end

% Select indexed images
if ~isequal(indices, -1)
    if iscell(I)
        I = I(indices);
        n = numel(I);
    else
        I = I(:,:,:,indices);
        n = size(I, 4);
    end
end

% Get the current figure
hFig = get(0, 'CurrentFigure');
if isempty(hFig)
    % Create a new figure
    hFig = figure;
elseif n > 1
    % Clear the figure
    hFig = clf(hFig, 'reset');
end

% Set the colormap
set(hFig, 'Colormap', map);

% Display the image(s)
if n == 0
    hIm = display_image([], gca, [0 1]);
    
    if nargout == 0
        clear hIm % Avoid printing this out
    end
    return
elseif n == 1
    % IMSHOW mode
    % Display the single image
    hAx = gca;
    if iscell(I)
        I = I{1};
    end
    hIm = display_image(I, hAx, lims);
    
    if nargout == 0
        clear hIm % Avoid printing this out
    end
    
    % Only resize image if it is alone in the figure
    if numel(findobj(get(hFig, 'Children'), 'Type', 'axes')) > 1
        return
    end
    % Could still be the first subplot - do another check
    axesPos = get(hAx, 'Position');
    newAxesPos = [gap(1) gap(end) 1-2*gap(1) 1-2*gap(end)];
    if isequal(axesPos, get(hFig, 'DefaultAxesPosition'))
        % Default position => not a subplot
        % Fill the window
        set(hAx, 'Units', 'normalized', 'Position', newAxesPos);
        axesPos = newAxesPos;
    end
    if ~isequal(axesPos, newAxesPos)
        % Figure not alone, so don't resize.
        return
    end
    layout = [1 1];
else
    % MONTAGE mode
    % Compute a good layout
    layout = choose_layout(n, y, x, layout);

    % Create a data structure to store the data in
    num = prod(layout);
    state.num = num * ceil(n / num);
    hIm = zeros(layout);
    hAx = zeros(layout);

    % Set the first lot of images
    index = mod(0:num-1, state.num) + 1;
    hw = 1 ./ layout;
    gap = gap ./ layout;
    dims = hw - 2 * gap;
    dims = dims([2 1]);
    for a = 1:layout(1)
        for b = 1:layout(2)
            c = index(b + (layout(1) - a) * layout(2));
            if c > n
                A = [];
            elseif iscell(I)
                A = I{c};
                if ischar(A)
                    A = imread_rgb(A);
                    I{c} = A;
                end
            else
                A = I(:,:,:,c);
            end
            hAx(a,b) = axes('Position', [(b-1)*hw(2)+gap(2) (a-1)*hw(1)+gap(1) dims], 'Units', 'normalized');
            hIm(a,b) = display_image(A, hAx(a,b), lims);
        end
    end
    
    % Check if we need to be able to scroll through images
    if n > num
        % Intialize rest of data structure
        state.hIm = hIm;
        state.hAx = hAx;
        state.index = 1;
        state.layout = layout;
        state.n = n;
        state.I = I;
        % Set the callback for image navigation, and save the image data in the figure
        set(hFig, 'KeyPressFcn', @keypress_callback, 'Interruptible', 'off', 'BusyAction', 'cancel', 'UserData', state);
    end
    
    % Flip hIm so it matches the layout
    hIm = hIm(end:-1:1,:);
    
    if nargout == 0
        clear hIm % Avoid printing this out
    end
end

if strcmp(get(hFig, 'WindowStyle'), 'docked')
    % Figure is docked, so can't resize
    return
end

% Set the figure size well
% Compute the image size
ImSz = layout([2 1]) .* [x y] ./ (1 - 2 * gap([end 1]));
    
% Get the size of the monitor we're on
figPosCur = get(hFig, 'Position');
% Monitor sizes
MonSz = get(0, 'MonitorPositions');
MonOn = size(MonSz, 1);
if MonOn > 1
    % Make the origin the top left corner of the primary monitor
    correction = 0;
    if ispc
        for a = 1:MonOn
            if isequal(MonSz(a,1:2), [1 1])
                correction = MonSz(a,4);
                break
            end
        end
    end
    % Determine which monitor the centre of the image is on
    figCenter = figPosCur(1:2) + figPosCur(3:4) / 2;
    figCenter = MonSz - repmat(figCenter, [MonOn 2]);
    MonOn = all(sign(figCenter) == repmat([-1 -1 1 1], [MonOn 1]), 2);
    MonOn(1) = MonOn(1) | ~any(MonOn);
    MonSz = MonSz(MonOn,:);
    % Correct the size
    MonSz(3:4) = MonSz(3:4) - MonSz(1:2) + 1;
    % Correct the origin
    if correction
        MonSz(2) = correction - MonSz(4) - MonSz(2) + 2;
    end
end

% Check if the window is maximized
% This is a hack which may only work on Windows! No matter, though.
if isequal(MonSz([1 3]), figPosCur([1 3]))
    % Leave maximized
    return
end

% Compute the size to set the window
MaxSz = MonSz(3:4) - [20 120];
RescaleFactor = min(MaxSz ./ ImSz);
if RescaleFactor > 1
    % Integer scale for enlarging, but don't make too big
    MaxSz = min(MaxSz, [1200 800]);
    RescaleFactor = max(floor(min(MaxSz ./ ImSz)), 1);
end
figPosNew = ceil(ImSz * RescaleFactor);

% Don't move the figure if the size isn't changing
if isequal(figPosCur(3:4), figPosNew)
    return
end

% Keep the centre of the figure stationary
figPosNew = [floor(figPosCur(1:2)+(figPosCur(3:4)-figPosNew)/2) figPosNew];

% Ensure the figure is in bounds
figPosNew(1:2) = min(max(figPosNew(1:2), MonSz(1:2)+6), MonSz(1:2)+MonSz(3:4)-[6 101]-figPosNew(3:4));

% Set the figure size and position
set(hFig, 'Position', figPosNew);
return

%% Keypress callback
% The function which does all the display stuff
function keypress_callback(fig, event_data)
% Check what key was pressed and update the image index as necessary
switch event_data.Character
    case 28 % Left
        up = -1; % Back a page
    case 29 % Right
        up = 1; % Forward a page
    case 30 % Up
        up = -0.1; % Back a row
    case 31 % Down
        up = 0.1; % Forward a row
    otherwise
        % Another key was pressed - ignore it
        return
end
% Use control and shift for faster scrolling
if ~isempty(event_data.Modifier)
    up = up * (2 ^ (strcmpi(event_data.Modifier, {'shift', 'control'}) * [1; 2]));
end
% Get the state data
state = get(fig, 'UserData');
% Get the current index
index = state.index;
% Get number of images
n = prod(state.layout);
% Generate 12 valid indices
if abs(up) < 1
    % Increment by row
    index = index + state.layout(2) * (up * 10) - 1;
else
    if state.layout(1) == 1
        % Increment by column
        index = index + up - 1;
    else
        % Increment by page
        index = index + n * up - 1;
    end
end
index = mod(index:index+n, state.num) + 1;
% Plot the images
figure(fig);
for a = 1:state.layout(1)
    for b = 1:state.layout(2)
        % Get the image
        c = index(b + (state.layout(1) - a) * state.layout(2));
        if c > state.n
            % Set the image data
            set(state.hIm(a,b), 'CData', []);
        elseif iscell(state.I)
            A = state.I{c};
            if ischar(A)
                % Filename - read the image from disk
                A = imread_rgb(A);
                state.I{c} = A;
            end
            % Set the image data
            set(state.hIm(a,b), 'CData', A);
            % Reset the axes limits
            if ~isempty(A)
                set(state.hAx(a,b), 'XLim', [0.5 size(A, 2)+0.5], 'YLim', [0.5 size(A, 1)+0.5]);
            end
        else
            % Set the image data
            set(state.hIm(a,b), 'CData', state.I(:,:,:,c));
        end
    end
end
drawnow;
% Save the current index
state.index = index(1);
set(fig, 'UserData', state);
return

%% Display the image
function hIm = display_image(A, hAx, lims)
if isempty(A)
    hIm = image(zeros(1, 1, 3));
    set(hIm, 'CData', []);
else
    hIm = image(A);
end
set(hAx, 'Visible', 'off', 'DataAspectRatio', [1 1 1], 'DrawMode', 'fast', 'CLim', lims);
set(get(hAx, 'XLabel'), 'Visible', 'on');
set(get(hAx, 'YLabel'), 'Visible', 'on');
set(get(hAx, 'Title'), 'Visible', 'on');
set(hIm, 'CDataMapping', 'scaled');
return

%% Choose a good layout for the images
function layout = choose_layout(n, y, x, layout)
v = numel(layout);
N = isnan(layout);
if v == 0 || all(N)
    % Compute approximate layout
    sz = get(0, 'ScreenSize');
    sz = sz(3:4) ./ [x y];
    layout = ceil(sz([2 1]) ./ sqrt(prod(sz) / n));
    % Remove superfluous rows or columns
    while 1
        switch ([prod(layout - [1 0]) prod(layout - [0 1])] >= n) * [2; 1]
            case 0
                break;
            case 1
                layout = layout - [0 1];
            case 2
                layout = layout - [1 0];
            case 3
                if min(sz .* (layout - [0 1])) > min(sz .* (layout - [1 0]))
                    layout = layout - [0 1];
                else
                    layout = layout - [1 0];
                end
        end
    end
elseif v == 1
    layout = layout([1 1]);
elseif any(N)
    layout(N) = ceil(n / layout(~N));
end
layout = reshape(layout, 1, 2);
return

%% Read image to uint8 rgb array
function A = imread_rgb(name)
try
    [A map alpha] = imread(name);
catch
    % Format not recognized by imread, so create a red cross (along diagonals)
    A = eye(101) | diag(ones(100, 1), 1) | diag(ones(100, 1), -1);
    A = (uint8(1) - uint8(A | flipud(A))) * uint8(255);
    A = cat(3, zeros(size(A), 'uint8')+uint8(255), A, A);
    return
end
A = A(:,:,:,1); % Keep only first frame of multi-frame files
if ~isempty(map)
    map = uint8(map * 256 - 0.5); % Convert to uint8 for storage
    A = reshape(map(uint32(A)+1,:), [size(A) size(map, 2)]); % Assume indexed from 0
elseif size(A, 3) == 4
    if lower(name(end)) == 'f'
        % TIFF in CMYK colourspace - convert to RGB
        if isfloat(A)
            A = A * 255;
        else
            A = single(A);
        end
        A = 255 - A;
        A(:,:,4) = A(:,:,4) / 255;
        A = uint8(A(:,:,1:3) .* A(:,:,[4 4 4]));
    else
        % Assume 4th channel is an alpha matte
        alpha = A(:,:,4);
        A = A(:,:,1:3);
    end
end
if ~isempty(alpha)
    % Apply transprency over a grey checkerboard pattern
    if isa(alpha, 'uint8')
        alpha = double(alpha) / 255;
    end
    A = double(A) .* alpha(:,:,ones(1, size(A, 3)));
    sqSz = max(size(alpha));
    sqSz = floor(max(log(sqSz / 100), 0) * 10 + 1 + min(sqSz, 100) / 20);
    grid = repmat(85, ceil(size(alpha) / sqSz));
    grid(2:2:end,1:2:end) = 171;
    grid(1:2:end,2:2:end) = 171;
    grid = kron(grid, ones(sqSz));
    alpha = grid(1:size(A, 1),1:size(A, 2)) .* (1 - alpha);
    A = uint8(A + alpha(:,:,ones(1, size(A, 3))));
end
return

%% Read (potentially) multi-frame image to uint8 rgb array
function A = imread_rgb_multi(name)
try
    % Get file info
    info = imfinfo(name);
catch
    % Revert to standard case
    A = imread_rgb(name);
    return
end
if numel(info) < 2
    % Single image
    A = imread_rgb(name);
    return
else
    % Multi-frame image
    switch lower(info(1).Format)
        case 'gif'
            [A map] = imread(name, 'frames', 'all');
            if ~isempty(map)
                map = uint8(map * 256 - 0.5); % Convert to uint8 for storage
                A = reshape(map(uint32(A)+1,:), [size(A) size(map, 2)]); % Assume indexed from 0
                A = permute(A, [1 2 5 4 3]);
            end
        case {'tif', 'tiff'}
            A = cell(numel(info), 1);
            for a = 1:numel(A)
                [A{a} map] = imread(name, 'Index', a, 'Info', info);
                if ~isempty(map)
                    map = uint8(map * 256 - 0.5); % Convert to uint8 for storage
                    A{a} = reshape(map(uint32(A{a})+1,:), [size(A) size(map, 2)]); % Assume indexed from 0
                end
                if size(A{a}, 3) == 4
                    % TIFF in CMYK colourspace - convert to RGB
                    if isfloat(A{a})
                        A{a} = A{a} * 255;
                    else
                        A{a} = single(A{a});
                    end
                    A{a} = 255 - A{a};
                    A{a}(:,:,4) = A{a}(:,:,4) / 255;
                    A{a} = uint8(A(:,:,1:3) .* A{a}(:,:,[4 4 4]));
                end
            end
        otherwise
            % Multi-frame not supported for this format
            A = imread_rgb(name);
            return
    end
end
return

%% Get the names of all images in a directory
function L = get_im_names
D = dir;
n = 0;
L = cell(size(D));
% Go through the directory list
for a = 1:numel(D)
    % Check if file is a supported image type
    if numel(D(a).name) > 4 && ~D(a).isdir && (any(strcmpi(D(a).name(end-3:end), {'.png', '.tif', '.jpg', '.bmp', '.ppm', '.pgm', '.pbm', '.gif', '.ras'})) || any(strcmpi(D(a).name(end-4:end), {'.tiff', '.jpeg'})))
        n = n + 1;
        L{n} = D(a).name;
    end
end
L = L(1:n);
return

%% Parse inputs
function [map layout gap indices lims] = parse_inputs(inputs)

% Set defaults
map = [];
layout = [];
gap = 0;
indices = -1;
lims = 0;

% Check for map and display range
for b = 1:numel(inputs)
    if ~isnumeric(inputs{b})
        b = b - 1;
        break;
    end
    if size(inputs{b}, 2) == 3
        map = inputs{b};
    elseif numel(inputs{b}) < 3
        lims = inputs{b};
    end
end

% Go through option pairs
for a = b+1:2:numel(inputs)
    switch lower(inputs{a})
        case 'map'
            map = inputs{a+1};
            if ischar(map)
                map = feval(map, 256);
            end
        case {'size', 'grid'}
            layout = inputs{a+1};
        case {'gap', 'border'}
            gap = inputs{a+1};
        case 'indices'
            indices = inputs{a+1};
        case {'lims', 'displayrange'}
            lims = inputs{a+1};
        otherwise
            error('Input option %s not recognized', inputs{a});
    end
end

if isempty(map)
   map = gray(256);
end
return

%% Return default limits for the image type
function lims = default_limits(A)
if size(A, 3) == 1
    lims = min_max(A);
else
    lims = [0 1];
    if ~isfloat(A)
        lims = lims * double(intmax(class(A)));
    end
end
return

%% Return minimum and maximum values
function lims = min_max(A)
M = isfinite(A);
lims = double([min(A(M)) max(A(M))]);
if isempty(lims)
    lims = [0 1];
elseif lims(1) == lims(2)
    lims(2) = lims(1) + 1;
end
return
