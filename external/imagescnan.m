function [H,HNAN] = imagescnan(varargin)
%IMAGESCNAN   Scale data and display as image with uncolored NaNs.
%
%   SYNTAX:
%                imagescnan(U)
%                imagescnan(U,...,'NanColor',CNAN)
%                imagescnan(U,...,'NanMask',MNAN)
%                imagescnan(U,...,IOPT)
%                imagescnan(X,Y,U,...)
%     [H,HNAN] = imagescnan(...);
%
%   INPUT:
%     U    - 2 dimensional N-by-M image or N-by-M-by-3 RGB image.
%     X    - 2 extrema X-axis data; or the M values; or the N-by-M values
%            as obtained from MESHGRID (see DESCRIPTION below). 
%            DEFAULT: [1 N]
%     Y    - 2 extrema X-axis data; or the N values; or the N-by-M values
%            as obtained from MESHGRID (see DESCRIPTION below). 
%            DEFAULT: [1 M]
%     CNAN - Color for the NaNs elements. May be a char specifier or an [R
%            G B] triplet specifying the color.
%            DEFAULT: invisible (axes background color)
%     MNAN - Elements to be ignored besides not finite values. May be an
%            scalar or a logical M-by-N matrix indicating the elements to
%            be ignored.
%            DEFAULT: []
%     IOPT - IMAGE function normal optional pair arguments like
%            ('Parent',H) or/and CLIM like optional last argument as in
%            IMAGESC. 
%            DEFAULT: none
%
%   OUTPUT (all optional):
%     H    - Image handle
%     HNAN - Handle of every ignored (NaN) value colored patch.
%
%   DESCRIPTION:
%     MATLAB function IMAGESC does not work properly with NaNs. This
%     programs deals with this problem by including colored patches over
%     this elements and maybe others specyfied by the user with MNAN. 
%
%     Besides, those functions does not work properly with X,Y values
%     variable interval, but this functions does it by generating a whole
%     new image of several rectangular patches, but whose centers may not
%     lay in the specified coordinate (see NOTE below). This functionality
%     is experimental and not recommended (see ADDITIONAL NOTES inside this
%     program).
%
%     In previous release, 2-dim input images were transformed into a
%     3-dim RGB image. This is not used anymore (see ADDITIONAL NOTES
%     inside this file).
%
%   NOTE:
%     * Optional inputs use its DEFAULT value when not given or [].
%     * Optional outputs may or not be called.
%     * If X is a two element vector, min(X) will be the coordinate of the
%       first column and max(X) of the last column.
%     * If Y is a two element vector, min(Y) will be the coordinate of the
%       first row and max(Y) of the last row.
%     * If vector X-axis is decreasing U=fliplr(U) will be used.
%     * If vector Y-axis is decreasing U=flipud(U) will be used.
%     * When X or Y do not have a constant increasing/decreasing step, the
%       vertices of the color rectangules are set in the middle of each
%       pair of coordinates. For this reason its center may not lay on the
%       specified coordinate, except on the coordinates at the edges where
%       it always lays on the center.
%     * To get a non-scaled image (IMAGE instead of IMAGESC) use:
%         >> H = imagescnan(...);
%         >> set(H,'CDataMapping','direct')
%     * ADDITIONAL NOTES are included inside this file.
%
%   EXAMPLE:
%     % Compares with normal IMAGESC:
%      N     = 100;
%      PNaNs = 0.10;
%      U     = peaks(N);
%      U(round(1 + (N^2-1).*rand(N^2*PNaNs,1))) = NaN;         % Adds NaNs
%      subplot(221), imagesc(U)
%       title('With IMAGESC: ugly NaNs')
%      subplot(222), imagescnan(U) 
%       title('With IMAGESCNAN: uncolored NaNs')
%     % Compares with SPY:
%      subplot(223), spy(isnan(U))
%       title('SPY(isnan(U))')
%      subplot(224), imagescnan(isnan(U),'NaNMask',0), axis equal tight
%       title('SPY with IMAGESCNAN')
%     
%   SEE ALSO:
%     IMAGE, IMAGESC, COLORBAR, IMREAD, IMWRITE
%     and
%     CMAPPING, CBFREEZE by Carlos Vargas
%     at http://www.mathworks.com/matlabcentral/fileexchange
%
%
%   ---
%   MFILE:   imagescnan.m
%   VERSION: 2.1 (Aug 20, 2009) (<a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/authors/11258')">download</a>) 
%   MATLAB:  7.7.0.471 (R2008b)
%   AUTHOR:  Carlos Adrian Vargas Aguilera (MEXICO)
%   CONTACT: nubeobscura@hotmail.com

%   ADDITIONAL NOTES:
%     * I keep getting a kind of BUG with the edges of the patched NaNs. I
%       added two NOTE inside this program that may fix this problem.
%       Another way is to convert the intensity matrix U into RGB colors by
%       using the CMAPPING function, as used by the first version of this
%       program.
%     * Besides, if the matrix is too large, sometimes there is an
%       undocumented failure while drawing the patch NaNs. Is recommended
%       to use U = cmapping(U,[],'k','discrete') instead, and change the
%       CLIM to [min(U(:)) max(U(:))].
%     * The use of not homogeneous step interval X,Y axes is not
%       recommended because the program tries to put its value in the
%       middle of the colored rectangule (as IMAGESC does) and soetimes the
%       result may not be what the user wants. So this is for experimental
%       use only.

%   REVISIONS:
%   1.0      Released. (Jun 30, 2008)
%   1.1      Fixed bug when CAXIS used. Colorbar freezed colormap. Fixed
%            bug in color vector input (Found by Greg King) and now 
%            accets RGB image as input. (Jul 14, 2008)
%   2.0      Totally rewritten code. Do not converts to RGB anymore. Do not
%            freezes the colormap anymore. Do not output any colorbar. New
%            X and Y variable steps accepted input. Now uses patches. (Jun
%            08, 2009)
%   2.1      Fixed bug with RGB input. Added a NOTE about the use of
%            CMAPPING. (Aug 20, 2009)

%   DISCLAIMER:
%   imagescnan.m is provided "as is" without warranty of any kind, under
%   the revised BSD license.

%   Copyright (c) 2008,2009 Carlos Adrian Vargas Aguilera


% INPUTS CHECK-IN
% -------------------------------------------------------------------------

% Initializes:
X    = [];
Y    = [];
CNAN = [];
MNAN = [];
ha   = [];

% Checks number of inputs:
if     nargin<1
 error('CVARGAS:imagescnan:notEnoughInputs',...
  'At least 1 input is required.')
elseif nargout>2
 error('CVARGAS:imagescnan:tooManyOutputs',...
  'At most 2 outputs are allowed.')
end

% Gets X,Y,U:
if ((nargin==1) || (nargin==2))
 U = varargin{1};
 varargin(1) = [];
else
 if (isnumeric(varargin{1}) && isnumeric(varargin{2}) && ...
   isnumeric(varargin{3}))
  X = varargin{1};
  Y = varargin{2};
  U = varargin{3};
  varargin(1:3) = [];
 else
  U = varargin{1};
  varargin(1) = [];
 end
end

% Check U:
ndim = ndims(U);
if     (ndim==2)
 [M,N]   = size(U);
 O = 1;
elseif (ndim==3)
 [M,N,O] = size(U);
 if (O~=3)
  error('CVARGAS:imagescnan:incorrectRgbImage',...
   'RGB image must be of size M-by-N-by-3.')
 end
else
 error('CVARGAS:imagescnan:incorrectImageSize',...
  'Image must be 2-dimensional or a 3-dim RGB image.')
end

% Check X:
aequal = true;    % Equal intervals on x-axis?
dX     = [];
if isempty(X)
 X = [1 N];
else
 if (ndims(X)>2)
  error('CVARGAS:imagescnan:incorrectXDims',...
   'X must be a vector or a matrix as a result of MESHGRID.')
 end
 if any(~isfinite(X(:)))
  error('CVARGAS:imagescnan:incorrectXValue',...
   'X elements must be numeric and finite.')
 end
 [Mx,Nx] = size(X);
 if ((Mx*Nx)==2)
  if X(2)<X(1)
   X = X([2 1]);
   for k = 1:O % Fixed bug Aug 2009
    U(:,:,k) = fliplr(U(:,:,k));
   end
  end 
 else
  if     ((Mx==M) && (Nx==N))
   % Checks if generated with MESHGRID:
   dX    = abs(X(2:M,:)-repmat(X(1,:),M-1,1));
   if any(abs(dX(:))>(eps*max(abs(dX(:)))*1000))
    error('CVARGAS:imagescnan:incorrectXMatrix',...
     'X matrix must be as generated by MESHGRID.')
   end
   X = X(1,:);
  elseif (~any([Mx Nx]==1) && ~((Mx*Nx)==N))
   error('CVARGAS:imagescnan:incorrectXSize',...
     'X must be an scalar or a matrix.')
  end     
  % Forces ascending x-axis:
  [X,I] = sort(X(:).');
  for k = 1:O % Fixed bug Aug 2009
   U(:,:,k) = U(:,I,k);
  end
  clear I
  % Checks equal intervals:
  dX = diff(X);
  if any(abs(dX(1)-dX(2:end))>(eps*max(dX)*1000))
   if aequal
    aequal = false;
   end
  else
   X  = [X(1) X(end)];
   dX = [];
  end
 end
end

% Check Y:
dY = [];
if isempty(Y)
 Y = [1 M];
else
 if (ndims(Y)>2)
  error('CVARGAS:imagescnan:incorrectYDims',...
   'Y must be a vector or a matrix as a result of MESHGRID.')
 end
 if any(~isfinite(Y(:)))
  error('CVARGAS:imagescnan:incorrectYValue',...
   'Y elements must be numeric and finite.')
 end
 [My,Ny] = size(Y);
 if ((My*Ny)==2)
  if Y(2)<Y(1)
   Y = Y([2 1]);
   for k = 1:O % Fixed bug Aug 2009
    U(:,:,k) = flipud(U(:,:,k));
   end
  end
 else
  if     ((My==M) && (Ny==N))
   % Checks if generated with MESHGRID:
   dY = abs(Y(:,2:N)-repmat(Y(:,1),1,N-1));
   if any(abs(dY(:))>(eps*max(abs(dY(:)))*1000))
    error('CVARGAS:imagescnan:incorrectYMatrix',...
     'Y matrix must be as generated by MESHGRID.')
   end
   Y = Y(:,1);
  elseif (~any([My Ny]==1) && ~((My*Ny)==M))
   error('CVARGAS:imagescnan:incorrectYSize',...
     'Y must be an scalar or a matrix.')
  end     
  % Forces ascending y-axis:
  [Y,I] = sort(Y(:).');
  for k = 1:O % Fixed bug Aug 2009
   U(:,:,k) = U(I,:,k);
  end
  clear I
  % Checks equal intervals:
  dY = diff(Y);
  if any(abs(dY(1)-dY(2:end))>(eps*max(dY)*1000))
   if aequal
    aequal = false;
   end
  else
   Y  = [Y(1) Y(end)];
   dY = [];
  end
 end
end

% Checks varargin:
ind  = [];
Nopt = length(varargin); 
for k = 1:Nopt-1
 if (~isempty(varargin{k}) && ischar(varargin{k}))
  if     strncmpi(varargin{k},'NanColor',4)
   CNAN = varargin{k+1};
   ind  = [ind k k+1];
  elseif strncmpi(varargin{k},'NanMask',4)
   MNAN = varargin{k+1};
   ind  = [ind k k+1];
  elseif (strncmpi(varargin{k},'Parent',2) && isempty(CNAN))
   try
    CNAN = get(varargin{k+1},'Color');
    ha   = varargin{k+1};
   catch
    error('CVARGAS:imagescnan:incorrectParentHandle',...
     '''Parent'' must be a valid axes handle.')
   end
  end
 end
end
varargin(ind) = [];
Nargin = length(varargin);

% Check ha:
if isempty(ha)
 ha = gca;
end

% Check CNAN:
if     isempty(CNAN)
 CNAN = get(ha,'Color');
elseif ischar(CNAN)
 switch lower(CNAN)
  case 'y', CNAN = [1 1 0];
  case 'm', CNAN = [1 0 0];
  case 'c', CNAN = [0 1 1];
  case 'r', CNAN = [1 0 0];
  case 'g', CNAN = [0 1 0];
  case 'b', CNAN = [0 0 1];
  case 'w', CNAN = [1 1 1];
  case 'k', CNAN = [0 0 0];
  otherwise
   error('CVARGAS:imagescnan:incorrectNancString',...
    'Color string must be a valid color identifier. One of ''ymcrgbwk''.')
 end
elseif isnumeric(CNAN) && (length(CNAN)==3)
 CNAN = CNAN(:).'; % Forces row vector.
else
 error('CVARGAS:imagescnan:incorrectNancInput',...
  'Not recognized CNAN input.')
end

% Check MNAN:
if isempty(MNAN)
 MNAN = any(~isfinite(U),3);
else
 if (ndims(MNAN)==2)
  [Mm,Nm] = size(MNAN);
  if     ((Mm*Nm)==1)
   MNAN = (any(~isfinite(U),3) | any(U==MNAN,3));
  elseif ((Mm==M) && (Nm==N) && islogical(MNAN))
   MNAN = (any(~isfinite(U),3) | MNAN);
  else
   error('CVARGAS:imagescnan:incorrectNanmSize',...
   'MNAN must be an scalar or a logical matrix of size M-by-N.')
  end
 else
  error('CVARGAS:imagescnan:incorrectNanmDims',...
   'MNAN must be an scalar or a matrix.')
 end
end

% -------------------------------------------------------------------------
% MAIN
% -------------------------------------------------------------------------

% Generates the image:
if aequal
 % IMAGESC way.
 H = imagesc(X,Y,U,varargin{:});

else
 % PATCH way.

 % Check clim:
 if (rem(Nargin,2)==1)
  clim          = varargin{end};
  varargin(end) = [];
  if ((length(clim)~=2) || (clim(1)>clim(2)))
   error('CVARGAS:imagescnan:incorrectClimInput',...
    'clim must be a 2 element increasing vector.')
  end
 else
  clim = [];
 end

 % Generates vertices between coordinates (coordinates may not be at the
 % center of these vertices): 
 if (length(X)~=N)
  X  = (0:N-1)*((X(2)-X(1))/(N-1)) + X(1);
 end
 if (length(Y)~=M)
  Y  = (0:M-1)*((Y(2)-Y(1))/(M-1)) + Y(1);
 end
 if isempty(dX)
  dX = diff(X);
 end
 if isempty(dY)
  dY = diff(Y);
 end
 [X,Y] = meshgrid([X(1)-dX(1)/2 X+dX([1:N-1 N-1])/2],...
                  [Y(1)-dY(1)/2 Y+dY([1:M-1 M-1])/2]);
 
 % Generates faces:
 ind              = (1:(M+1)*N)';
 ind(M+1:M+1:end) = [];
 
 % Generates patches:
 H = patch(...
  'Vertices'       ,[X(:) Y(:)],...
  'Faces'          ,[ind ind+1 ind+M+2 ind+M+1],...
  'FaceVertexCData',U(:),...
  'FaceColor'      ,'flat',...
  'EdgeColor'      ,'none',... % NOTE: Sometimes this is not required.
  varargin{:});
 set(ha,...
  'YDir' ,'reverse',...
  'View' ,[0 90],...
  'Box'  ,'on',...
  'Layer','top')
 axis(ha,'tight')
 
 % Sets clim:
 if ~isempty(clim)
  set(ha,'CLim',clim)
 else
  set(ha,'CLimMode','auto')
 end
 
end

% Adds NaNs patches:
if any(MNAN(:))
 if aequal
  % dX and dY is constant:
  [MNAN,NNAN] = ind2sub([M,N],find(MNAN));
  Nnan        = length(MNAN);
  dX   = (X(2)-X(1))/(N-1)/2;
  dY   = (Y(2)-Y(1))/(M-1)/2;
  HNAN = patch(repmat((X(1)+(NNAN(:)'-1)*(2*dX)),4,1) + ...
                                       (repmat([-1 1 1 -1]'*dX,1,Nnan)),...
               repmat((Y(1)+(MNAN(:)'-1)*(2*dY)),4,1) + ...                                       
                                       (repmat([1 1 -1 -1]'*dY,1,Nnan)),...
               CNAN,...
               'EdgeColor',CNAN,... 'EdgeColor','none',... 
               varargin{1:Nargin-rem(Nargin,2)});
 else
  % dX and/or dY is not constant:
  MNAN = find(MNAN);
  HNAN = patch(...
  'Vertices'       ,[X(:) Y(:)],...
  'Faces'          ,[ind(MNAN) ind(MNAN)+1 ind(MNAN)+M+2 ind(MNAN)+M+1],...
  'FaceColor'      ,CNAN,...
  'EdgeColor'      ,'none',... 'EdgeColor',CNAN,... % NOTE: may be better?
  varargin{:});
 end
else
 HNAN = [];
end

 
% OUTPUTS CHECK-OUT
% -------------------------------------------------------------------------

% Clears outputs?:
if (nargout==0)
 clear H
end


% [EOF]   imagescnan.m