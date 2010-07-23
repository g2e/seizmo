function y = inpaintn(x,n,y0)

% INPAINTN Inpaint over missing data in N-D array
%   Y = INPAINTN(X) replaces the missing data in X by extra/interpolating
%   the non-missing elements. The non finite values (NaN or Inf) in X are
%   considered as missing data. X can be any N-D array.
%
%   INPAINTN (no input/output argument) runs the following 3-D example.
%
%   Important note:
%   --------------
%   INPAINTN uses an iterative process that converges toward the solution.
%   Y = INPAINTN(X,N) uses N iterations. By default, N = 100. If you
%   estimate that INPAINTN did not totally converge, increase N:
%   Y = INPAINTN(X,1000);
%
%   Y = INPAINTN(X,N,Y0) uses Y0 as initial guess. This could be useful if
%   you want to run the process a second time or if you have a GOOD guess
%   of the final result. By default, INPAINTN makes a nearest neighbor
%   interpolation (by using BWDIST) to obtain a rough guess.
%
%   Notes:
%   -----
%   <a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/4551')">INPAINT_NANS</a> and <a href="matlab:web('http://www.mathworks.com/matlabcentral/fileexchange/21214')">INPAINT_NANS3</a> are much faster than INPAINTN when
%   the number of NaN elements is (relatively) small. However, because
%   INPAINT_NANS and INPAINT_NANS3 both involve huge matrices, they can be
%   limited with large datasets.
%
%   Examples
%   --------
%
%     %% ---- Image ---- %%
%     onion = imread('onion.png');
%     I = randperm(numel(onion));
%     onionNaN = double(onion); onionNaN(I(1:round(numel(I)*0.5))) = NaN;
%     subplot(211), imshow(uint8(onionNaN)), title('Corrupted image - 50%')
%     for k=1:3, onion(:,:,k) = inpaintn(onionNaN(:,:,k)); end
%     subplot(212), imshow(uint8(onion)), title('Inpainted image')
%
%     %% ---- 3-D data ---- %%
%     load wind
%     xmin = min(x(:)); xmax = max(x(:));
%     zmin = min(z(:)); ymax = max(y(:));
%     %-- wind velocity
%     vel0 = interp3(sqrt(u.^2+v.^2+w.^2),1,'cubic');
%     x = interp3(x,1); y = interp3(y,1); z = interp3(z,1);
%     %-- remove randomly 90% of the data
%     I = randperm(numel(vel0));
%     velNaN = vel0;
%     velNaN(I(1:round(numel(I)*.9))) = NaN;
%     %-- inpaint using INPAINTN
%     vel = inpaintn(velNaN);
%     %-- display the results
%     subplot(221), imagesc(velNaN(:,:,15)), axis equal off
%     title('Corrupt plane, z = 15')
%     subplot(222), imagesc(vel(:,:,15)), axis equal off
%     title('Reconstructed plane, z = 15')    
%     subplot(223)
%     hsurfaces = slice(x,y,z,vel0,[xmin,100,xmax],ymax,zmin);
%     set(hsurfaces,'FaceColor','interp','EdgeColor','none')
%     hcont = contourslice(x,y,z,vel0,[xmin,100,xmax],ymax,zmin);
%     set(hcont,'EdgeColor',[.7,.7,.7],'LineWidth',.5)
%     view(3), daspect([2,2,1]), axis tight
%     title('Original data compared with...')
%     subplot(224)
%     hsurfaces = slice(x,y,z,vel,[xmin,100,xmax],ymax,zmin);
%     set(hsurfaces,'FaceColor','interp','EdgeColor','none')
%     hcont = contourslice(x,y,z,vel,[xmin,100,xmax],ymax,zmin);
%     set(hcont,'EdgeColor',[.7,.7,.7],'LineWidth',.5)
%     view(3), daspect([2,2,1]), axis tight
%     title('... reconstructed data')
%
%     %% --- 4-D data --- %%
%     [x1,x2,x3,x4] = ndgrid(-2:0.2:2);
%     z0 = x2.*exp(-x1.^2-x2.^2-x3.^2-x4.^2);
%     I = randperm(numel(z0));
%     % remove 50% of the data
%     zNaN = z0; zNaN(I(1:round(numel(I)*.5))) = NaN;
%     % reconstruct the data using INPAINTN
%     z = inpaintn(zNaN);
%     % display the results (for x4 = 0)
%     subplot(211)
%     zNaN(isnan(zNaN)) = 0.5;
%     slice(x2(:,:,:,1),x1(:,:,:,1),x3(:,:,:,1),zNaN(:,:,:,11),...
%        [-1.2 0.8 2],2,[-2 0.2])
%     title('Corrupt data, x4 = 0')
%     subplot(212)
%     slice(x2(:,:,:,1),x1(:,:,:,1),x3(:,:,:,1),z(:,:,:,11),...
%        [-1.2 0.8 2],2,[-2 0.2])
%     title('Reconstructed data')
%
%   See also GRIDDATAN, INPAINT_NANS, INPAINT_NANS3
%
%   -- Damien Garcia -- 2010/06
%   website: <a
%   href="matlab:web('http://www.biomecardio.com')">www.BiomeCardio.com</a>

test4DCTNandIDCTN

if nargin==0&&nargout==0, RunTheExample, return, end

x = double(x);
if nargin==1, n = 100; end

sizx = size(x);
d = ndims(x);
Lambda = zeros(sizx);
for i = 1:d
    siz0 = ones(1,d);
    siz0(i) = sizx(i);
    Lambda = bsxfun(@plus,Lambda,...
        cos(pi*(reshape(1:sizx(i),siz0)-1)/sizx(i)));
end
Lambda = -2*(d-Lambda);

% Initial condition
W = isfinite(x);
if nargin==3
    y = y0;
    s0 = 0;
else
    if any(~W(:))
        [y,s0] = InitialGuess(x,isfinite(x));
    else
        y = x;
        return
    end
end
x(~W) = 0;

% Smoothness parameters: from high to negligible values
s = logspace(s0,-3,n);

RF = 2; % relaxation factor
Lambda = Lambda.^2;

h = waitbar(0,'Inpainting...');
for i = 1:n
        Gamma = 1./(1+s(i)*Lambda);
        y = RF*idctn(Gamma.*dctn(W.*(x-y)+y)) + (1-RF)*y;
        waitbar(i/n,h)
end
close(h)

y(W) = x(W);

end

%% Test for DCTN and IDCTN
function test4DCTNandIDCTN
    if ~exist('dctn','file')
        error('MATLAB:smoothn:MissingFunction',...
            ['DCTN and IDCTN are required. Download DCTN <a href="matlab:web(''',...
            'http://www.biomecardio.com/matlab/dctn.html'')">here</a>.'])
    elseif ~exist('idctn','file')
        error('MATLAB:smoothn:MissingFunction',...
            ['DCTN and IDCTN are required. Download IDCTN <a href="matlab:web(''',...
            'http://www.biomecardio.com/matlab/idctn.html'')">here</a>.'])
    end
end

%% Initial Guess
function [z,s0] = InitialGuess(y,I)

if license('test','image_toolbox')
    %-- nearest neighbor interpolation
    [z,L] = bwdist(I);
    z = y;
    z(~I) = y(L(~I));
    s0 = 3;
else
    warning('MATLAB:inpaintn:InitialGuess',...
        ['BWDIST (Image Processing Toolbox) does not exist. ',...
        'The initial guess may not be optimal; additional',...
        ' iterations can thus be required to ensure complete',...
        ' convergence. Increase N value if necessary.'])
    z = y;
    z(~I) = mean(y(I));
    s0 = 6;
end

end

%% Example (3-D)
function RunTheExample
      load wind
      xmin = min(x(:)); xmax = max(x(:));
      zmin = min(z(:)); ymax = max(y(:));
      %-- wind velocity
      vel0 = interp3(sqrt(u.^2+v.^2+w.^2),1,'cubic');
      x = interp3(x,1); y = interp3(y,1); z = interp3(z,1);
      %-- remove randomly 90% of the data
      I = randperm(numel(vel0));
      velNaN = vel0;
      velNaN(I(1:round(numel(I)*.9))) = NaN;
      %-- inpaint using INPAINTN
      vel = inpaintn(velNaN);
      %-- display the results
      subplot(221), imagesc(velNaN(:,:,15)), axis equal off
      title('Corrupt plane, z = 15')
      subplot(222), imagesc(vel(:,:,15)), axis equal off
      title('Reconstructed plane, z = 15')    
      subplot(223)
      hsurfaces = slice(x,y,z,vel0,[xmin,100,xmax],ymax,zmin);
      set(hsurfaces,'FaceColor','interp','EdgeColor','none')
      hcont = contourslice(x,y,z,vel0,[xmin,100,xmax],ymax,zmin);
      set(hcont,'EdgeColor',[.7,.7,.7],'LineWidth',.5)
      view(3), daspect([2,2,1]), axis tight
      title('Original data compared with...')
      subplot(224)
      hsurfaces = slice(x,y,z,vel,[xmin,100,xmax],ymax,zmin);
      set(hsurfaces,'FaceColor','interp','EdgeColor','none')
      hcont = contourslice(x,y,z,vel,[xmin,100,xmax],ymax,zmin);
      set(hcont,'EdgeColor',[.7,.7,.7],'LineWidth',.5)
      view(3), daspect([2,2,1]), axis tight
      title('... reconstructed data')
end
