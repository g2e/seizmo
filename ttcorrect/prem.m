function [mout]=prem(varargin)
%PREM    Returns the PREM Earth Model
%
%    Usage:    model=prem
%              model=prem(...,'depths',depths,...)
%              model=prem(...,'range',[top bottom],...)
%              model=prem(...,'crust',true|false,...)
%              model=prem(...,'ocean',true|false,...)
%
%    Description: MODEL=PREM returns a struct containing the isotropic
%     version of the 1D radial Earth model PREM.  Note that by default the
%     ocean has been replaced with an extended upper crust.  The struct has
%     the following fields:
%      MODEL.name  - model name ('PREM')
%           .depth - km depths from 0 to 6371
%           .vp    - isotropic p-wave velocity at 1s periods (km/s)
%           .vs    - isotropic s-wave velocity at 1s periods (km/s)
%           .rho   - density (g/cm^3)
%           .qp    - p-wave quality factor
%           .qs    - s-wave quality factor
%     Note that the model includes repeated depths at discontinuities.
%
%     MODEL=PREM(...,'DEPTHS',DEPTHS,...) returns the model parameters
%     only at the depths in DEPTHS.  DEPTHS is assumed to be in km.  The
%     model parameters are found by linear interpolation between known
%     values.  DEPTHS at discontinuities return values from the under side
%     of the discontinuity (this currently cannot be adjusted).
%
%     MODEL=PREM(...,'RANGE',[TOP BOTTOM],...) specifies the range of
%     depths that known model parameters are returned.  [TOP BOTTOM] must
%     be a 2 element array in km.  Note this does not block depths given by
%     the DEPTHS option.
%
%     MODEL=PREM(...,'CRUST',TRUE|FALSE,...) indicates if the crust of
%     PREM is to be removed or not.  Setting CRUST to FALSE will return a
%     crustless model (the mantle is extended to the surface using linear
%     interpolation).
%
%     MODEL=PREM(...,'OCEAN',TRUE|FALSE,...) indicates if the ocean of
%     PREM is to be removed or not.  Setting CRUST to FALSE will return a
%     oceanless model (the crust (or mantle when CRUST is FALSE) is
%     extended to the surface using linear interpolation).  The default is
%     FALSE (no ocean).
%
%    Notes:
%     - Note that this is not a perfect match to PREM
%
%    Examples:
%     Plot parameters for the CMB region:
%      model=prem('r',[2600 3400]);
%      figure;
%      plot(model.depth,model.vp,'r',...
%           model.depth,model.vs,'g',...
%           model.depth,model.rho,'b','linewidth',2);
%      title('PREM')
%      legend({'Vp' 'Vs' '\rho'})
%
%    See also: AK135, IASP91

%     Version History:
%        May  19, 2010 - initial version
%        May  20, 2010 - discon on edge handling, quicker
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  20, 2010 at 16:25 GMT

% todo:

% check nargin
if(mod(nargin,2))
    error('seizmo:prem:badNumInputs',...
        'Unpaired Option/Value!');
end

% option defaults
varargin=[{'d' [] 'c' true 'o' false 'r' [0 6371]} varargin];

% check options
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:prem:badOption',...
        'All Options must be specified with a string!');
end
for i=1:2:numel(varargin)
    % skip empty
    skip=false;
    if(isempty(varargin{i+1})); skip=true; end

    % check option is available
    switch lower(varargin{i})
        case {'o' 'ocean'}
            if(skip); continue; end
            if(~islogical(varargin{i+1}) || ~isscalar(varargin{i+1}))
                error('seizmo:prem:badOCEAN',...
                    'OCEAN must be a TRUE or FALSE!');
            end
            ocean=varargin{i+1};
        case {'d' 'dep' 'depth' 'depths'}
            if(~isempty(varargin{i+1}) && (~isreal(varargin{i+1}) ...
                    || any(varargin{i+1}<0 | varargin{i+1}>6371)))
                error('seizmo:prem:badDEPTHS',...
                    ['DEPTHS must be real-valued km depths within ' ...
                    'the range [0 6371] in km!']);
            end
            depths=varargin{i+1}(:);
        case {'c' 'cru' 'crust'}
            if(skip); continue; end
            if(~islogical(varargin{i+1}) || ~isscalar(varargin{i+1}))
                error('seizmo:prem:badCRUST',...
                    'CRUST must be a TRUE or FALSE!');
            end
            crust=varargin{i+1};
        case {'r' 'rng' 'range'}
            if(skip); continue; end
            if(~isreal(varargin{i+1}) || numel(varargin{i+1})~=2)
                error('seizmo:prem:badRANGE',...
                    ['RANGE must be a 2 element vector specifying ' ...
                    '[TOP BOTTOM] in km!']);
            end
            range=sort(varargin{i+1});
        otherwise
            error('seizmo:prem:badOption',...
                'Unknown Option: %s',varargin{i});
    end
end

% the prem model
model=[
    0.00     1.45000   0.00000   1.02000   57822.0       0.0
    3.00     1.45000   0.00000   1.02000   57822.0       0.0
    3.00     5.80000   3.20000   2.60000    1456.0     600.0
   15.00     5.80000   3.20000   2.60000    1456.0     600.0
   15.00     6.80000   3.90000   2.90000    1350.0     600.0
   24.40     6.80000   3.90000   2.90000    1350.0     600.0
   24.40     8.11061   4.49094   3.38076    1446.0     600.0
   40.00     8.10119   4.48486   3.37906    1446.0     600.0
   60.00     8.08907   4.47715   3.37688    1447.0     600.0
   80.00     8.07688   4.46953   3.37471     195.0      80.0
  115.00     8.05540   4.45643   3.37091     195.0      80.0
  150.00     8.03370   4.44361   3.36710     195.0      80.0
  185.00     8.01180   4.43108   3.36330     195.0      80.0
  220.00     7.98970   4.41885   3.35950     195.0      80.0
  220.00     8.55896   4.64391   3.43578     362.0     143.0
  265.00     8.64552   4.67540   3.46264     365.0     143.0
  310.00     8.73209   4.70690   3.48951     367.0     143.0
  355.00     8.81867   4.73840   3.51639     370.0     143.0
  400.00     8.90522   4.76989   3.54325     372.0     143.0
  400.00     9.13397   4.93259   3.72378     366.0     143.0
  450.00     9.38990   5.07842   3.78678     365.0     143.0
  500.00     9.64588   5.22428   3.84980     364.0     143.0
  550.00     9.90185   5.37014   3.91282     363.0     143.0
  600.00    10.15782   5.51602   3.97584     362.0     143.0
  635.00    10.21203   5.54311   3.98399     362.0     143.0
  670.00    10.26622   5.57020   3.99214     362.0     143.0
  670.00    10.75131   5.94508   4.38071     759.0     312.0
  721.00    10.91005   6.09418   4.41241     744.0     312.0
  771.00    11.06557   6.24046   4.44317     730.0     312.0
  871.00    11.24490   6.31091   4.50372     737.0     312.0
  971.00    11.41560   6.37813   4.56307     743.0     312.0
 1071.00    11.57828   6.44232   4.62129     750.0     312.0
 1171.00    11.73357   6.50370   4.67844     755.0     312.0
 1271.00    11.88209   6.56250   4.73460     761.0     312.0
 1371.00    12.02445   6.61891   4.78983     766.0     312.0
 1471.00    12.16126   6.67317   4.84422     770.0     312.0
 1571.00    12.29316   6.72548   4.89783     775.0     312.0
 1671.00    12.42075   6.77606   4.95073     779.0     312.0
 1771.00    12.54466   6.82512   5.00299     784.0     312.0
 1871.00    12.66550   6.87289   5.05469     788.0     312.0
 1971.00    12.78389   6.91957   5.10590     792.0     312.0
 2071.00    12.90045   6.96538   5.15669     795.0     312.0
 2171.00    13.01579   7.01053   5.20713     799.0     312.0
 2271.00    13.13055   7.05525   5.25729     803.0     312.0
 2371.00    13.24532   7.09974   5.30724     807.0     312.0
 2471.00    13.36074   7.14423   5.35706     811.0     312.0
 2571.00    13.47742   7.18892   5.40681     815.0     312.0
 2671.00    13.59597   7.23403   5.45657     819.0     312.0
 2741.00    13.68041   7.26597   5.49145     822.0     312.0
 2771.00    13.68753   7.26575   5.50642     823.0     312.0
 2871.00    13.71168   7.26486   5.55641     826.0     312.0
 2891.00    13.71660   7.26466   5.56645     826.0     312.0
 2891.00     8.06482   0.00000   9.90349   57822.0       0.0
 2971.00     8.19939   0.00000  10.02940   57822.0       0.0
 3071.00     8.36019   0.00000  10.18134   57822.0       0.0
 3171.00     8.51298   0.00000  10.32726   57822.0       0.0
 3271.00     8.65805   0.00000  10.46727   57822.0       0.0
 3371.00     8.79573   0.00000  10.60152   57822.0       0.0
 3471.00     8.92632   0.00000  10.73012   57822.0       0.0
 3571.00     9.05015   0.00000  10.85321   57822.0       0.0
 3671.00     9.16752   0.00000  10.97091   57822.0       0.0
 3771.00     9.27867   0.00000  11.08335   57822.0       0.0
 3871.00     9.38418   0.00000  11.19067   57822.0       0.0
 3971.00     9.48409   0.00000  11.29298   57822.0       0.0
 4017.00     9.57881   0.00000  11.39042   57822.0       0.0
 4171.00     9.66865   0.00000  11.48311   57822.0       0.0
 4271.00     9.75393   0.00000  11.57119   57822.0       0.0
 4371.00     9.83496   0.00000  11.65478   57822.0       0.0
 4471.00     9.91206   0.00000  11.73401   57822.0       0.0
 4571.00     9.98554   0.00000  11.80900   57822.0       0.0
 4671.00    10.05572   0.00000  11.87990   57822.0       0.0
 4771.00    10.12291   0.00000  11.94682   57822.0       0.0
 4871.00    10.18743   0.00000  12.00989   57822.0       0.0
 4971.00    10.24959   0.00000  12.06924   57822.0       0.0
 5071.00    10.30971   0.00000  12.12500   57822.0       0.0
 5149.50    10.35568   0.00000  12.16634   57822.0       0.0
 5149.50    11.02827   3.50432  12.76360     445.0      85.0
 5171.00    11.03643   3.51002  12.77493     445.0      85.0
 5271.00    11.07249   3.53522  12.82501     443.0      85.0
 5371.00    11.10542   3.55823  12.87073     440.0      85.0
 5471.00    11.13521   3.57905  12.91211     439.0      85.0
 5571.00    11.16186   3.59767  12.94912     437.0      85.0
 5671.00    11.18538   3.61411  12.98178     436.0      85.0
 5771.00    11.20576   3.62835  13.01009     434.0      85.0
 5871.00    11.22301   3.64041  13.03404     433.0      85.0
 5971.00    11.23712   3.65027  13.05364     432.0      85.0
 6071.00    11.24809   3.65794  13.06888     432.0      85.0
 6171.00    11.25593   3.66342  13.07977     431.0      85.0
 6271.00    11.26064   3.66670  13.08630     431.0      85.0
 6371.00    11.26220   3.66780  13.08848     431.0      85.0];

% remove ocean if desired
if(~ocean)
    % linearly extrapolated to the surface
    model(1,:)=[0.00     5.80000   3.20000   2.60000    1456.0     600.0];
    model(2:3,:)=[];
end

% remove crust if desired
if(~crust)
    % linearly extrapolated to the surface
    model(1,:)=[0 8.1253 4.5004 3.3834 1446.0 600.0];
    model(2:(4+2*ocean),:)=[];
end

% interpolate depths if desired
if(~isempty(depths))
    %depths=depths(depths>=range(1) & depths<=range(2));
    model=interpdc1(model(:,1),model(:,2:end),depths);
    model=[depths model];
else
    % get index range (assumes depths are always non-decreasing in model)
    idx1=find(model(:,1)>range(1),1);
    idx2=find(model(:,1)<range(2),1,'last');
    
    % are range points amongst the knots?
    tf=ismember(range,model(:,1));
    
    % if they are, just use the knot, otherwise interpolate
    if(tf(1))
        idx1=idx1-1;
    else
        vtop=interp1q(model(idx1-1:idx1,1),model(idx1-1:idx1,2:end),range(1));
    end
    if(tf(2))
        idx2=idx2+1;
    else
        vbot=interp1q(model(idx2:idx2+1,1),model(idx2:idx2+1,2:end),range(2));
    end
    
    % clip model
    model=model(idx1:idx2,:);
    
    % pad range knots if not there
    if(~tf(1)); model=[range(1) vtop; model]; end
    if(~tf(2)); model=[model; range(2) vbot]; end
end

% array to struct
mout.depth=model(:,1);
mout.vp=model(:,2);
mout.vs=model(:,3);
mout.rho=model(:,4);
mout.qp=model(:,5);
mout.qs=model(:,6);

end
