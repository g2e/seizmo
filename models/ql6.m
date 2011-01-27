function [mout]=ql6(varargin)
%QL6    Quality factor model of the Earth by Durek & Ekstrom 1996
%
%    Usage:    model=ql6()
%              model=ql6(...,'model',model,...)
%              model=ql6(...,'depths',depths,...)
%              model=ql6(...,'dcbelow',true|false,...)
%              model=ql6(...,'range',[top bottom],...)
%              model=ql6(...,'crust',true|false,...)
%              model=ql6(...,'ocean',true|false,...)
%
%    Description:
%     MODEL=QL6() returns the 1D radial Q factor model of Durek & Ekstrom
%     for the whole Earth.  The struct has the following fields:
%      MODEL.name      - model name ('QL6')
%           .ocean     - true if OCEAN & CRUST are both TRUE
%           .crust     - true if CRUST is TRUE
%           .isotropic - always true here
%           .refperiod - always 1sec here
%           .flattened - always false here (see FLATTEN_1DMODEL)
%           .depth     - km depths from 0 to 6371
%           .qk        - bulk moduli quality factor
%           .qu        - shear moduli quality factor
%     Note that the model includes repeated depths at discontinuities.
%
%     MODEL=QL6(...,'MODEL',MODEL,...) replaces the values of the Qu & Qk
%     fields in the input model MODEL.  The depths in MODEL.depth are used
%     to determine the new values.  The .ocean and .crust fields are also
%     taken into account.  MODEL must not be flattened!
%
%     MODEL=QL6(...,'DEPTHS',DEPTHS,...) returns the model parameters only
%     at the depths in DEPTHS.  DEPTHS is assumed to be in km.  DEPTHS at
%     discontinuities return values from the deeper (bottom) side of the
%     discontinuity for the first time and from the top side for the second
%     time.  Depths can not be repeated more than twice and must be
%     monotonically non-decreasing.
%
%     MODEL=QL6(...,'DCBELOW',TRUE|FALSE,...) returns values from the
%     shallow (top) side of the discontinuity the first time a depth is
%     given at one (using the DEPTHS option) if DCBELOW is FALSE.  The
%     default is TRUE (returns value from bottom-side the first time).  The
%     second time a depth is used, the opposite side is given.
%
%     MODEL=QL6(...,'RANGE',[TOP BOTTOM],...) specifies the range of
%     depths that known model parameters are returned.  [TOP BOTTOM] must
%     be a 2 element array in km.  Note this does not block depths given by
%     the DEPTHS option.  It does limit the depths replaced in a model
%     passed in with the MODEL option.
%
%     MODEL=QL6(...,'CRUST',TRUE|FALSE,...) indicates if the crust of
%     PREM is to be removed or not.  Setting CRUST to FALSE will return a
%     crustless model (the mantle is extended to the surface).  This will
%     also remove the ocean.
%
%     MODEL=QL6(...,'OCEAN',TRUE|FALSE,...) indicates if the ocean of
%     PREM is to be removed or not.  Setting CRUST to FALSE will return a
%     oceanless model (the crust (or mantle when CRUST is FALSE) is
%     extended to the surface).  The default is FALSE (no ocean).
%
%    Notes:
%     - QL6 reference:
%        Durek & Ekstrom 1996, A radial model of anelasticity
%        consistent with long-period surface-wave attenuation
%
%    Examples:
%     % Add QL6 to PREM, AK135, IASP91:
%     prem_ql6=ql6('model',prem);
%     ak135_ql6=ql6('model',ak135);
%     iasp91_ql6=ql6('model',iasp91);
%
%    See also: QLM9, PREM, PREM_PERFECT

%     Version History:
%        Sep. 19, 2010 - initial version
%        Jan. 25, 2011 - fix nan/inf bug when model given as input
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 25, 2011 at 10:00 GMT

% todo:

% check nargin
if(mod(nargin,2))
    error('seizmo:ql6:badNumInputs',...
        'Unpaired Option/Value!');
end

% option defaults
varargin=[{'m' [] 'd' [] 'b' true 'c' true 'o' false 'r' [0 6371]} ...
    varargin];

% check options
modelgiven=false;
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
        case {'m' 'mod' 'model'}
            if(~isempty(varargin{i+1}) ...
                    && ~isempty(chk1dmodel(varargin{i+1})))
                error(chk1dmodel(varargin{i+1}));
            end
            if(~isempty(varargin{i+1}))
                mout=varargin{i+1};
                modelgiven=true;
                ocean=mout.ocean;
                crust=mout.crust;
                depths=mout.depth;
                if(mout.flattened)
                    error('seizmo:ql6:badInput',...
                        'MODEL can not be flattened!');
                end
            else
                clear mout;
                mout.name='';
            end
        case {'o' 'ocean'}
            if(skip); continue; end
            if(modelgiven)
                warning('seizmo:ql6:oddInput',...
                    'Setting OCEAN with a given model (bad idea)!');
            end
            if(~islogical(varargin{i+1}) || ~isscalar(varargin{i+1}))
                error('seizmo:ql6:badOCEAN',...
                    'OCEAN must be a TRUE or FALSE!');
            end
            ocean=varargin{i+1};
        case {'d' 'dep' 'depth' 'depths'}
            if(~isempty(varargin{i+1}))
                if(~isreal(varargin{i+1}) || any(varargin{i+1}<0 ...
                        | varargin{i+1}>6371) || any(isnan(varargin{i+1})))
                    error('seizmo:ql6:badDEPTHS',...
                        ['DEPTHS must be real-valued km depths within ' ...
                        'the range [0 6371] in km!']);
                elseif(any(diff(varargin{i+1})<0))
                    error('seizmo:ql6:badDEPTHS',...
                        'DEPTHS must be monotonically non-increasing!');
                elseif(any(histc(varargin{i+1},...
                        varargin{i+1}([find(diff(varargin{i+1}));end]))>3))
                    error('seizmo:ql6:badDEPTHS',...
                        'DEPTHS has values repeated 3+ times!');
                end
            end
            depths=varargin{i+1}(:);
        case {'dcb' 'dc' 'below' 'b' 'dcbelow'}
            if(skip); continue; end
            if(~islogical(varargin{i+1}) || ~isscalar(varargin{i+1}))
                error('seizmo:ql6:badDCBELOW',...
                    'DCBELOW must be a TRUE or FALSE!');
            end
            dcbelow=varargin{i+1};
        case {'c' 'cru' 'crust'}
            if(skip); continue; end
            if(modelgiven)
                warning('seizmo:ql6:oddInput',...
                    'Setting CRUST with a given model (bad idea)!');
            end
            if(~islogical(varargin{i+1}) || ~isscalar(varargin{i+1}))
                error('seizmo:ql6:badCRUST',...
                    'CRUST must be a TRUE or FALSE!');
            end
            crust=varargin{i+1};
        case {'r' 'rng' 'range'}
            if(skip); continue; end
            if(~isreal(varargin{i+1}) || numel(varargin{i+1})~=2)
                error('seizmo:ql6:badRANGE',...
                    ['RANGE must be a 2 element vector specifying ' ...
                    '[TOP BOTTOM] in km!']);
            end
            range=sort(varargin{i+1});
        otherwise
            error('seizmo:prem:badOption',...
                'Unknown Option: %s',varargin{i});
    end
end

% the ql6 model
% - includes the ocean
model=[
   0 inf inf
   3 inf inf
   3 300 inf
24.4 300 inf
24.4 191 943
  80 191 943
  80  70 943
 220  70 943
 220 165 943
 670 165 943
 670 355 inf
2891 355 inf
2891 inf inf
5150 inf inf
5150 104 inf
6371 104 inf];

% remove ocean if desired
if(~ocean)
    model(1,:)=[0 300 inf];
    model(2:3,:)=[];
end

% remove crust & ocean if desired
if(~crust)
    model(1,:)=[0 191 943];
    model(2:(4+2*ocean),:)=[];
end

% interpolate depths if desired
if(~isempty(depths))
    if(~isempty(mout.name))
        % isolate depths that matter
        idx=depths>=range(1) & depths<=range(2);
        depths=depths(idx);
        
        % get new qu values
        [bot,top]=interpdc1(model(:,1),model(:,2:end),depths);
        
        % replace the first of repeated values with top
        [tidx,tidx]=unique(depths,'first');
        bot(tidx,:)=top(tidx,:);
        
        % insert info
        mout.name=[mout.name '+QL6'];
        mout.qu(idx)=bot(:,1);
        mout.qk(idx)=bot(:,2);
        
        % Inf values become NaNs in interpolation, so change back
        mout.qu(isnan(mout.qu))=inf;
        mout.qk(isnan(mout.qk))=inf;
        
        % update ocean & crust (possibly changed - we warned earlier)
        mout.ocean=ocean & crust;
        mout.crust=crust;
    else
        [bot,top]=interpdc1(model(:,1),model(:,2:end),depths);
        if(dcbelow)
            [tidx,tidx]=unique(depths);
            top(tidx,:)=bot(tidx,:);
            model=[depths top];
        else
            [tidx,tidx]=unique(depths,'first');
            bot(tidx,:)=top(tidx,:);
            model=[depths bot];
        end
    end
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
        vtop=interp1q(...
            model(idx1-1:idx1,1),model(idx1-1:idx1,2:end),range(1));
    end
    if(tf(2))
        idx2=idx2+1;
    else
        vbot=interp1q(...
            model(idx2:idx2+1,1),model(idx2:idx2+1,2:end),range(2));
    end
    
    % clip model
    model=model(idx1:idx2,:);
    
    % pad range knots if not there
    if(~tf(1)); model=[range(1) vtop; model]; end
    if(~tf(2)); model=[model; range(2) vbot]; end
end

% array to struct
if(isempty(mout.name))
    % fix interpolation of infinity
    model(isnan(model))=inf;
    
    % make output struct
    mout.name='QL6';
    mout.ocean=ocean & crust;
    mout.crust=crust;
    mout.isotropic=true;
    mout.refperiod=1;
    mout.flattened=false;
    mout.depth=model(:,1);
    mout.qu=model(:,2);
    mout.qk=model(:,3);
end

end
