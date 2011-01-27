function [mout]=qlm9(varargin)
%QLM9    Quality factor model of the Earth by Lawrence & Wysession 2006
%
%    Usage:    model=qlm9()
%              model=qlm9(...,'model',model,...)
%              model=qlm9(...,'depths',depths,...)
%              model=qlm9(...,'dcbelow',true|false,...)
%              model=qlm9(...,'range',[top bottom],...)
%
%    Description:
%     MODEL=QLM9() returns the 1D radial Q factor model of Lawrence &
%     Wysession for the whole mantle.  The struct has the following fields:
%      MODEL.name      - model name ('QL6')
%           .ocean     - always false here
%           .crust     - always false here
%           .isotropic - always true here
%           .refperiod - always 1sec here
%           .flattened - always false here (see FLATTEN_1DMODEL)
%           .depth     - km depths from 0 to 2891
%           .qu        - shear moduli quality factor
%     Note that the model includes repeated depths at discontinuities.
%
%     MODEL=QLM9(...,'MODEL',MODEL,...) replaces the values of the Qu
%     field in the input model MODEL.  The depths in MODEL.depth are used
%     to determine the new values.  Note that the .ocean and .crust fields
%     are not accounted for (any Qu values in those layers are replaced
%     with the QLM9 lithospheric value).  Depths below the mantle (2891km)
%     are not replaced.
%
%     MODEL=QLM9(...,'DEPTHS',DEPTHS,...) returns the model parameters only
%     at the depths in DEPTHS.  DEPTHS is assumed to be in km.  DEPTHS at
%     discontinuities return values from the deeper (bottom) side of the
%     discontinuity for the first time and from the top side for the second
%     time.  Depths can not be repeated more than twice and must be
%     monotonically non-decreasing.
%
%     MODEL=QLM9(...,'DCBELOW',TRUE|FALSE,...) returns values from the
%     shallow (top) side of the discontinuity the first time a depth is
%     given at one (using the DEPTHS option) if DCBELOW is FALSE.  The
%     default is TRUE (returns value from bottom-side the first time).  The
%     second time a depth is used, the opposite side is given.
%
%     MODEL=QLM9(...,'RANGE',[TOP BOTTOM],...) specifies the range of
%     depths that known model parameters are returned.  [TOP BOTTOM] must
%     be a 2 element array in km.  Note this does not block depths given by
%     the DEPTHS option.  It does limit the depths replaced in a model
%     passed in with the MODEL option.
%
%    Notes:
%     - QLM9 reference:
%        Lawrence & Wysession 2006, QLM9: A new radial quality
%        factor (Qu) model for the lower mantle
%
%    Examples:
%     % Add QLM9 to PREM, AK135, IASP91:
%     prem_qlm9=qlm9('model',prem);
%     ak135_qlm9=qlm9('model',ak135);
%     iasp91_qlm9=qlm9('model',iasp91);
%
%    See also: QL6, PREM, PREM_PERFECT

%     Version History:
%        Sep. 19, 2010 - initial version
%        Jan. 25, 2011 - minor code cleanup for inf values
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 25, 2011 at 10:00 GMT

% todo:

% check nargin
if(mod(nargin,2))
    error('seizmo:qlm9:badNumInputs',...
        'Unpaired Option/Value!');
end

% option defaults
varargin=[{'m' [] 'd' [] 'b' true 'r' [0 2891]} varargin];

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
        case {'m' 'mod' 'model'}
            if(~isempty(varargin{i+1}) ...
                    && ~isempty(chk1dmodel(varargin{i+1})))
                error(chk1dmodel(varargin{i+1}));
            end
            if(~isempty(varargin{i+1}))
                mout=varargin{i+1};
                depths=mout.depth;
                if(mout.flattened)
                    error('seizmo:qlm9:badInput',...
                        'MODEL can not be flattened!');
                end
            else
                clear mout;
                mout.name='';
            end
        case {'d' 'dep' 'depth' 'depths'}
            if(~isempty(varargin{i+1}))
                if(~isreal(varargin{i+1}) || any(varargin{i+1}<0 ...
                        | varargin{i+1}>6371) || any(isnan(varargin{i+1})))
                    error('seizmo:qlm9:badDEPTHS',...
                        ['DEPTHS must be real-valued km depths within ' ...
                        'the range [0 6371] in km!']);
                elseif(any(diff(varargin{i+1})<0))
                    error('seizmo:qlm9:badDEPTHS',...
                        'DEPTHS must be monotonically non-increasing!');
                elseif(any(histc(varargin{i+1},...
                        varargin{i+1}([find(diff(varargin{i+1}));end]))>3))
                    error('seizmo:qlm9:badDEPTHS',...
                        'DEPTHS has values repeated 3+ times!');
                end
            end
            depths=varargin{i+1}(:);
        case {'dcb' 'dc' 'below' 'b' 'dcbelow'}
            if(skip); continue; end
            if(~islogical(varargin{i+1}) || ~isscalar(varargin{i+1}))
                error('seizmo:qlm9:badDCBELOW',...
                    'DCBELOW must be a TRUE or FALSE!');
            end
            dcbelow=varargin{i+1};
        case {'r' 'rng' 'range'}
            if(skip); continue; end
            if(~isreal(varargin{i+1}) || numel(varargin{i+1})~=2)
                error('seizmo:qlm9:badRANGE',...
                    ['RANGE must be a 2 element vector specifying ' ...
                    '[TOP BOTTOM] in km!']);
            end
            range=sort(varargin{i+1});
        otherwise
            error('seizmo:prem:badOption',...
                'Unknown Option: %s',varargin{i});
    end
end

% the qlm9 model
model=[
   0 600
  80 600
  80  80
 220  80
 220 143
 400 143
 400 276
 670 276
 670 362
1000 362
1000 325
1350 325
1350 287
1700 287
1700 307
2050 307
2050 383
2400 383
2400 459
2700 459
2700 452
2800 452
2800 278
2891 278];

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
        bot(tidx)=top(tidx);
        
        % insert info
        mout.name=[mout.name '+QLM9'];
        mout.qu(idx)=bot;
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

% fix interpolation of infinity
model(isnan(model))=inf;

% array to struct
if(isempty(mout.name))
    % fix interpolation of infinity
    model(isnan(model))=inf;
    
    % make output struct
    mout.name='QLM9';
    mout.ocean=false;
    mout.crust=false;
    mout.isotropic=true;
    mout.refperiod=1;
    mout.flattened=false;
    mout.depth=model(:,1);
    mout.qu=model(:,2);
end

end
