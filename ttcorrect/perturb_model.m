function [model]=perturb_model(model,varargin)
%PERTURB_MODEL    Perturbs 1D Earth models
%
%    Usage:    model=perturb_model(model,prop,value,...)
%
%    Description: MODEL=PERTURB_MODEL(MODEL,PROP,VALUE,...) will perturb
%     the property PROP of the 1D Earth model MODEL using the matrix VALUE
%     to create new layers in the model.  VALUE is a Nx5 matrix with the
%     following layout:
%       [depth1 value1 vtype1 gtype12 npts12;
%        depth2 value2 vtype2 gtype23 npts23;
%        ...
%        depthN valueN vtypeN 0       0     ]
%     Each row in VALUE defines a depth & value as well as the gradient and
%     npts between that depth and the next.  VALUE must have at least 2
%     rows.  The last row's gtype and npts entries must be 0 (for integrity
%     checking).  The control parameters are defined as follows:
%      depth -- depth in kilometers of new layer top/bottom
%      value -- value at layer top/bottom (see vtype for value type)
%      vtype -- value type:  0 -- true value
%                            1 -- value is relative to row above
%                            2 -- percent relative to row above
%                            3 -- value is relative to model
%                            4 -- percent relative to model
%               Note: 1 & 2 for 1st row are relative to model (== 3 & 4)
%      gtype -- gradient type:  0 -- linear
%                               + -- error function (strong grad at bottom)
%                               - -- error function (strong grad at top)
%               Note: higher values give stronger gradients
%      npts  -- number of points from layer top to bottom
%               Note: must be 2+ if not discontinuity or last row
%     Depths in the input model that are within the newly defined layers
%     are lost except for discontinuities (so that discontinuities for
%     properties other than those being perturbed are still preserved).
%     Properties not being perturbed by the new layer are linear
%     interpolated at the new depth points.
%
%     Multiple PROP/VALUE pairs may be passed to PERTURB_MODEL so that
%     multiple properties can be adjusted.  Note that the adjustments are
%     performed in the order that they are given so it is possible to
%     adjust a property and then readjust it several times.
%
%    Notes:
%
%    Examples:
%     PREM with a D" discontinuity in Vs.  There is a 100km thick layer
%     above the discontinuity at 2700km that decays from PREM to -1% of
%     PREM.  The discontinuity is a 3% jump.  Afterwards the D" layer
%     decays like an error function (that was given values from 0 to 2) to
%     the CMB which is -0.2km/s slower than the top of the D" layer:
%      mod=prem;
%      newmod=perturb_model(mod,'vs',[2600    0 1 1 10;
%                                     2700   -1 4 0  0;
%                                     2700    2 4 2 10;
%                                     2891 -0.2 1 0  0]);
%      figure; plot(newmod.depth,newmod.vs,'.');
%
%     Same but adding a 5-point, 25km thick ULVZ with a 10% drop in Vs:
%      mod=prem;
%      newmod=perturb_model(mod,'vs',[2600    0 1 1 10;
%                                     2700   -1 4 0  0;
%                                     2700    2 4 2 10;
%                                     2891 -0.2 1 0  0],...
%                               'vs',[2891-25 -10 4 1 5;
%                                     2891    -10 4 0 0]);
%      figure; plot(newmod.depth,newmod.vs,'.');
%
%    See also: PREM, AK135, IASP91, PREM_PERFECT, IASP91_PERFECT

%     Version History:
%        May  24, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  24, 2010 at 00:45 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));
if(mod(nargin-1,2))
    error('seizmo:perturb_model:unpairedOption',...
        'All Option/Values must be paired!');
end

% check model & properties
if(~isstruct(model) || ~isfield(model,'depth'))
    error('seizmo:perturb_model:badInput',...
        'MODEL appears not be a Earth model!');
elseif(~iscellstr(varargin(1:2:end)))
    error('seizmo:perturb_model:badInput',...
        'PROPERTY must be a string!');
elseif(any(~ismember(varargin(1:2:end),fieldnames(model))))
    error('seizmo:perturb_model:badInput',...
        'MODEL is missing 1 or more properties to be adjusted!');
end

% assemble model into a matrix
depths=model.depth;
ndep=numel(depths);
fields=fieldnames(model);
fields(strcmp(fields,'depth'))=[];
nf=numel(fields);
modmat=nan(ndep,nf);
for i=1:numel(fields)
    % could break if field doesn't have ndep entries
    modmat(:,i)=model.(fields{i});
end

% loop over prop/value pairs
for i=1:2:nargin-1
    % check value matrix is ok
    pidx=strcmp(fields,varargin{i});
    value=varargin{i+1};
    sz=size(value);
    if(~isreal(value))
        error('seizmo:perturb_model:badInput',...
            'VALUE matrix must be a real-valued Nx5 array!');
    elseif(numel(sz)~=2 || sz(2)~=5)
        error('seizmo:perturb_model:badInput',...
            'VALUE matrix must be [DEPTH VALUE VTYPE GTYPE NPTS]!');
    elseif(sz(1)<2 || any(value(end,4:5)~=0))
        error('seizmo:perturb_model:badInput',...
            ['VALUE matrix must have 2+ rows & ' ...
            'the last row must end with two 0s!']);
    elseif(any(value(:,1)<0 | value(:,1)>6371) || ~issorted(value(:,1)))
        error('seizmo:perturb_model:badInput',...
            'DEPTHS are out of range or out of order!');
    elseif(any(~ismember(value(:,3),0:4)) || any(value(:,5)<0) ...
            || any(value(:,5)~=fix(value(:,5))))
        error('seizmo:perturb_model:badInput',...
            'VTYPE must be 0 to 4 and NPTS must be positive integer!');
    end
    
    % get all discontinuities
    dci=find(diff(depths)==0);
    dcd=depths(dci);
    
    % pull model above and below
    dtop=value(1,1);
    dbot=value(end,1);
    if(any(dcd==dtop))
        topidx=dci(dcd==dtop);
    else
        topidx=find(depths<dtop,1,'last');
    end
    if(any(dcd==dbot))
        botidx=dci(dcd==dbot)+1;
    else
        botidx=find(depths>dbot,1);
    end
    modtop=modmat(1:topidx,:);
    modbot=modmat(botidx:end,:);
    deptop=depths(1:topidx);
    depbot=depths(botidx:end);
    
    % loop over each layer
    modnew=cell(sz(1)-1,1); depnew=modnew;
    v0=getvalue(depths,modmat,value(1,1),pidx,false);
    for j=1:sz(1)-1
        % depth range to adjust
        d1=value(j,1);
        d2=value(j+1,1);
        
        % values on either end
        switch value(j,3)
            case 0
                v1=value(j,2);
            case 1
                v1=value(j,2)+v0;
            case 2
                v1=(value(j,2)+100)/100*v0;
            case 3
                v1=value(j,2)...
                    +getvalue(depths,modmat,value(j,1),pidx,false);
            case 4
                v1=(value(j,2)+100)/100 ...
                    *getvalue(depths,modmat,value(j,1),pidx,false);
        end
        v0=v1;
        switch value(j+1,3)
            case 0
                v2=value(j+1,2);
            case 1
                v2=value(j+1,2)+v0;
            case 2
                v2=(value(j+1,2)+100)/100*v0;
            case 3
                v2=value(j+1,2)...
                    +getvalue(depths,modmat,value(j+1,1),pidx,true);
            case 4
                v2=(value(j+1,2)+100)/100 ...
                    *getvalue(depths,modmat,value(j+1,1),pidx,true);
        end
        v0=v2;
        
        % special handling of imposed discontinuity
        if(d1==d2)
            depnew{j}=d2;
            modnew{j}=nan(1,nf);
            modnew{j}(1,pidx)=v2;
            % this will give props at bottom of discon (what we want)
            modnew{j}(1,~pidx)=interpdc1(depths,modmat(:,~pidx),d2);
            continue;
        end
        
        % get new depths
        % - skip top if not first layer
        npts=value(j,5);
        if(j==1); begin=0; else begin=1; end
        dpts=(d1+(d2-d1)*(begin:npts-1)/(npts-1)).';
        
        % find discontinuities already existing in this range
        dcin=dcd>d1 & dcd<d2;
        dctop=dcd==d1;
        dcbot=dcd==d2;
        
        % if any new points are at the discontinuity, drop them
        dpts(ismember(dpts,dcd))=[];
        
        % combined discon + new depths
        dcds=[dcd(dctop); dcd(dcin); dcd(dcin); dcd(dcbot)];
        ndcds=numel(dcds);
        [depnew{j},lidx]=sort([dcds; dpts]);
        lidx(lidx)=1:numel(depnew{j});
        lidx1=lidx(1:ndcds);
        lidx2=lidx(ndcds+1:end);
        modnew{j}=nan(numel(depnew{j}),nf);
        
        % get adjusted property at all depths
        if(value(j,4)==0)
            % linear
            modnew{j}(:,pidx)=interp1q([d1; d2],[v1; v2],depnew{j});
        elseif(value(j,4)<0)
            % error function decaying from bottom
            x=-value(j,4)*(depnew{j}-d1)/(d2-d1);
            y=erf(x);
            y=y/max(y);
            modnew{j}(:,pidx)=v1+y*(v2-v1);
        else
            % error function decaying from top
            x=value(j,4)*(d2-depnew{j})/(d2-d1);
            y=erf(x);
            y=y/max(y);
            modnew{j}(:,pidx)=v2+y*(v1-v2);
        end
        
        % get other properties at new depths
        modnew{j}(lidx2,~pidx)=interpdc1(depths,modmat(:,~pidx),dpts);
        
        % get other properties at discons
        dcis=[dci(dctop)+1; dci(dcin); dci(dcin)+1; dci(dcbot)];
        modnew{j}(lidx1,~pidx)=modmat(dcis,~pidx);
    end
    
    % now reassemble
    modmat=[modtop; cat(1,modnew{:}); modbot];
    depths=[deptop; cat(1,depnew{:}); depbot];
end

% back into struct
model.depth=depths;
for i=1:nf
    model.(fields{i})=modmat(:,i);
end

end

function v=getvalue(depths,model,depth,col,top)
% value is known?
if(any(depth==depths))
    % yes, get it
    didx=depth==depths;
    if(sum(didx)==2)
        if(top)
            didx=find(didx,1);
        else
            didx=find(didx,1,'last');
        end
    end
    v=model(didx,col);
else
    % no, interpolate it
    dup=find(depths<depth,1,'last');
    v=interp1q(depths([dup dup+1]),model([dup dup+1],col),depth);
end
end
