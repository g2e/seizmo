function [corr,ttmod,ttref]=crucor(lat,lon,rayp,wtype,varargin)
%CRUCOR    Returns crustal travel time corrections
%
%    Usage:    corr=crucor(lat,lon,rayp,wavetype)
%              corr=crucor(...,'elev',elev,...)
%              corr=crucor(...,'hole',depth,...)
%              corr=crucor(...,'stretch',logical,...)
%              corr=crucor(...,'model',crustal_model,...)
%              corr=crucor(...,'refmod',1dmodel,...)
%              corr=crucor(...,'bottom',1ddepth,...)
%              corr=crucor(...,'top',1ddepth,...)
%              [corr,ttmod,ttref]=crucor(...)
%
%    Description:
%     CORR=CRUCOR(LAT,LON,RAYP,WAVETYPE) calculates the difference in
%     travel time through the default 3D crustal model (see the MODEL
%     option for more info) compared to PREM for a seismic phase with ray
%     parameter RAYP (in sec/deg) and of type WAVETYPE ('p' or 's').  The
%     travel time through the crustal model goes from the moho to the
%     elevation at LAT/LON for SRTM30_PLUS.  To account for the difference
%     between the crustal model's elevation and SRTM30_PLUS, the model is
%     stretched or squished to match while preserving the moho depth.  The
%     travel time for PREM is for a section going from sealevel (PREM has
%     no topography) to the moho depth of the crustal model.  The
%     correction CORR is in seconds and gives TT3D=TT1D+CORR.  Combining
%     CRUCOR with ELLCOR and MANCOR will provide a more complete 3D
%     travel time correction.
%
%     CORR=CRUCOR(...,'ELEV',ELEV,...) sets the local elevations used in
%     the crustal corrections to ELEV if and only if STRETCH=TRUE.  If
%     STRETCH is set to TRUE (the default - see below) then the crust model
%     layers are stretched to match this elevation while preserving the
%     depth to the moho discontinuity.  By default, if STRETCH=TRUE then
%     ELEV is the topography from SRTM30_PLUS (via TOPO_POINTS).  If
%     STRETCH=FALSE, ELEV is the topography from the crustal model and
%     cannot be changed (ELEV option is ignored).  Units are in km.
%
%     CORR=CRUCOR(...,'HOLE',DEPTH,...) sets the seismometer depth from the
%     local elevation.  Units are in km.  Default is 0.
%
%     CORR=CRUCOR(...,'STRETCH',LOGICAL,...) indicates if the crust is to
%     be stretched to the elevation value supplied through ELEV while
%     preserving the moho discontinuity depth.  If STRETCH is false, ELEV
%     is ignored (it is replaced by the crustal model elevation).  Default
%     is TRUE.
%
%     CORR=CRUCOR(...,'MODEL',CRUSTAL_MODEL,...) sets the crustal model
%     used to generate the crustal travel time corrections from the
%     reference 1D model.  The default model and available models are
%     defined by GETCRUST.
%
%     CORR=CRUCOR(...,'REFMOD',1DMODEL,...) sets the reference 1D Earth
%     model to correct.  Can be any listed by AVAILABLE_1DMODELS.  Default
%     is 'PREM'.
%
%     CORR=CRUCOR(...,'BOTTOM',1DDEPTH,...)
%     CORR=CRUCOR(...,'TOP',1DDEPTH,...) are for setting bounds on the
%     depths of the models used in the crustal correction.  This is for
%     source side corrections where the earthquake is located within the
%     crust.  It is currently an error to set both HOLE & TOP to nonzero.
%     The default is 0 for TOP and the moho depth for BOTTOM.
%
%     [CORR,TTMOD,TTREF]=CRUCOR(...) also returns the total travel time
%     through both the crustal model (TTMOD) and the 1D reference model
%     (TTREF).
%
%    Notes:
%     - Latitudes are assumed to be geographic.
%
%    Examples:
%     % Generate a map of crustal corrections for Pdiff (no distance
%     % dependance to Pdiff ray parameter which allows us to do this):
%     rayp=4.428;
%     [lon,lat]=meshgrid(-179.5:179.5,-89.5:89.5);
%     corr=crucor(lat,lon,rayp,'p','s',false);
%     figure; imagesc(-179.5:179.5,-89.5:89.5,reshape(corr,size(lon)));
%     axis xy equal tight;
%     xlabel('Longitude')
%     ylabel('Latitude')
%     title('Pdiff Crustal Travel Time Corrections');
%     hc=colorbar; ylabel(hc,'sec');
%
%    See also: ELLCOR, MANCOR, GETCRUST, PREM, AK135, IASP91, TOPO_POINTS

%     Version History:
%        May  20, 2010 - initial version
%        Jan. 14, 2011 - improved verbose message, fixed message bug
%        Apr.  2, 2012 - minor doc update
%        Aug.  6, 2012 - better handling of stretch=false, doc update
%        Jan. 23, 2014 - update for crustal model changes, model option
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 23, 2014 at 15:45 GMT

% todo

% check nargin
error(nargchk(4,inf,nargin));

% constants
d2r=pi/180;

% check required inputs
if(ischar(wtype)); wtype=cellstr(wtype); end
if(~isreal(lat) || ~isreal(lon) || ~isreal(rayp))
    error('seizmo:crucor:badInput',...
        'LAT, LON, & RAYP must be real valued arrays!');
elseif(~iscellstr(wtype) || any(~ismember(lower(wtype),{'p' 's'})))
    error('seizmo:crucor:badInput',...
        'WAVETYPE must be ''P'' or ''S''!');
elseif(~isequalsizeorscalar(lat,lon,rayp,wtype))
    error('seizmo:crucor:badInput',...
        'LAT, LON, RAYP, & WAVETYPE should be scalar or equal sized!');
end

% expand scalars
[lat,lon,rayp,wtype]=expandscalars(lat,lon,rayp,wtype);
sz=size(lat); npts=prod(sz);

% valid 1d reference models
valid.refmod=available_1dmodels;

% option defaults
varargin=[{'e' [] 'h' 0 't' 0 'b' [] 's' true 'm' [] 'r' 'prem'} varargin];

% go through optional inputs
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:crucor:badOption',...
        'All Options must be specified with a string!');
end
for i=1:2:numel(varargin)
    % skip empty
    skip=false;
    if(isempty(varargin{i+1})); skip=true; end
    
    % check option is available
    switch lower(varargin{i})
        case {'e' 'elev' 'stel'}
            % local elevation in km
            % - not allowing anything outside +/-12km (m vs km check)
            if(~isempty(varargin{i+1}) && (~isreal(varargin{i+1}) ...
                    || any(abs(varargin{i+1})>12)))
                error('seizmo:crucor:badELEV',...
                    'ELEV must be a real-valued array in km!');
            end
            elev=varargin{i+1};
        case {'h' 'hole' 'stdp'}
            if(skip); continue; end
            % depth of seismometer
            % - relative to surface in km
            % - values < 0 (airborne) are set to 0
            % - not allowing > 10km (m vs km check)
            if(~isreal(varargin{i+1}) || any(abs(varargin{i+1})>10))
                error('seizmo:crucor:badTOP',...
                    'HOLE must be a real-valued seismometer depth in km!');
            end
            hole=varargin{i+1};
        case {'t' 'top'}
            if(skip); continue; end
            % depth of event/bounce in km
            % - relative to sea level in km
            % - not allowing > 100km (m vs km check)
            if(~isreal(varargin{i+1}) || any(abs(varargin{i+1})>100))
                error('seizmo:crucor:badTOP',...
                    'TOP must be a real-valued array in km!');
            end
            top=varargin{i+1};
        case {'b' 'bot' 'bottom'}
            % bottom of crust to use in correction
            % - relative to sealevel in km
            % - values > local moho are set to moho
            % - not allowing > 100km moho (m vs km check)
            if(~isempty(varargin{i+1}) && (~isreal(varargin{i+1}) ...
                    || any(abs(varargin{i+1})>100)))
                error('seizmo:crucor:badBOTTOM',...
                    'BOTTOM must be a real-valued array in km!');
            end
            bot=varargin{i+1};
        case {'s' 'str' 'stretch'}
            if(skip); continue; end
            % stretch crust to local elevation (set via elev)
            % - true or false
            if(~islogical(varargin{i+1}))
                error('seizmo:crucor:badSTRETCH',...
                    'STRETCH must be a logical array!');
            end
            stretch=varargin{i+1};
        case {'m' 'cm' 'mod' 'cmod' 'crustal_model' 'model'}
            % 3D crustal model
            if(~isempty(varargin{i+1}) && (~ischar(varargin{i+1}) ...
                    || ndims(varargin{i+1})~=2 ...
                    || size(varargin{i+1},1)~=1))
                error('seizmo:crucor:badMOD',...
                    'MODEL must be a string!');
            end
            model=varargin{i+1};
        case {'r' 'rm' 'ref' 'refmod'}
            if(skip); continue; end
            % reference model to calc against
            % should call a function that has a list
            if(ischar(varargin{i+1}))
                varargin{i+1}=cellstr(varargin{i+1});
            end
            if(~iscellstr(varargin{i+1}) ...
                    || any(~ismember(varargin{i+1},valid.refmod)))
                error('seizmo:crucor:badREFMOD',...
                    ['REFMOD must be one of the following:\n' ...
                    sprintf('''%s'' ',valid.refmod{:})]);
            end
            refmod=varargin{i+1};
    end
end

% get crust model
crust=getcrust(geographic2geocentriclat(lat),lon,model);

% extract moho depth & elevation
moho=reshape(-crust.top(:,end),sz);
celev=reshape(crust.top(:,2),sz);

% set bottom to moho if empty
if(isempty(bot)); bot=moho; end

% set elev (if not set) to that of the crustal model
eflag=false;
if(isempty(elev))
    elev=celev;
    if(isscalar(stretch))
        eflag=stretch;
    else
        eflag=true;
    end
end

% check sizes & expand
if(~isequalsizeorscalar(lat,elev,hole,top,bot,stretch,refmod))
    error('seizmo:crucor:badInput',...
        'All inputs must be scalar or equal sized!');
end
[elev,hole,top,bot,stretch,refmod]=...
    expandscalars(elev,hole,top,bot,stretch,refmod,lat);

% get highres local elevation if stretching and no elevation given
% - this is really slow
if(eflag); elev(stretch)=topo_points(lat(stretch),lon(stretch))/1000; end

% force crustal model elevation if not stretching
% - already set if eflag
if(~eflag); elev(~stretch)=celev(~stretch); end

% do not allow hole & top ~= 0
if(any(hole & top))
    error('seizmo:crucor:badInput',...
        'HOLE & TOP can not be non-zero at the same time!');
end

% set hole to be relative to sea level
hole(hole<0)=0; % no airborne
hole=hole-elev;

% keep track of source side
ss=top~=0;

% replace out-of-range values with more valid ones
bot(bot>moho)=moho(bot>moho); % cannot go below moho
top(top>bot)=bot(top>bot); % we will shortcut for top==bot, corr=0
bot(bot<top)=top(bot<top); % we will shortcut for top==bot, corr=0
top=max(top,-elev); % no airborne event/bounce

% set top to hole if receiver side
top(~ss)=hole(~ss);

% convert to radii
rtop=6371-top;
rbot=6371-bot;

% convert ray parameter from sec/deg to sec/radian
rayp=rayp/d2r;

% get reference models ahead of time
[mod1dnames,refidx,refidx]=unique(refmod);
ref=cell(numel(mod1dnames),1);
for i=1:numel(mod1dnames)
    modfun=str2func(mod1dnames{i});
    ref{i}=modfun();
end

% verbosity
verbose=seizmoverbose;
if(verbose)
    disp(['Getting Crustal Correction(s) for ' crust.model ' vs ' ...
        upper(joinwords(mod1dnames,' & '))]);
    print_time_left(0,npts);
end

% loop over points
ttmod=zeros(npts,1); ttref=ttmod;
for i=1:npts
    % special skip if top==bot
    if(top(i)==bot(i)); continue; end
    
    % stretch crustal model if wanted
    if(stretch(i))
        factor=(moho(i)+elev(i))/(moho(i)+celev(i));
        crust.thk(i,:)=crust.thk(i,:)*factor;
    end
    
    % get radii of layers
    rlbot=6371+elev(i)-cumsum(crust.thk(i,2:end));
    rltop=[6371+elev(i) rlbot(1:end-1)];
    %rlmid=(rlbot+rltop)/2;
    
    % get fractions of each layer that are in correction
    frac=max(min(rltop,rtop(i))-max(rlbot,rbot(i)),0)...
        ./crust.thk(i,2:end);
    frac(crust.thk(i,2:end)==0)=0;
    
    % velocity based on wtype
    switch lower(wtype{i})
        case 'p'
            v=crust.vp(i,2:end-1);
            v(v==0)=1; % avoid divide by zero
        case 's'
            v=crust.vs(i,2:end-1);
            v(v==0)=1; % avoid divide by zero
    end
    
    % exact formula for tt in constant velocity spherical shells
    % sqrt(r ^2/v^2-p^2)-sqrt(r ^2/v^2-p^2)
    %       1                  2
    %
    % r  >  r
    %  1     2
    ttmod(i)=sum(frac.*(sqrt((rltop./v).^2-rayp(i)^2) ...
                      -sqrt((rlbot./v).^2-rayp(i)^2)));
    
    % ALTERNATIVE
    % numerically evaluate integral using adaptive Simpson quadrature
    % (for tt in model with constant velocity)
    % - using equation (16) on page 159 in Stein & Wysession
    % - requires loop
    %for j=1:numel(v)
    %    ttmod(i)=ttmod(i)+quad(...
    %        @(r)(r/v(j)).^2./(r.*sqrt((r/v(j)).^2-rayp(i)^2)),...
    %        rlbot(j),rltop(j));
    %end
    
    % ALTERNATIVE
    % my old approximate formula for tt
    % - get path length in layer based on rayparameter
    %len=crust.thk(i,2:end)./cos(asin(rayp(i)*v./rlmid));
    %ttmod(i)=sum(frac.*len./v);
    
    % allocate reference model
    radii=6371-ref{refidx(i)}.depth;
    switch lower(wtype{i})
        case 'p'
            v=ref{refidx(i)}.vp;
        case 's'
            v=ref{refidx(i)}.vs;
    end
    
    % clipping reference model
    % - don't clip top if receiver side
    if(ss(i)) % source side
        % knots that are definitely included
        ktop=find(radii<rtop(i),1);
        kbot=find(radii>rbot(i),1,'last');
        
        % do we need to interpolate top/bottom knot
        tf=ismember([rtop(i) rbot(i)],radii);
        if(tf(1))
            % no, it is already there
            ktop=ktop-1;
        else
            % yes, get it
            vtop=interp1q(radii(ktop:-1:ktop-1),v(ktop:-1:ktop-1),rtop(i));
        end
        if(tf(2))
            % no, it is already there
            kbot=kbot+1;
        else
            % yes, get it
            vbot=interp1q(radii(kbot+1:-1:kbot),v(kbot+1:-1:kbot),rbot(i));
        end
        
        % clip model
        radii=radii(ktop:kbot);
        v=v(ktop:kbot);
        
        % pad on top/bottom if not there
        if(~tf(1))
            radii=[rtop(i); radii]; %#ok mlint is stupid
            v=[vtop; v]; %#ok mlint is stupid
        end
        if(~tf(2))
            radii=[radii; rbot(i)]; %#ok mlint is stupid
            v=[v; vbot]; %#ok mlint is stupid
        end
    else % receiver side
        % knots that are definitely included
        ktop=1;
        kbot=find(radii>rbot(i),1,'last');
        
        % do we need to interpolate bottom knot
        tf=any(rbot(i)==radii);
        if(tf)
            % no, it is already there
            kbot=kbot+1;
        else
            % yes, get it
            vbot=interp1q(radii(kbot+1:-1:kbot),v(kbot+1:-1:kbot),rbot(i));
        end
        
        % clip model
        radii=radii(ktop:kbot);
        v=v(ktop:kbot);
        
        % pad on bottom if not there
        if(~tf)
            radii=[radii; rbot(i)]; %#ok mlint is stupid
            v=[v; vbot]; %#ok mlint is stupid
        end
    end
    
    % TT CALC FOR MODEL WITH LINEAR VELOCITY GRADIENT
    % numerically evaluate integral using adaptive Simpson quadrature
    % - using equation (16) on page 159 in Stein & Wysession
    % - requires loop over each section between knots
    % - we assume the model is linear between knots
    %   v(r)=a*r+b
    % - we could use the actual model polynomials instead ....
    nlay=numel(radii)-1;
    a=diff(v)./diff(radii);
    b=v(1:nlay)-a.*radii(1:nlay);
    good=find(diff(radii))';
    for j=good
        ttref(i)=ttref(i)+quad(...
            @(r)(r./(a(j)*r+b(j))).^2./(r.*sqrt((r./(a(j)*r+b(j))).^2-rayp(i)^2)),...
            radii(j+1),radii(j));
    end
    
    % detail message
    if(verbose); print_time_left(i,npts); end
end

% corrections
corr=ttmod-ttref;

end
