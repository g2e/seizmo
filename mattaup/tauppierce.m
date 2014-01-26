function [varargout]=tauppierce(varargin)
%TAUPPIERCE    Calculate ray path pierce points using the TauP toolkit
%
%    Usage:    tauppierce(...)
%           tt=tauppierce(...)
%           tt=tauppierce(...,'m|mod|model',model,...)
%           tt=tauppierce(...,'h|z|dep|evdp|depth',depth,...)
%           tt=tauppierce(...,'p|ph|phases',phases,...)
%           tt=tauppierce(...,'d|deg|gcarc|degrees',degdist,...)
%           tt=tauppierce(...,'k|km|dist|kilometers',kmdist,...)
%           tt=tauppierce(...,'s|st|sta|station',[lat lon],...)
%           tt=tauppierce(...,'e|ev|evt|event',[lat lon],...)
%           tt=tauppierce(...,'a|az|azi|azimuth',azi,...)
%           tt=tauppierce(...,'b|baz|bazi|backazimuth',bazi,...)
%           tt=tauppierce(...,'t|turn|turning',true|false,...)
%           tt=tauppierce(...,'u|under|underside',true|false,...)
%           tt=tauppierce(...,'r|rev|reverse',true|false,...)
%           tt=tauppierce(...,'n|nd|nodiscon',true|false,...)
%           tt=tauppierce(...,'x|pnts|pierce',depths,...)
%
%    Description:
%     TAUPPIERCE(...) (no outputs, w/ or w/o inputs) displays a formatted
%     list of information on many seismic phases found at a random distance
%     from an event with a depth of 0km.  This is the same as for TAUPTIME
%     (see the second form for the difference between TAUPTIME and
%     TAUPPIERCE).  The 1D Earth model utilized by default is IASP91 (see
%     the MODEL option to adjust this).  The default phase list is
%     'ttbasic' (equivalent to setting phases to BASIC in Brian Kennett's
%     TTIMES program) which lists many common phases by default.  See
%     option PHASES to adjust the phase list.  A variety of options
%     (DEG, KM, STA & EVT) allow for changing the distance from the event.
%
%     TT=TAUPPIERCE(...) (1 output, w/ or w/o inputs) returns a structure
%     array with the following layout:
%       TT(index).modelname        - name of velocity model
%                .event            - event location if known ([lat lon])
%                .station          - stn location if known ([lat lon])
%                .depth            - depth of earthquake (km)
%                .distance         - true distance (deg)
%                .mindistance      - least distance (deg)
%                .phase            - seismic phase name
%                .puristphase      - verbose seismic phase name
%                .time             - travel time (sec)
%                .rayparameter     - ray parameter (sec/deg)
%                .pierce.time      - travel times (sec) at pierce points
%                .pierce.distance  - distance (deg) of pierce points
%                .pierce.depth     - depth (km) of pierce points
%                .pierce.latitude  - latitude (deg) of pierce points
%                .pierce.longitude - longitude (deg) of pierce points
%     Each phase has its own indice in the struct array TT.  Use TT(index)
%     to access individual phase information.  The sub-structure 'pierce'
%     contains several fields populated with arrays of info on the pierce
%     points for the phase.
%
%     *********************************************************
%     All the following field/value pair options may be entered
%     at the command line in any order.
%     *********************************************************
%
%     TT=TAUPPIERCE(...,'M|MOD|MODEL',MODEL,...) sets the 1D Earth model to
%     MODEL.  Accepts a variety of common models like 'prem', 'iasp91',
%     'ak135'.  See the TauP program/documentation for more.  You may also
%     give a 1DMODEL struct like from CMB_1DMODEL_LIBRARY.  The default
%     model is 'iasp91'.
%
%     TT=TAUPPIERCE(...,'H|Z|DEP|EVDP|DEPTH',DEPTH,...) sets the event
%     depth to DEPTH.  DEPTH should be in kilometers!  The default is 0.
%
%     TT=TAUPPIERCE(...,'P|PH|PHASES',PHASES,...) sets the phase list to
%     PHASES.  PHASES must be a comma separated list of seismic phases.  
%     See the TauP documentation for conventions and valid phase names.
%     The default phase list is 'ttbasic'.
%
%     TT=TAUPPIERCE(...,'D|DEG|GCARC|DEGREES',DEGDIST,...) sets the event
%     distance to DEGDIST.  DEGDIST is expected to be in angular degrees.
%     The default value is random (somewhere between 0 and 180 degrees).
%     DEGDIST must be scalar.
%
%     TT=TAUPPIERCE(...,'K|KM|DIST|KILOMETERS',KMDIST,...) sets the event
%     distance to KMDIST.  KMDIST is expected to be in kilometers.  There
%     is no default value for this option.  KMDIST must be scalar.
%
%     TT=TAUPPIERCE(...,'S|ST|STA|STATION',[LAT LON],...) sets the seismic
%     station location to [LAT LON].  LAT and LON must be in degrees.  This
%     option must be paired with the EVT option to calculate the event
%     distance.  There is no default value for this option.
%
%     TT=TAUPPIERCE(...,'E|EV|EVT|EVENT',[LAT LON],...) sets the earthquake
%     location to [LAT LON].  LAT and LON must be in degrees.  This option
%     must be paired with the STA option to calculate the event distance.
%     There is no default value for this option.
%
%     TT=TAUPPIERCE(...,'A|AZ|AZI|AZIMUTH',AZI,...) sets the azimuth from
%     the earthquake to the seismic station.  Note that this option really
%     isn't useful unless paired with the EVT option and one of the DEG or
%     KM options so that the station location can be found.  There is no
%     default value for this option.
%
%     TT=TAUPPIERCE(...,'B|BAZ|BAZI|BACKAZIMUTH',BAZI,...) sets the
%     backazimuth from the seismic station to the earthquake.  This option
%     needs to be paired with one of the DEG and KM options plus the STA
%     option to allow for calculating the earthquake location.  There is no
%     default value for this option.
%
%     TT=TAUPPIERCE(...,'T|TURN|TURNING',TRUE|FALSE,...) setting TURN to
%     TRUE will just return turning (aka bottoming) points.  This option
%     may be combined with UNDER, REV, and NODISCON (ie they do not exclude
%     one another - although REV is a super-set of TURN).  By default TURN
%     is set to FALSE.
%
%     TT=TAUPPIERCE(...,'U|UNDER|UNDERSIDE',TRUE|FALSE,...) toggles whether
%     or not to just return underside reflection points.  This option may
%     be combined with the TURN, REV, and NODISCON options (no mutual
%     exclusion).  By default UNDER is set to FALSE.
%
%     TT=TAUPPIERCE(...,'R|REV|REVERSE',TRUE|FALSE,...) toggles the
%     returning of only points where a reversal occurs (turning and
%     reflection points).  By default REV is set to FALSE.
%
%     TT=TAUPPIERCE(...,'N|ND|NODISCON',TRUE|FALSE,...) toggles returning
%     just the points in the PIERCE option.  May be combined with the TURN,
%     UNDER, and REV options without excluding each others results.  By
%     default NODISCON is set to FALSE.
%
%     TT=TAUPPIERCE(...,'X|PNTS|PIERCE',DEPTHS,...) adds the depths in
%     DEPTHS to the ray path pierce points to calculate.  Depths must be a
%     numeric array. There is no default value for this option.
%
%    Notes:
%     - TauP toolkit developed by:
%        H. Philip Crotwell, Thomas J. Owens, Jeroen Ritsema
%        Department of Geological Sciences
%        University of South Carolina
%        http://www.seis.sc.edu
%        crotwell@seis.sc.edu
%
%    Examples:
%     % P'P' pierce points for a 44km deep event at some random distance:
%     tauppierce('p','PKIKPPKIKP','z',44)
%
%     % Add some other interesting piercing depths:
%     tauppierce('p','PKIKPPKIKP','z',44,'x',[1800 2580])
%
%     % Return underside reflection points for pS^640S (note 3 nonzeros):
%     tauppierce('p','pS^640S','z',33,'d',120,'under',true)
%
%    See also: TAUPPATH, TAUPCURVE, TAUPTIME, TAUPCREATE, TAUP

%     Version History:
%        Sep.  2, 2009 - major revision of script, name change to avoid
%                        breakage due to input/output changes
%        Sep.  5, 2009 - minor doc update
%        Sep. 30, 2009 - changed abssawmod to abslatmod
%        Nov. 13, 2009 - dropped abslatmod for getModuloDistDeg, dropped
%                        some import calls
%        Jan.  6, 2011 - add matTaup.jar to dynamic java classpath if
%                        necessary
%        Feb. 24, 2012 - switch to native taup
%        Jan. 26, 2014 - minor fix necessary for update to TauP 2.1.1, no
%                        longer need to update jar filenames
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 26, 2014 at 17:15 GMT

% todo:

% check nargin
if(mod(nargin,2))
    error('TauP:tauppierce:badNumOptions','Unpaired option(s)!');
end

% try adding *.jar if no TauP class exists
if(~exist('edu.sc.seis.TauP.TauP_Pierce','class'))
    fs=filesep;
    mypath=fileparts(mfilename('fullpath'));
    jars=dir([mypath fs 'lib' fs '*.jar']);
    for i=1:numel(jars)
        if(~ismember([mypath fs 'lib' fs jars(i).name],javaclasspath))
            javaaddpath([mypath fs 'lib' fs jars(i).name]);
        end
    end
end

% default options
model='iasp91'; modelname='iasp91';
depth=0; phases={'ttbasic'};
deg=[]; km=[];
st=[]; ev=[];
az=[]; baz=[];
turn=false;
under=false;
rev=false;
nodiscon=false;
pierce=[]; spierce='';

% check options
if(~iscellstr(varargin(1:2:end)))
    error('TauP:tauppierce:badInput',...
        'All options must be specified with strings!');
end
for i=1:2:nargin
    if(isempty(varargin{i+1})); continue; end
    switch lower(varargin{i})
        case {'m' 'mod' 'model'}
            if(ischar(varargin{i+1})) % name of model on TauP path
                model=varargin{i+1};
                modelname=model;
            elseif(isempty(chk1dmodel(varargin{i+1}))) % 1dmodel struct
                model=taupcreate(varargin{i+1});
                modelname=char(model.getModelName);
            elseif(strcmp(class(varargin{i+1}),... % velocity model obj
                    'edu.sc.seis.TauP.VelocityModel'))
                model=taupcreate(varargin{i+1});
                modelname=char(model.getModelName);
            elseif(strcmp(class(varargin{i+1}),... % tau model obj
                    'edu.sc.seis.TauP.TauModel'))
                model=varargin{i+1};
                modelname=char(model.getModelName);
            else
                error('TauP:tauppierce:badInput',...
                    'MODEL must be a string (like ''prem'')!');
            end
        case {'h' 'z' 'dep' 'evdp' 'depth'}
            if(~isscalar(varargin{i+1}) ...
                    || ~isreal(varargin{i+1}) || varargin{i+1}<0)
                error('TauP:tauppierce:badInput',...
                    'DEPTH must be a positive number (in km)!');
            end
            depth=varargin{i+1};
        case {'p' 'ph' 'phases'}
            if(~ischar(varargin{i+1}))
                error('TauP:tauppierce:badInput',...
                    'PHASES must be a string (like ''P,S'')!');
            end
            phases=varargin{i+1};
            phases=strtrim(getwords(phases,','));
        case {'d' 'deg' 'gcarc' 'degrees'}
            if(~isscalar(varargin{i+1}) ...
                    || ~isreal(varargin{i+1}) || varargin{i+1}<0)
                error('TauP:tauppierce:badInput',...
                    'DEG must be a positive number (in degrees)!');
            end
            deg=varargin{i+1};
        case {'k' 'km' 'dist' 'kilometers'}
            if(~isscalar(varargin{i+1}) ...
                    || ~isreal(varargin{i+1}) || varargin{i+1}<0)
                error('TauP:tauppierce:badInput',...
                    'KM must be a positive number (in kilometers)!');
            end
            km=varargin{i+1};
        case {'s' 'st' 'sta' 'station'}
            if(~isequal([1 2],size(varargin{i+1})) ...
                    || ~all(isreal(varargin{i+1})))
                error('TauP:tauppierce:badInput',...
                    'STA must be 2 real numbers ([LAT LON] in degrees)!');
            end
            st=varargin{i+1};
        case {'e' 'ev' 'evt' 'event'}
            if(~isequal([1 2],size(varargin{i+1})) ...
                    || ~all(isreal(varargin{i+1})))
                error('TauP:tauppierce:badInput',...
                    'EVT must be 2 real numbers ([LAT LON] in degrees)!');
            end
            ev=varargin{i+1};
        case {'a' 'az' 'azi' 'azimuth'}
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}))
                error('TauP:tauppierce:badInput',...
                    'AZ must be a scalar number (in degrees)!');
            end
            az=varargin{i+1};
        case {'b' 'baz' 'bazi' 'backazimuth'}
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}))
                error('TauP:tauppierce:badInput',...
                    'BAZ must be a scalar number (in degrees)!');
            end
            baz=varargin{i+1};
        case {'t' 'turn' 'turning'}
            if(~isscalar(varargin{i+1}) || ~islogical(varargin{i+1}))
                error('TauP:tauppierce:badInput',...
                    'TURN must be a scalar logical!');
            end
            turn=varargin{i+1};
        case {'u' 'under' 'underside'}
            if(~isscalar(varargin{i+1}) || ~islogical(varargin{i+1}))
                error('TauP:tauppierce:badInput',...
                    'UNDER must be a scalar logical!');
            end
            under=varargin{i+1};
        case {'r' 'rev' 'reverse'}
            if(~isscalar(varargin{i+1}) || ~islogical(varargin{i+1}))
                error('TauP:tauppierce:badInput',...
                    'REV must be a scalar logical!');
            end
            rev=varargin{i+1};
        case {'n' 'nd' 'nodiscon'}
            if(~isscalar(varargin{i+1}) || ~islogical(varargin{i+1}))
                error('TauP:tauppierce:badInput',...
                    'NODISCON must be a scalar logical!');
            end
            nodiscon=varargin{i+1};
        case {'x' 'pnts' 'pierce'}
            if(any(~isreal(varargin{i+1}(:))) || any(varargin{i+1}(:)<0))
                error('TauP:tauppierce:badInput',...
                    'PIERCE must be positive depths in kilometers!');
            end
            pierce=varargin{i+1}(:);
            spierce=strcat(num2str(pierce),',')'; % comma delimit
            spierce=spierce(:)';                  % make row vector
            spierce=spierce(spierce~=32);         % remove spaces
            spierce=spierce(1:end-1);             % trim off trailing comma
        otherwise
            error('TauP:tauppierce:badOption',...
                'Unknown Option: %s',varargin{i});
    end
end

% check geometry
sc=javaObject('edu.sc.seis.TauP.SphericalCoords');
d=~isempty(deg); k=~isempty(km); s=~isempty(st);
e=~isempty(ev); a=~isempty(az); b=~isempty(baz);
if(d && k)
    % over specified
    error('TauP:tauppierce:badInput',...
        'DEG & KM options cannot both be used!');
elseif(k)
    deg=km/(6371*pi/180);
elseif(d)
    %km=deg*(6371*pi/180);
elseif(~(e && s))
    % default to random distance
    disp('Distance not given or could not be determined - using random!')
    deg=rand*180;
    %km=deg*(6371*pi/180);
    d=true;
end
if((e && s) || (e && a && (d || k)) || (s && b && (d || k)))
    % specific geometry
    generic=false;
    if(e && s)
        deg=sc.distance(ev(1),ev(2),st(1),st(2));
        az=sc.azimuth(ev(1),ev(2),st(1),st(2));
        %baz=sc.azimuth(st(1),st(2),ev(1),ev(2));
    elseif(s)
        ev(1)=sc.latFor(st(1),st(2),deg,baz);
        ev(2)=sc.lonFor(st(1),st(2),deg,baz);
        az=sc.azimuth(ev(1),ev(2),st(1),st(2));
    elseif(e)
        st(1)=sc.latFor(ev(1),ev(2),deg,az);
        st(2)=sc.lonFor(ev(1),ev(2),deg,az);
        %baz=sc.azimuth(st(1),st(2),ev(1),ev(2));
    end
elseif((d || k) && ~((d && k) || s || e || a || b))
    % generic geometry
    generic=true;
else
    error('TauP:tauppierce:badInput',...
        'Event-Station geometry under/over specified!');
end

% debug
%[generic ev st deg az baz]

% create path object for velocity model
txobj=javaObject('edu.sc.seis.TauP.TauP_Pierce',model);

% calculate pierce points for the specified parameters
txobj.setSourceDepth(depth);
txobj.setPhaseNames(phases);
%txobj.setOnlyTurnPoints(turn);    % These do not do
%txobj.setOnlyUnderPoints(under);  % anything useful
%txobj.setOnlyRevPoints(rev);      % here so I account
%txobj.setOnlyAddPoints(nodiscon); % for them below. 
if(~isempty(pierce)); txobj.setAddDepths(spierce); end
txobj.calculate(deg);
narr=txobj.getNumArrivals;

% conversion
R2D=180/pi;

% struct output
tt(1:narr,1)=struct('modelname',[],'event',[],'station',[],'depth',[],...
    'distance',[],'mindistance',[],'phase',[],'puristphase',[],...
    'time',[],'rayparameter',[],'pierce',[]);
for ii=1:narr
    % get phase info
    arr=txobj.getArrival(ii-1);
    tt(ii).modelname=modelname;
    tt(ii).event=ev;
    tt(ii).station=st;
    tt(ii).depth=arr.getSourceDepth;
    tt(ii).distance=arr.getDistDeg;
    tt(ii).mindistance=arr.getModuloDistDeg;
    tt(ii).phase=char(arr.getName);
    tt(ii).puristphase=char(arr.getPuristName);
    tt(ii).time=arr.getTime;
    tt(ii).rayparameter=arr.getRayParam/R2D;
    
    % pierce info (preallocate)
    pts=arr.getPierce; npts=numel(pts);
    tt(ii).pierce.depth=nan(npts,1);
    tt(ii).pierce.distance=nan(npts,1);
    tt(ii).pierce.time=nan(npts,1);
    if(~generic)
        tt(ii).pierce.latitude=nan(npts,1);
        tt(ii).pierce.longitude=nan(npts,1);
    end
    for jj=1:npts
        % this gets all points
        tt(ii).pierce.depth(jj)=pts(jj).depth;
        %tt(ii).pierce.distance(jj)=pts(jj).dist*R2D;      % For  < 2.X
        tt(ii).pierce.distance(jj)=pts(jj).distRadian*R2D; % For >= 2.X
        tt(ii).pierce.time(jj)=pts(jj).time;
        if(~generic)
            tt(ii).pierce.latitude(jj)=sc.latFor(ev(1),ev(2),...
                tt(ii).pierce.distance(jj),az);
            tt(ii).pierce.longitude(jj)=sc.lonFor(ev(1),ev(2),...
                tt(ii).pierce.distance(jj),az);
        end
    end
    
    % limiting output
    if(turn || under || rev || nodiscon)
        % what are the points?
        t=[false; diff(tt(ii).pierce.depth)>=0 ...
            & [diff(tt(ii).pierce.depth(2:end)); 0]<=0];
        u=[0; diff(tt(ii).pierce.depth)]<=0 ...
            & [diff(tt(ii).pierce.depth); 0]>=0;
        r=[0; diff(tt(ii).pierce.depth)].*[diff(tt(ii).pierce.depth); 0]<0;
        n=ismember(tt(ii).pierce.depth,pierce);
        %[t u r n]
        
        % isolate
        ok=false(npts,1);
        if(turn); ok=ok | t; end
        if(under); ok=ok | u; end
        if(rev); ok=ok | r; end
        if(nodiscon); ok=ok | n; end
        nok=~ok;
        tt(ii).pierce.depth(nok)=[];
        tt(ii).pierce.distance(nok)=[];
        tt(ii).pierce.time(nok)=[];
        if(~generic)
            tt(ii).pierce.latitude(nok)=[];
            tt(ii).pierce.longitude(nok)=[];
        end
    end
end

% formatted listing
if(nargout==0)
    % header
    disp(' ')
    disp(['Model: ' model])
    disp('Distance   Depth   Phase        Travel    Ray Param   Purist    Purist')
    disp('  (deg)     (km)   Name         Time (s)  p (s/deg)  Distance   Name  ')
    disp('--------------------------------------------------------------------------')
    
    % list phase info
    for ii=1:narr
        disp(' ')
        fprintf(' %7.2f  %6.1f   %-10s   %7.2f   %7.3f    %7.2f  = %-10s\n',...
            tt(ii).mindistance,tt(ii).depth,...
            tt(ii).phase,tt(ii).time,...
            tt(ii).rayparameter,tt(ii).distance,...
            tt(ii).puristphase);
        fprintf('%16.1f %22.2f %20.2f\n',[tt(ii).pierce.depth ...
            tt(ii).pierce.time tt(ii).pierce.distance]');
    end
else
    varargout{1}=tt;
end

end
