function [varargout]=tauptime(varargin)
%TAUPTIME    Calculate travel times using the TauP toolkit
%
%    Usage:    tauptime(...)
%           tt=tauptime(...)
%           tt=tauptime(...,'m|mod|model',model,...)
%           tt=tauptime(...,'h|z|dep|evdp|depth',depth,...)
%           tt=tauptime(...,'p|ph|phases',phases,...)
%           tt=tauptime(...,'d|deg|gcarc|degrees',degdist,...)
%           tt=tauptime(...,'k|km|dist|kilometers',kmdist,...)
%           tt=tauptime(...,'s|st|sta|station',[lat lon],...)
%           tt=tauptime(...,'e|ev|evt|event',[lat lon],...)
%
%    Description:
%     TAUPTIME(...) (no outputs, w/ or w/o inputs) displays a formatted
%     list of information on many seismic phases found at a random distance
%     (0 to 180 deg) from an event with a depth of 0km.  The 1D Earth model
%     utilized by default is IASP91 (see the MODEL option to adjust this).
%     The default phase list is 'ttbasic' (equivalent to setting phases to
%     BASIC in Brian Kennett's TTIMES program) which lists many common
%     phases by default.  See option PHASES to adjust the phase list.  A
%     variety of options (DEG, KM, STA & EVT) allow for changing the
%     distance from the event.
%
%     TT=TAUPTIME(...) (1 output, w/ or w/o inputs) returns a structure
%     array with the following fields:
%             TT(index).modelname    - name of velocity model
%                      .event        - event location if known ([lat lon])
%                      .station      - stn location if known ([lat lon])
%                      .depth        - depth of earthquake (km)
%                      .distance     - true distance (deg)
%                      .mindistance  - least distance (deg)
%                      .phase        - seismic phase name
%                      .puristphase  - verbose seismic phase name
%                      .time         - travel time (sec)
%                      .rayparameter - ray parameter (sec/deg)
%     Each phase has its own indice in the struct array TT.  Use TT(index)
%     to access individual phase information.
%
%     *********************************************************
%     All the following field/value pair options may be entered
%     at the command line in any order.
%     *********************************************************
%
%     TT=TAUPTIME(...,'M|MOD|MODEL',MODEL,...) sets the 1D Earth model to
%     MODEL.  Accepts a variety of common models like 'prem', 'iasp91',
%     'ak135'.  See the TauP program/documentation for more.  You may also
%     give a 1DMODEL struct like from CMB_1DMODEL_LIBRARY.  The default
%     model is 'iasp91'.
%
%     TT=TAUPTIME(...,'H|Z|DEP|EVDP|DEPTH',DEPTH,...) sets the event depth
%     to DEPTH.  DEPTH should be in kilometers!  The default is 0.
%
%     TT=TAUPTIME(...,'P|PH|PHASES',PHASES,...) sets the phase list to
%     PHASES.  PHASES must be a comma separated list of seismic phases.  
%     See the TauP documentation for conventions and valid phase names.
%     The default phase list is 'ttbasic'.
%
%     TT=TAUPTIME(...,'D|DEG|GCARC|DEGREES',DEGDIST,...) sets the event
%     distance to DEGDIST.  DEGDIST is expected to be in angular degrees.
%     The default value is random (somewhere between 0 and 180 degrees).
%     DEGDIST must be scalar.
%
%     TT=TAUPTIME(...,'K|KM|DIST|KILOMETERS',KMDIST,...) sets the event
%     distance to KMDIST.  KMDIST is expected to be in kilometers.  There
%     is no default value for this option.  KMDIST must be scalar.
%
%     TT=TAUPTIME(...,'S|ST|STA|STATION',[LAT LON],...) sets the seismic
%     station location to [LAT LON].  LAT and LON must be in degrees.  This
%     option must be paired with the EVT option to calculate the event
%     distance.  There is no default value for this option.
%
%     TT=TAUPTIME(...,'E|EV|EVT|EVENT',[LAT LON],...) sets the earthquake
%     location to [LAT LON].  LAT and LON must be in degrees.  This option
%     must be paired with the STA option to calculate the event distance.
%     There is no default value for this option.
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
%     % Return info on several P arrivals
%     % expected at 25 degrees from an event:
%     tauptime('p','ttp','d',25)
%
%     % Some other valid examples:
%     tauptime('dep',50,'ph','P,S','deg',45.6)
%     tauptime('mod','prem','dep',50,'ph','Pdiff,PKP',...
%         'sta',[-40 -100],'evt',[30,50])
%
%    See also: TAUPPATH, TAUPCURVE, TAUPPIERCE, TAUPCREATE, TAUP

%     Version History:
%        June 29, 2009 - added defaults for depth, phase, model
%                        also did some minor code cleaning
%        Aug. 27, 2009 - major script update, name change to avoid breakage
%                        due to input/output changes
%        Aug. 30, 2009 - preserve ordering of distance options, add az &
%                        baz (not useful here, but good for others), doc
%                        cleaning
%        Sep.  2, 2009 - allow multiple calls to PHASES
%        Sep.  5, 2009 - minor doc update
%        Sep. 30, 2009 - changed abssawmod to abslatmod
%        Nov. 13, 2009 - dropped abslatmod for getModuloDistDeg, dropped
%                        some import calls
%        Jan.  6, 2011 - add matTaup.jar to dynamic java classpath if
%                        necessary
%        Feb. 16, 2012 - works with Octave+Octave-Java
%        Feb. 24, 2012 - switch to native taup
%        Jan. 26, 2014 - no longer need to update jar filenames
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 26, 2014 at 17:15 GMT

% todo:

% check nargin
if(mod(nargin,2))
    error('matTaup:tauptime:badNumOptions','Unpaired option(s)!');
end

% try adding *.jar if no TauP class exists
if(~exist('edu.sc.seis.TauP.TauP_Time','class'))
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

% check options
if(~iscellstr(varargin(1:2:end)))
    error('TauP:tauptime:badInput',...
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
                error('TauP:tauptime:badInput',...
                    'MODEL must be a string (like ''prem'')!');
            end
        case {'h' 'z' 'dep' 'evdp' 'depth'}
            if(~isscalar(varargin{i+1}) ...
                    || ~isreal(varargin{i+1}) || varargin{i+1}<0)
                error('TauP:tauptime:badInput',...
                    'DEPTH must be a positive number (in km)!');
            end
            depth=varargin{i+1};
        case {'p' 'ph' 'phases'}
            if(~ischar(varargin{i+1}))
                error('TauP:tauptime:badInput',...
                    'PHASES must be a string (like ''P,S'')!');
            end
            phases=varargin{i+1};
            phases=strtrim(getwords(phases,','));
        case {'d' 'deg' 'gcarc' 'degrees'}
            if(~isscalar(varargin{i+1}) ...
                    || ~isreal(varargin{i+1}) || varargin{i+1}<0)
                error('TauP:tauptime:badInput',...
                    'DEG must be a positive number (in degrees)!');
            end
            deg=varargin{i+1};
        case {'k' 'km' 'dist' 'kilometers'}
            if(~isscalar(varargin{i+1}) ...
                    || ~isreal(varargin{i+1}) || varargin{i+1}<0)
                error('TauP:tauptime:badInput',...
                    'KM must be a positive number (in kilometers)!');
            end
            km=varargin{i+1};
        case {'s' 'st' 'sta' 'station'}
            if(~isequal([1 2],size(varargin{i+1})) ...
                    || ~all(isreal(varargin{i+1})))
                error('TauP:tauptime:badInput',...
                    'STA must be 2 real numbers ([LAT LON] in degrees)!');
            end
            st=varargin{i+1};
        case {'e' 'ev' 'evt' 'event'}
            if(~isequal([1 2],size(varargin{i+1})) ...
                    || ~all(isreal(varargin{i+1})))
                error('TauP:tauptime:badInput',...
                    'EVT must be 2 real numbers ([LAT LON] in degrees)!');
            end
            ev=varargin{i+1};
        case {'a' 'az' 'azi' 'azimuth'}
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}))
                error('TauP:tauptime:badInput',...
                    'AZ must be a scalar number (in degrees)!');
            end
            az=varargin{i+1};
        case {'b' 'baz' 'bazi' 'backazimuth'}
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}))
                error('TauP:tauptime:badInput',...
                    'BAZ must be a scalar number (in degrees)!');
            end
            baz=varargin{i+1};
        otherwise
            error('TauP:tauptime:badOption',...
                'Unknown Option: %s',varargin{i});
    end
end

% check geometry
sc=javaObject('edu.sc.seis.TauP.SphericalCoords');
d=~isempty(deg); k=~isempty(km); s=~isempty(st);
e=~isempty(ev); a=~isempty(az); b=~isempty(baz);
if(d && k)
    % over specified
    error('TauP:tauptime:badInput',...
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
    %generic=false;
    if(e && s)
        deg=sc.distance(ev(1),ev(2),st(1),st(2));
        %az=sc.azimuth(ev(1),ev(2),st(1),st(2));
        %baz=sc.azimuth(st(1),st(2),ev(1),ev(2));
    elseif(s)
        ev(1)=sc.latFor(st(1),st(2),deg,baz);
        ev(2)=sc.lonFor(st(1),st(2),deg,baz);
        %az=sc.azimuth(ev(1),ev(2),st(1),st(2));
    elseif(e)
        st(1)=sc.latFor(ev(1),ev(2),deg,az);
        st(2)=sc.lonFor(ev(1),ev(2),deg,az);
        %baz=sc.azimuth(st(1),st(2),ev(1),ev(2));
    end
elseif((d || k) && ~((d && k) || s || e || a || b))
    % generic geometry
    %generic=true;
else
    error('TauP:tauptime:badInput',...
        'Event-Station geometry under/over specified!');
end

% debug
%[generic ev st deg az baz]

% create time object for velocity model
ttobj=javaObject('edu.sc.seis.TauP.TauP_Time',model);

% calculate times for the specified depth, phase & distance
ttobj.setSourceDepth(depth);
ttobj.setPhaseNames(phases);
ttobj.calculate(deg);
narr=ttobj.getNumArrivals;

% conversion
R2D=180/pi;

% struct output
tt(1:narr,1)=struct('modelname',[],'event',[],'station',[],'depth',[],...
    'distance',[],'mindistance',[],'phase',[],'puristphase',[],...
    'time',[],'rayparameter',[]);
for ii=1:narr
    % get phase info
    arr=ttobj.getArrival(ii-1);
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
end

% formatted listing
if(nargout==0)
    % header
    disp(' ')
    disp(['Model: ' modelname])
    disp('Distance   Depth   Phase        Travel    Ray Param   Purist    Purist')
    disp('  (deg)     (km)   Name         Time (s)  p (s/deg)  Distance   Name  ')
    disp('--------------------------------------------------------------------------')
    
    % loop over arrivals
    for ii=1:narr
        % list phase info
        fprintf(' %7.2f  %6.1f   %-10s   %7.2f   %7.3f    %7.2f  = %-10s\n',...
            tt(ii).mindistance,tt(ii).depth,...
            tt(ii).phase,tt(ii).time,...
            tt(ii).rayparameter,tt(ii).distance,...
            tt(ii).puristphase);
    end
else
    varargout{1}=tt;
end

end
