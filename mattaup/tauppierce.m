function tt=tauppierce(varargin)
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
%    Description: TAUPPIERCE(...) (no outputs, w/ or w/o inputs) displays a
%     formatted list of information on many seismic phases found at a
%     random distance from an event with a depth of 0km.  This is the same
%     as for TAUPTIME (see the second form for the difference between
%     TAUPTIME and TAUPPIERCE).  The 1D Earth model utilized by default is
%     IASP91 (see the MODEL option to adjust this).  The default phase list
%     is 'ttbasic' (equivalent to setting phases to BASIC in Brian
%     Kennett's TTIMES program) which lists many common phases by default.
%     See option PHASES to adjust the phase list.  A variety of options
%     (DEG, KM, STA & EVT) allow for changing the distance from the event.
%
%     TT=TAUPPIERCE(...) (1 output, w/ or w/o inputs) returns a structure
%     array with the following layout:
%       TT(index).phase               - seismic phase name
%                .distance            - purist distance (deg)
%                .depth               - depth of earthquake (km)
%                .time                - travel time (sec)
%                .rayparameter        - slowness (sec/deg)
%                .pierce.rayparameter - slowness values at pierce points
%                .pierce.time         - travel times (sec) at pierce points
%                .pierce.distance     - distance (deg) of pierce points
%                .pierce.depth        - depth (km) of pierce points
%                .pierce.latitude     - latitude (deg) of pierce points
%                .pierce.longitude    - longitude (deg) of pierce points
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
%     'ak135'.  See the TauP program/documentation for more.  The default
%     model is 'iasp91'.
%
%     TT=TAUPPIERCE(...,'H|Z|DEP|EVDP|DEPTH',DEPTH,...) sets the event
%     depth to DEPTH.  DEPTH should be in kilometers!  The default is 0.
%
%     TT=TAUPPIERCE(...,'P|PH|PHASES',PHASES,...) sets the phase list to
%     PHASES.  PHASES must be a comma separated list of seismic phases.  
%     See the TauP documentation for conventions and valid phase names.
%     Multiple calls are honored.  The default phase list is 'ttbasic'.
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
%     TRUE will just return turning (aka bottoming) points.  Info for other
%     pierce points is set to zero (blame TauP?).  This option may be
%     combined with UNDER, REV, and NODISCON (ie they do not exclude one
%     another - although REV is a super-set of TURN).  By default TURN is
%     set to FALSE.
%
%     TT=TAUPPIERCE(...,'U|UNDER|UNDERSIDE',TRUE|FALSE,...) toggles whether
%     or not to just return underside reflection points.  This also returns
%     the start and end points of the ray paths (blame TauP - the REV
%     option returns underside points with the start/end points excluded
%     but also returns turning points).  Info for other pierce points is
%     set to zero (blame TauP?).  This option may be combined with the
%     TURN, REV, and NODISCON options (no mutual exclusion).  By default
%     UNDER is set to FALSE.
%
%     TT=TAUPPIERCE(...,'R|REV|REVERSE',TRUE|FALSE,...) toggles the
%     returning of only points where a reversal occurs (turning and
%     reflection points).  Note that setting REV to TRUE is not the same as
%     setting both the TURN and UNDER options to TRUE, as the UNDER option
%     returns the start and end points of the ray path (blame TauP) and REV
%     does not.  Info for other pierce points is set to zero (blame TauP?).
%     By default REV is set to FALSE.
%
%     TT=TAUPPIERCE(...,'N|ND|NODISCON',TRUE|FALSE,...) toggles returning
%     just the points in the PIERCE option.  Info for other pierce points
%     is set to zero (blame TauP?).  May be combined with the TURN, UNDER,
%     and REV options without excluding each others results.  By default
%     NODISCON is set to FALSE.
%
%     TT=TAUPPIERCE(...,'X|PNTS|PIERCE',DEPTHS,...) adds the depths in
%     DEPTHS to the ray path pierce points to calculate.  Honors all calls
%     when called multiple times.  Depths must be a numeric array. There is
%     no default value for this option.
%
%    Notes:
%     - These scripts require the included file:
%          mattaup/lib/matTaup.jar
%       to be added in Matlab's javaclasspath.  You may use the functions
%       javaaddpath and javarmpath to alter the dynamic portion of the path
%       or you will need to add the jar file to the Matlab system file
%       'classpath.txt'.  Use the command 'edit classpath.txt' in Matlab to
%       add the the jar file (be careful and use the full path!).  This may
%       require administrator privileges.
%
%     - This script passes all position/distance options in the order they
%       are given.  This leaves the interpretation of conflicting inputs to
%       TauP.
%
%     - MatTauP is only a wrapping program for TauP toolkit, which is
%       developed by:
%        H. Philip Crotwell, Thomas J. Owens, Jeroen Ritsema
%        Department of Geological Sciences
%        University of South Carolina
%        http://www.seis.sc.edu
%        crotwell@seis.sc.edu
%
%     - MatTauP was written by:
%        Qin Li 
%        Unverisity of Washington
%        qinli@u.washington.edu
%        Nov, 2002
%
%    Examples:
%     P'P' pierce points for a 44km deep event at some random distance:
%      tauppierce('p','PKIKPPKIKP','z',44)
%
%     Add some other interesting piercing depths:
%      tauppierce('p','PKIKPPKIKP','z',44,'x',[1800 2580])
%
%     Return underside reflection points for pS^640S (note 3 nonzeros):
%      tauppierce('p','pS^640S','z',33,'d',120,'under',true)
%
%    See also: TAUP, TAUPPATH, TAUPCURVE, TAUPTIME

%     Version History:
%        Sep.  2, 2009 - major revision of script, name change to avoid
%                        breakage due to input/output changes
%        Sep.  5, 2009 - minor doc update
%        Sep. 30, 2009 - changed abssawmod to abslatmod
%        Nov. 13, 2009 - dropped abslatmod for getModuloDistDeg, dropped
%                        some import calls
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 13, 2009 at 17:25 GMT

% todo:

% check nargin
if(mod(nargin,2))
    error('matTaup:tauppierce:badNumOptions','Unpaired option(s)!');
end

% initialize java code
import edu.sc.seis.TauP.*;

% default options
model='iasp91';
depth='0';
phases='ttbasic';

% check options
pargs=cell(0); np=0;
xargs=cell(0); nx=0;
iargs=cell(0); ni=0;
dargs=cell(0); nd=0;
d=false; k=d; s=d; e=d; a=d; b=d;
for i=1:2:nargin
    switch lower(varargin{i})
        case {'m' 'mod' 'model'}
            if(isempty(varargin{i+1})); continue; end
            if(~ischar(varargin{i+1}))
                error('matTaup:tauppierce:badInput',...
                    'MODEL must be a string (like ''prem'')!');
            end
            model=varargin{i+1}(:)';
        case {'h' 'z' 'dep' 'evdp' 'depth'}
            if(isempty(varargin{i+1})); continue; end
            if(~isscalar(varargin{i+1}) ...
                    || ~isreal(varargin{i+1}) || varargin{i+1}<0)
                error('matTaup:tauppierce:badInput',...
                    'DEPTH must be a positive number (in km)!');
            end
            depth=num2str(varargin{i+1});
        case {'p' 'ph' 'phases'}
            if(isempty(varargin{i+1})); continue; end
            if(~ischar(varargin{i+1}))
                error('matTaup:tauppierce:badInput',...
                    'PHASES must be a string (like ''P,S'')!');
            end
            np=np+2;
            pargs(1,np-1:np)={'-ph' varargin{i+1}};
        case {'d' 'deg' 'gcarc' 'degrees'}
            if(isempty(varargin{i+1})); continue; end
            if(~isscalar(varargin{i+1}) ...
                    || ~isreal(varargin{i+1}) || varargin{i+1}<0)
                error('matTaup:tauppierce:badInput',...
                    'DEG must be a positive number (in degrees)!');
            end
            nd=nd+2;
            dargs(1,nd-1:nd)={'-deg' num2str(varargin{i+1})};
            d=true;
        case {'k' 'km' 'dist' 'kilometers'}
            if(isempty(varargin{i+1})); continue; end
            if(~isscalar(varargin{i+1}) ...
                    || ~isreal(varargin{i+1}) || varargin{i+1}<0)
                error('matTaup:tauppierce:badInput',...
                    'KM must be a positive number (in kilometers)!');
            end
            nd=nd+2;
            dargs(1,nd-1:nd)={'-km' num2str(varargin{i+1})};
            k=true;
        case {'s' 'st' 'sta' 'station'}
            if(isempty(varargin{i+1})); continue; end
            if(numel(varargin{i+1})~=2 ...
                    || ~all(isreal(varargin{i+1})))
                error('matTaup:tauppierce:badInput',...
                    'STA must be 2 real numbers ([LAT LON] in degrees)!');
            end
            nd=nd+3;
            dargs(1,nd-2:nd)={'-sta' num2str(varargin{i+1}(1)) ...
                num2str(varargin{i+1}(2))};
            s=true;
        case {'e' 'ev' 'evt' 'event'}
            if(isempty(varargin{i+1})); continue; end
            if(numel(varargin{i+1})~=2 ...
                    || ~all(isreal(varargin{i+1})))
                error('matTaup:tauppierce:badInput',...
                    'EVT must be 2 real numbers ([LAT LON] in degrees)!');
            end
            nd=nd+3;
            dargs(1,nd-2:nd)={'-evt' num2str(varargin{i+1}(1)) ...
                num2str(varargin{i+1}(2))};
            e=true;
        case {'a' 'az' 'azi' 'azimuth'}
            if(isempty(varargin{i+1})); continue; end
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}))
                error('matTaup:tauppierce:badInput',...
                    'AZ must be a scalar number (in degrees)!');
            end
            nd=nd+2;
            dargs(1,nd-1:nd)={'-az' num2str(varargin{i+1})};
            a=true;
        case {'b' 'baz' 'bazi' 'backazimuth'}
            if(isempty(varargin{i+1})); continue; end
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}))
                error('matTaup:tauppierce:badInput',...
                    'BAZ must be a scalar number (in degrees)!');
            end
            nd=nd+2;
            dargs(1,nd-1:nd)={'-baz' num2str(varargin{i+1})};
            b=true;
        case {'t' 'turn' 'turning'}
            if(isempty(varargin{i+1})); continue; end
            if(~isscalar(varargin{i+1}) || ~islogical(varargin{i+1}))
                error('matTaup:tauppierce:badInput',...
                    'TURN must be a scalar logical!');
            end
            if(varargin{i+1})
                ni=ni+1;
                iargs(1,ni)={'-turn'};
            end
        case {'u' 'under' 'underside'}
            if(isempty(varargin{i+1})); continue; end
            if(~isscalar(varargin{i+1}) || ~islogical(varargin{i+1}))
                error('matTaup:tauppierce:badInput',...
                    'UNDER must be a scalar logical!');
            end
            if(varargin{i+1})
                ni=ni+1;
                iargs(1,ni)={'-under'};
            end
        case {'r' 'rev' 'reverse'}
            if(isempty(varargin{i+1})); continue; end
            if(~isscalar(varargin{i+1}) || ~islogical(varargin{i+1}))
                error('matTaup:tauppierce:badInput',...
                    'REV must be a scalar logical!');
            end
            if(varargin{i+1})
                ni=ni+1;
                iargs(1,ni)={'-rev'};
            end
        case {'n' 'nd' 'nodiscon'}
            if(isempty(varargin{i+1})); continue; end
            if(~isscalar(varargin{i+1}) || ~islogical(varargin{i+1}))
                error('matTaup:tauppierce:badInput',...
                    'NODISCON must be a scalar logical!');
            end
            if(varargin{i+1})
                ni=ni+1;
                iargs(1,ni)={'-nodiscon'};
            end
        case {'x' 'pnts' 'pierce'}
            if(isempty(varargin{i+1})); continue; end
            if(any(~isreal(varargin{i+1}(:))) || any(varargin{i+1}(:)<0))
                error('matTaup:tauppierce:badInput',...
                    'PIERCE must be positive numbers (in kilometers!');
            end
            nx=nx+2;
            tmp=strcat(num2str(varargin{i+1}(:)),',')'; % comma delimited
            tmp=tmp(:)';      % assure row vector
            tmp=tmp(tmp~=32); % remove spaces
            tmp=tmp(1:end-1); % trim off trailing comma
            xargs(1,nx-1:nx)={'-pierce' tmp};
        otherwise
            error('matTaup:tauppierce:badOption',...
                'Unknown Option: %s',varargin{i});
    end
end

% set up inputs
inArgs{1}='-mod';
inArgs{2}=model;
inArgs{3}='-h';
inArgs{4}=depth;
if(np>0)
    inArgs=[inArgs pargs];
else
    inArgs{5}='-ph';
    inArgs{6}=phases;
end
if(ni>0); inArgs=[inArgs iargs]; end
if(nx>0); inArgs=[inArgs xargs]; end
if(d || k || (s && e) || (s && b && (d || k)) || (e && a && (d || k)))
    inArgs=[inArgs dargs];
else
    % default to random distance
    disp('Distance not given or could not be determined - using random!')
    inArgs=[inArgs dargs {'-deg' num2str(rand*180)}];
end

% debug
%disp(inArgs);
%arrivals=MatTauP_Pierce.run_pierce(inArgs);

% attempt run
try
    arrivals=MatTauP_Pierce.run_pierce(inArgs);
catch
    % oops!
    error('matTaup:tauppierce:runFailed',...
        ['Java exception occurred! Please check input options and\n'...
        'make sure your classpath.txt file includes matTaup.jar!']);
end

% conversion
R2D=180/pi;

% formatted listing
if(nargout==0)
    % header
    disp(' ')
    disp(['Model: ' model])
    disp('Distance   Depth   Phase        Travel    Ray Param   Purist    Purist')
    disp('  (deg)     (km)   Name         Time (s)  p (s/deg)  Distance   Name  ')
    disp('--------------------------------------------------------------------------')
    
    % list phase info
    for ii=1:arrivals.length
        disp(' ')
        fprintf(' %7.2f  %6.1f   %-10s   %7.2f   %7.3f    %7.2f  = %-10s\n',...
            arrivals(ii).getModuloDistDeg,arrivals(ii).getSourceDepth,...
            char(arrivals(ii).getName),arrivals(ii).getTime,...
            arrivals(ii).getRayParam/R2D,arrivals(ii).getDistDeg,...
            char(arrivals(ii).getPuristName));
        disp(char(strcat(...
            {'  '},...
            num2str(arrivals(ii).getMatPath.dist,'%7.2f'),...
            {'  '},num2str(arrivals(ii).getMatPath.depth,'%6.1f'),...
            {'                '},...
            num2str(arrivals(ii).getMatPath.time,'%7.2f'),...
            {'   '},num2str(arrivals(ii).getMatPath.p/R2D,'%7.3f'),...
            {'    '},...
            num2str(arrivals(ii).getMatPath.dist,'%7.2f'))))
    end
    
    return
end

% struct output
tt(1:arrivals.length)=struct('phase',[],'distance',[],'depth',[],...
    'time',[],'rayparameter',[]);
for ii=1:arrivals.length
    tt(ii).time=arrivals(ii).getTime;
    tt(ii).distance=arrivals(ii).getDistDeg;
    tt(ii).depth=arrivals(ii).getSourceDepth;
    tt(ii).phase=char(arrivals(ii).getName);
    tt(ii).rayparameter=arrivals(ii).getRayParam/R2D;
    tt(ii).pierce.rayparameter=arrivals(ii).getMatPath.p/R2D;
    tt(ii).pierce.time=arrivals(ii).getMatPath.time;
    tt(ii).pierce.distance=arrivals(ii).getMatPath.dist;
    tt(ii).pierce.depth=arrivals(ii).getMatPath.depth;
    tt(ii).pierce.latitude=arrivals(ii).getMatPath.lat;
    tt(ii).pierce.longitude=arrivals(ii).getMatPath.lon;
end

end
