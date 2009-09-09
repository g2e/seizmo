function tt=tauptime(varargin)
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
%    Description: TAUPTIME(...) (no outputs, w/ or w/o inputs) displays a
%     formatted list of information on many seismic phases found at a
%     random distance (0 to 180 deg) from an event with a depth of 0km.  
%     The 1D Earth model utilized by default is IASP91 (see the MODEL
%     option to adjust this).  The default phase list is 'ttbasic'
%     (equivalent to setting phases to BASIC in Brian Kennett's TTIMES
%     program) which lists many common phases by default.  See option
%     PHASES to adjust the phase list.  A variety of options (DEG, KM, STA
%     & EVT) allow for changing the distance from the event.
%
%     TT=TAUPTIME(...) (1 output, w/ or w/o inputs) returns a structure
%     array with the following fields:
%             TT(index).phase        - seismic phase name
%                      .distance     - purist distance (deg)
%                      .depth        - depth of earthquake (km)
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
%     'ak135'.  See the TauP program/documentation for more.  The default
%     model is 'iasp91'.
%
%     TT=TAUPTIME(...,'H|Z|DEP|EVDP|DEPTH',DEPTH,...) sets the event depth
%     to DEPTH.  DEPTH should be in kilometers!  The default is 0.
%
%     TT=TAUPTIME(...,'P|PH|PHASES',PHASES,...) sets the phase list to
%     PHASES.  PHASES must be a comma separated list of seismic phases.  
%     See the TauP documentation for conventions and valid phase names.
%     Multiple calls are honored.  The default phase list is 'ttbasic'.
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
%     Return info on several P arrivals expected at 25 degrees from an
%     event:
%      tauptime('p','ttp','d',25)
%
%     Some other valid examples:
%      tauptime('dep',50,'ph','P,S','deg',45.6)
%      tauptime('mod','prem','dep',50,'ph','Pdiff,PKP',...
%          'sta',[-40 -100],'evt',[30,50])
%
%    See also: taup, tauppath, taupcurve, tauppierce

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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep.  5, 2009 at 21:00 GMT

% todo:

% check nargin
if(mod(nargin,2))
    error('matTaup:tauptime:badNumOptions','Unpaired option(s)!');
end

% initialize java code
import edu.sc.seis.TauP.*;
import java.io.*;
import java.lang.*;
import java.util.*;
import java.util.zip.*;

% default options
model='iasp91';
depth='0';
phases='ttbasic';

% check options
pargs=cell(0); np=0;
dargs=cell(0); nd=0;
d=false; k=d; s=d; e=d; a=d; b=d;
for i=1:2:nargin
    switch lower(varargin{i})
        case {'m' 'mod' 'model'}
            if(isempty(varargin{i+1})); continue; end
            if(~ischar(varargin{i+1}))
                error('matTaup:tauptime:badInput',...
                    'MODEL must be a string (like ''prem'')!');
            end
            model=varargin{i+1}(:)';
        case {'h' 'z' 'dep' 'evdp' 'depth'}
            if(isempty(varargin{i+1})); continue; end
            if(~isscalar(varargin{i+1}) ...
                    || ~isreal(varargin{i+1}) || varargin{i+1}<0)
                error('matTaup:tauptime:badInput',...
                    'DEPTH must be a positive number (in km)!');
            end
            depth=num2str(varargin{i+1});
        case {'p' 'ph' 'phases'}
            if(isempty(varargin{i+1})); continue; end
            if(~ischar(varargin{i+1}))
                error('matTaup:tauptime:badInput',...
                    'PHASES must be a string (like ''P,S'')!');
            end
            np=np+2;
            pargs(1,np-1:np)={'-ph' varargin{i+1}};
        case {'d' 'deg' 'gcarc' 'degrees'}
            if(isempty(varargin{i+1})); continue; end
            if(~isscalar(varargin{i+1}) ...
                    || ~isreal(varargin{i+1}) || varargin{i+1}<0)
                error('matTaup:tauptime:badInput',...
                    'DEG must be a positive number (in degrees)!');
            end
            nd=nd+2;
            dargs(1,nd-1:nd)={'-deg' num2str(varargin{i+1})};
            d=true;
        case {'k' 'km' 'dist' 'kilometers'}
            if(isempty(varargin{i+1})); continue; end
            if(~isscalar(varargin{i+1}) ...
                    || ~isreal(varargin{i+1}) || varargin{i+1}<0)
                error('matTaup:tauptime:badInput',...
                    'KM must be a positive number (in kilometers)!');
            end
            nd=nd+2;
            dargs(1,nd-1:nd)={'-km' num2str(varargin{i+1})};
            k=true;
        case {'s' 'st' 'sta' 'station'}
            if(isempty(varargin{i+1})); continue; end
            if(numel(varargin{i+1})~=2 ...
                    || ~all(isreal(varargin{i+1})))
                error('matTaup:tauptime:badInput',...
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
                error('matTaup:tauptime:badInput',...
                    'EVT must be 2 real numbers ([LAT LON] in degrees)!');
            end
            nd=nd+3;
            dargs(1,nd-2:nd)={'-evt' num2str(varargin{i+1}(1)) ...
                num2str(varargin{i+1}(2))};
            e=true;
        case {'a' 'az' 'azi' 'azimuth'}
            if(isempty(varargin{i+1})); continue; end
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}))
                error('matTaup:tauptime:badInput',...
                    'AZ must be a scalar number (in degrees)!');
            end
            nd=nd+2;
            dargs(1,nd-1:nd)={'-az' num2str(varargin{i+1})};
            a=true;
        case {'b' 'baz' 'bazi' 'backazimuth'}
            if(isempty(varargin{i+1})); continue; end
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}))
                error('matTaup:tauptime:badInput',...
                    'BAZ must be a scalar number (in degrees)!');
            end
            nd=nd+2;
            dargs(1,nd-1:nd)={'-baz' num2str(varargin{i+1})};
            b=true;
        otherwise
            error('matTaup:tauptime:badOption',...
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
if(d || k || (s && e) || (s && b && (d || k)) || (e && a && (d || k)))
    inArgs=[inArgs dargs];
else
    % default to random distance
    disp('Distance not given or could not be determined - using random!')
    inArgs=[inArgs dargs {'-deg' num2str(rand*180)}];
end

% debug
%disp(inArgs);
%arrivals=MatTauP_Time.run_time(inArgs);

% attempt run
try
    arrivals=MatTauP_Time.run_time(inArgs);
catch
    % oops!
    error('matTaup:tauptime:runFailed',...
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
        fprintf(' %7.2f  %6.1f   %-10s   %7.2f   %7.3f    %7.2f  = %-10s\n',...
            abssawmod(arrivals(ii).getDistDeg,180),...
            arrivals(ii).getSourceDepth,...
            char(arrivals(ii).getName),arrivals(ii).getTime,...
            arrivals(ii).getRayParam/R2D,arrivals(ii).getDistDeg,...
            char(arrivals(ii).getPuristName));
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
end

end

function [c]=abssawmod(a,b)
% returns distance range always as the minor arc length (0-180deg)

n=round(0.5.*a./b);
s=1-2*mod(n,2);
c=s.*(a-2.*n.*b);
if(isscalar(b))
    if(b==0)
        c=a;
    end
else
    d=b==0;
    c(d)=a(d);
end

c=abs(c);

end
