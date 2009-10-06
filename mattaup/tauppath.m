function [tt]=tauppath(varargin)
%TAUPPATH    Calculate ray paths using the TauP toolkit
%
%    Usage:    tauppath(...)
%           tt=tauppath(...)
%           tt=tauppath(...,'m|mod|model',model,...)
%           tt=tauppath(...,'h|z|dep|evdp|depth',depth,...)
%           tt=tauppath(...,'p|ph|phases',phases,...)
%           tt=tauppath(...,'d|deg|gcarc|degrees',degdist,...)
%           tt=tauppath(...,'k|km|dist|kilometers',kmdist,...)
%           tt=tauppath(...,'s|st|sta|station',[lat lon],...)
%           tt=tauppath(...,'e|ev|evt|event',[lat lon],...)
%           tt=tauppath(...,'a|az|azi|azimuth',azi,...)
%           tt=tauppath(...,'b|baz|bazi|backazimuth',bazi,...)
%
%    Description: TAUPPATH(...) (no outputs, w/ or w/o inputs) displays a
%     formatted list of information on many seismic phases found at a
%     random distance from an event with a depth of 0km.  This is the same
%     as for TAUPTIME, but TAUPPATH also plots the ray paths for the listed
%     phases.  The 1D Earth model utilized by default is IASP91 (see the
%     MODEL option to adjust this).  The default phase list is 'ttbasic'
%     (equivalent to setting phases to BASIC in Brian Kennett's TTIMES
%     program) which lists many common phases by default.  See option
%     PHASES to adjust the phase list.  A variety of options (DEG, KM, STA
%     & EVT) allow for changing the distance from the event.
%
%     TT=TAUPPATH(...) (1 output, w/ or w/o inputs) returns a structure
%     array with the following layout:
%        TT(index).phase             - seismic phase name
%                 .distance          - purist distance (deg)
%                 .depth             - depth of earthquake (km)
%                 .time              - travel time (sec)
%                 .rayparameter      - ray parameter (sec/deg)
%                 .path.rayparameter - rayparameter values for path points
%                 .path.time         - travel times (sec) at path points
%                 .path.distance     - distance (deg) of path points
%                 .path.depth        - depth (km) of path points
%                 .path.latitude     - latitude (deg) of path points
%                 .path.longitude    - longitude (deg) of path points
%     Each phase has its own indice in the struct array TT.  Use TT(index)
%     to access individual phase information.  The sub-structure 'path'
%     contains several fields populated with arrays of info on the ray path
%     for the phase.
%
%     *********************************************************
%     All the following field/value pair options may be entered
%     at the command line in any order.
%     *********************************************************
%
%     TT=TAUPPATH(...,'M|MOD|MODEL',MODEL,...) sets the 1D Earth model to
%     MODEL.  Accepts a variety of common models like 'prem', 'iasp91',
%     'ak135'.  See the TauP program/documentation for more.  The default
%     model is 'iasp91'.
%
%     TT=TAUPPATH(...,'H|Z|DEP|EVDP|DEPTH',DEPTH,...) sets the event depth
%     to DEPTH.  DEPTH should be in kilometers!  The default is 0.
%
%     TT=TAUPPATH(...,'P|PH|PHASES',PHASES,...) sets the phase list to
%     PHASES.  PHASES must be a comma separated list of seismic phases.  
%     See the TauP documentation for conventions and valid phase names.
%     Multiple calls are honored.  The default phase list is 'ttbasic'.
%
%     TT=TAUPPATH(...,'D|DEG|GCARC|DEGREES',DEGDIST,...) sets the event
%     distance to DEGDIST.  DEGDIST is expected to be in angular degrees.
%     The default value is random (somewhere between 0 and 180 degrees).
%     DEGDIST must be scalar.
%
%     TT=TAUPPATH(...,'K|KM|DIST|KILOMETERS',KMDIST,...) sets the event
%     distance to KMDIST.  KMDIST is expected to be in kilometers.  There
%     is no default value for this option.  KMDIST must be scalar.
%
%     TT=TAUPPATH(...,'S|ST|STA|STATION',[LAT LON],...) sets the seismic
%     station location to [LAT LON].  LAT and LON must be in degrees.  This
%     option must be paired with the EVT option to calculate the event
%     distance.  There is no default value for this option.
%
%     TT=TAUPPATH(...,'E|EV|EVT|EVENT',[LAT LON],...) sets the earthquake
%     location to [LAT LON].  LAT and LON must be in degrees.  This option
%     must be paired with the STA option to calculate the event distance.
%     There is no default value for this option.
%
%     TT=TAUPPATH(...,'A|AZ|AZI|AZIMUTH',AZI,...) sets the azimuth from the
%     earthquake to the seismic station.  Note that this option really
%     isn't useful unless paired with the EVT option and one of the DEG or
%     KM options so that the station location can be found.  There is no
%     default value for this option.
%
%     TT=TAUPPATH(...,'B|BAZ|BAZI|BACKAZIMUTH',BAZI,...) sets the
%     backazimuth from the seismic station to the earthquake.  This option
%     needs to be paired with one of the DEG and KM options plus the STA
%     option to allow for calculating the earthquake location.  There is no
%     default value for this option.
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
%      tauppath('p','ttp','d',25)
%
%     Some other valid examples:
%      tauppath('dep',50,'ph','P,S','deg',45.6)
%      tauppath('mod','prem','dep',50,'ph','Pdiff,PKP',...
%          'sta',[-40 -100],'evt',[30,50])
%
%    See also: taup, tauptime, taupcurve, tauppierce

%     Version History:
%        Aug. 30, 2009 - major revision of script, name change to avoid
%                        breakage due to input/output changes
%        Sep.  1, 2009 - change default phase to ttbasic, add notes on
%                        input arg order preservation, add entries for az
%                        and baz, add figure name and title for figures
%        Sep.  2, 2009 - allow multiple calls to PHASES
%        Sep.  5, 2009 - minor doc update
%        Sep. 30, 2009 - changed abssawmod to abslatmod
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 30, 2009 at 15:40 GMT

% todo:

% check nargin
if(mod(nargin,2))
    error('matTaup:tauppath:badNumOptions','Unpaired option(s)!');
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
phase=cell(0); np1=0;
pargs=cell(0); np2=0;
dargs=cell(0); nd=0;
d=false; k=d; s=d; e=d; a=d; b=d;
for i=1:2:nargin
    switch lower(varargin{i})
        case {'m' 'mod' 'model'}
            if(isempty(varargin{i+1})); continue; end
            if(~ischar(varargin{i+1}))
                error('matTaup:tauppath:badInput',...
                    'MODEL must be a string (like ''prem'')!');
            end
            model=varargin{i+1}(:)';
        case {'h' 'z' 'dep' 'evdp' 'depth'}
            if(isempty(varargin{i+1})); continue; end
            if(~isscalar(varargin{i+1}) ...
                    || ~isreal(varargin{i+1}) || varargin{i+1}<0)
                error('matTaup:tauppath:badInput',...
                    'DEPTH must be a positive number (in km)!');
            end
            depth=num2str(varargin{i+1});
        case {'p' 'ph' 'phases'}
            if(isempty(varargin{i+1})); continue; end
            if(~ischar(varargin{i+1}))
                error('matTaup:tauppath:badInput',...
                    'PHASES must be a string (like ''P,S'')!');
            end
            np2=np2+2; np1=np1+1;
            pargs(1,np2-1:np2)={'-ph' varargin{i+1}};
            phase{np1}=varargin{i+1}(:)';
        case {'d' 'deg' 'gcarc' 'degrees'}
            if(isempty(varargin{i+1})); continue; end
            if(~isscalar(varargin{i+1}) ...
                    || ~isreal(varargin{i+1}) || varargin{i+1}<0)
                error('matTaup:tauppath:badInput',...
                    'DEG must be a positive number (in degrees)!');
            end
            nd=nd+2;
            dargs(1,nd-1:nd)={'-deg' num2str(varargin{i+1})};
            d=true;
            dist=varargin{i+1};
            sdist=[num2str(dist) 'deg'];
        case {'k' 'km' 'dist' 'kilometers'}
            if(isempty(varargin{i+1})); continue; end
            if(~isscalar(varargin{i+1}) ...
                    || ~isreal(varargin{i+1}) || varargin{i+1}<0)
                error('matTaup:tauppath:badInput',...
                    'KM must be a positive number (in kilometers)!');
            end
            nd=nd+2;
            dargs(1,nd-1:nd)={'-km' num2str(varargin{i+1})};
            k=true;
            dist=varargin{i+1};
            sdist=[num2str(dist) 'km'];
        case {'s' 'st' 'sta' 'station'}
            if(isempty(varargin{i+1})); continue; end
            if(numel(varargin{i+1})~=2 ...
                    || ~all(isreal(varargin{i+1})))
                error('matTaup:tauppath:badInput',...
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
                error('matTaup:tauppath:badInput',...
                    'EVT must be 2 real numbers ([LAT LON] in degrees)!');
            end
            nd=nd+3;
            dargs(1,nd-2:nd)={'-evt' num2str(varargin{i+1}(1)) ...
                num2str(varargin{i+1}(2))};
            e=true;
        case {'a' 'az' 'azi' 'azimuth'}
            if(isempty(varargin{i+1})); continue; end
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}))
                error('matTaup:tauppath:badInput',...
                    'AZ must be a scalar number (in degrees)!');
            end
            nd=nd+2;
            dargs(1,nd-1:nd)={'-az' num2str(varargin{i+1})};
            a=true;
        case {'b' 'baz' 'bazi' 'backazimuth'}
            if(isempty(varargin{i+1})); continue; end
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}))
                error('matTaup:tauppath:badInput',...
                    'BAZ must be a scalar number (in degrees)!');
            end
            nd=nd+2;
            dargs(1,nd-1:nd)={'-baz' num2str(varargin{i+1})};
            b=true;
        otherwise
            error('matTaup:tauppath:badOption',...
                'Unknown Option: %s',varargin{i});
    end
end

% set up inputs
inArgs{1}='-mod';
inArgs{2}=model;
inArgs{3}='-h';
inArgs{4}=depth;
if(np1>0)
    inArgs=[inArgs pargs];
    phases=char(strcat(phase,','))';    % comma delimited
    phases=phases(:)';                  % assure row vector
    phases=phases(phases~=32);          % remove spaces
    phases=phases(1:end-1);             % trim off trailing comma
else
    inArgs{5}='-ph';
    inArgs{6}=phases;
end
if(d || k || (s && e) || (s && b && (d || k)) || (e && a && (d || k)))
    inArgs=[inArgs dargs];
else
    % default to random distance
    disp('Distance not given or could not be determined - using random!')
    dist=rand*180;
    d=true;
    sdist=[num2str(dist) 'deg'];
    inArgs=[inArgs dargs {'-deg' num2str(dist)}];
end

% debug
%disp(inArgs);
%arrivals=MatTauP_Path.run_path(inArgs);

% attempt run
try
    arrivals=MatTauP_Path.run_path(inArgs);
catch
    % oops!
    error('matTaup:tauppath:runFailed',...
        ['Java exception occurred! Please check input options and\n'...
        'make sure your classpath.txt file includes matTaup.jar!']);
end

% get distance
if(~isempty(arrivals))
    dist=arrivals(1).getDistDeg;
    sdist=[num2str(dist) 'deg'];
elseif(~d && ~k)
    dist=nan;
    sdist=[num2str(dist) ' deg'];
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
    
    % initialize plot
    fh=figure('color','k','name','TauP Ray Paths');
    hold on
    
    % plot grid
    [cx,cy]=circle(6871);
    figure(fh);
    plot(cx,cy,'color',[0.2 0.2 0.2],'linewidth',2)
    set(gca,'position',[0.025 0.05 0.95 0.9])
    [cx1,cy1]=circle2(6871,180);
    [cx2,cy2]=circle2(6771,180);
    figure(fh);
    plot([cx1; cx2],[cy1; cy2],'color',[0.2 0.2 0.2],'linewidth',1)
    [cx1,cy1]=circle2(6871,36);
    [cx2,cy2]=circle2(6671,36);
    figure(fh);
    plot([cx1; cx2],[cy1; cy2],'color',[0.2 0.2 0.2],'linewidth',2)
    [cx,cy]=circle2(6371,4);
    figure(fh);
    plot(cx(1:2:3),cy(1:2:3),'color',[0.2 0.2 0.2],'linewidth',2)
    plot(cx(2:2:4),cy(2:2:4),'color',[0.2 0.2 0.2],'linewidth',2)
    
    % plot major discontinuities
    for i=[6371 3480 1220]
        [cx,cy]=circle(i);
        figure(fh);
        plot(cx,cy,'w','linewidth',2)
    end
    
    % plot minor discontinuities
    for i=[5961 5711 3780]
        [cx,cy]=circle(i);
        figure(fh);
        plot(cx,cy,'color',[0.5 0.5 0.5],'linewidth',1)
    end
    axis off
    axis equal
    
    % set default earthquake,station location
    cx=6371*sin((-45+[0 dist])/180*pi);
    cy=6371*cos((-45+[0 dist])/180*pi);
    
    % loop over phases
    colors=hsv(arrivals.length);
    ph=nan(1,arrivals.length); pn=cell(1,arrivals.length);
    for ii=1:arrivals.length
        % list phase info
        fprintf(' %7.2f  %6.1f   %-10s   %7.2f   %7.3f    %7.2f  = %-10s\n',...
            abslatmod(arrivals(ii).getDistDeg,180),...
            arrivals(ii).getSourceDepth,...
            char(arrivals(ii).getName),arrivals(ii).getTime,...
            arrivals(ii).getRayParam/R2D,arrivals(ii).getDistDeg,...
            char(arrivals(ii).getPuristName));
        
        % plot phase path
        temp.path.depth=arrivals(ii).getMatPath.depth;
        temp.path.distance=arrivals(ii).getMatPath.dist;
        cx=(6371-temp.path.depth).*sin((-45+temp.path.distance)/180*pi);
        cy=(6371-temp.path.depth).*cos((-45+temp.path.distance)/180*pi);
        figure(fh);
        ph(ii)=plot(cx,cy,'color',colors(ii,:));
        pn{ii}=char(arrivals(ii).getName);
    end
    
    % plot event and station markers
    figure(fh);
    h1=plot(cx(1),cy(1),'yp','markersize',10,'markerfacecolor','y');
    h2=plot(cx(end),cy(end),'r^','markersize',8,'markerfacecolor','r');
    hold off
    
    % legend
    lh=legend([h1 h2 ph],[{'event' 'station'} pn],...
        'location','westoutside','textcolor','w');
    set(lh,'color','none','edgecolor','w',...
        'fontsize',6,'interpreter','none')
    
    % title
    title({['MODEL: ' model '  PHASES: ' phases]...
        ['EVENT DEPTH: ' depth 'km  DISTANCE: ' ...
        sdist]},'color','w','interpreter','none')
    
    return
end

% struct output
tt(1:arrivals.length)=struct('phase',[],'distance',[],'depth',[],...
    'time',[],'rayparameter',[],'path',[]);
for ii=1:arrivals.length
    tt(ii).time=arrivals(ii).getTime;
    tt(ii).distance=arrivals(ii).getDistDeg;
    tt(ii).depth=arrivals(ii).getSourceDepth;
    tt(ii).phase=char(arrivals(ii).getName);
    tt(ii).rayparameter=arrivals(ii).getRayParam/R2D;
    tt(ii).path.rayparameter=arrivals(ii).getMatPath.p;
    tt(ii).path.time=arrivals(ii).getMatPath.time;
    tt(ii).path.distance=arrivals(ii).getMatPath.dist;
    tt(ii).path.depth=arrivals(ii).getMatPath.depth;
    tt(ii).path.latitude=arrivals(ii).getMatPath.lat;
    tt(ii).path.longitude=arrivals(ii).getMatPath.lon;
end

end

function [cx,cy]=circle(r)
    ang=0:0.002:pi*2;
    cx=sin(ang)*r;
    cy=cos(ang)*r;
end

function [cx,cy]=circle2(r,steps)
    ang=0:2*pi/steps:pi*2;
    cx=sin(-pi/4+ang)*r;
    cy=cos(-pi/4+ang)*r;
end

function [c]=abslatmod(a,b)
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
