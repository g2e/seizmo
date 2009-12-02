function tt=taupcurve(varargin)
%TAUPCURVE    Calculate travel time curves using the TauP toolkit
%
%    Usage:    taupcurve(...)
%           tt=taupcurve(...)
%           tt=taupcurve(...,'m|mod|model',model,...)
%           tt=taupcurve(...,'h|z|dep|evdp|depth',depth,...)
%           tt=taupcurve(...,'p|ph|phases',phases,...)
%           tt=taupcurve(...,'rd|rdeg|reddeg',velo,...)
%           tt=taupcurve(...,'rk|rkm|redkm',velo,...)
%
%    Description: TAUPCURVE(...) (no outputs, w/ or w/o inputs) displays a
%     formatted list of information on many seismic phases for an
%     earthquake with a depth of 0km.  The 1D Earth model utilized by
%     default is IASP91 (see the MODEL option to adjust this).  The default
%     phase list is 'ttbasic' (equivalent to setting phases to BASIC in
%     Brian Kennett's TTIMES program) which lists many common phases by
%     default.  See option PHASES to adjust the phase list.
%
%     TT=TAUPCURVE(...) (1 output, w/ or w/o inputs) returns a structure
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
%     TT=TAUPCURVE(...,'M|MOD|MODEL',MODEL,...) sets the 1D Earth model to
%     MODEL.  Accepts a variety of common models like 'prem', 'iasp91',
%     'ak135'.  See the TauP program/documentation for more.  The default
%     model is 'iasp91'.
%
%     TT=TAUPCURVE(...,'H|Z|DEP|EVDP|DEPTH',DEPTH,...) sets the event depth
%     to DEPTH.  DEPTH should be in kilometers!  The default is 0.
%
%     TT=TAUPCURVE(...,'P|PH|PHASES',PHASES,...) sets the phase list to
%     PHASES.  PHASES must be a comma separated list of seismic phases.  
%     See the TauP documentation for conventions and valid phase names.
%     Multiple calls are honored.  The default phase list is 'ttbasic'.
%
%     TT=TAUPCURVE(...,'RD|RDEG|REDDEG',VELO,...) applies a reduction
%     velocity of VELO degrees per second.  Note that this is in reciprocal
%     units of that of the ray parameter returned.  Also note that this is
%     only applied to the travel time, not the ray parameter.  There is no
%     default value for this option.
%
%     TT=TAUPCURVE(...,'RK|RKM|REDKM',VELO,...) applies a reduction
%     velocity of VELO kilometers per second.  This is only applied to the
%     travel time, but not the ray parameter (blame TauP).  There is no
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
%     Travel time curves for several P phases from a 300km deep event:
%      taupcurve('mod','prem','dep',300,'ph','ttp+')
%
%     Apply a reduction velocity of 25km/s
%      taupcurve('mod','prem','dep',300,'ph','ttp+','redkm',25)
%
%    See also: TAUP, TAUPPATH, TAUPTIME, TAUPPIERCE

%     Version History:
%        Sep.  2, 2009 - major revision of script, name change to avoid
%                        breakage due to input/output changes, multiple
%                        calls to PHASES allowed
%        Sep.  5, 2009 - minor doc update
%        Nov. 13, 2009 - dropped some import calls
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 13, 2009 at 17:25 GMT

% todo:

% check nargin
if(mod(nargin,2))
    error('matTaup:taupcurve:badNumOptions','Unpaired option(s)!');
end

% initialize java code
import edu.sc.seis.TauP.*;

% default options
model='iasp91';
depth='0';
phases='ttbasic';

% check options
phase=cell(0); np1=0;
pargs=cell(0); np2=0;
dargs=cell(0); nd=0; d=false; k=d;
for i=1:2:nargin
    switch lower(varargin{i})
        case {'m' 'mod' 'model'}
            if(isempty(varargin{i+1})); continue; end
            if(~ischar(varargin{i+1}))
                error('matTaup:taupcurve:badInput',...
                    'MODEL must be a string (like ''prem'')!');
            end
            model=varargin{i+1}(:)';
        case {'h' 'z' 'dep' 'evdp' 'depth'}
            if(isempty(varargin{i+1})); continue; end
            if(~isscalar(varargin{i+1}) ...
                    || ~isreal(varargin{i+1}) || varargin{i+1}<0)
                error('matTaup:taupcurve:badInput',...
                    'DEPTH must be a positive number (in km)!');
            end
            depth=num2str(varargin{i+1});
        case {'p' 'ph' 'phases'}
            if(isempty(varargin{i+1})); continue; end
            if(~ischar(varargin{i+1}))
                error('matTaup:taupcurve:badInput',...
                    'PHASES must be a string (like ''P,S'')!');
            end
            np2=np2+2; np1=np1+1;
            pargs(1,np2-1:np2)={'-ph' varargin{i+1}};
            phase{np1}=varargin{i+1}(:)';
        case {'rd' 'rdeg' 'reddeg'}
            if(isempty(varargin{i+1})); continue; end
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}))
                error('matTaup:taupcurve:badInput',...
                    'REDDEG must be a scalar number (in deg/sec)!');
            end
            nd=nd+2;
            dargs(1,nd-1:nd)={'-reddeg' num2str(varargin{i+1})};
            d=true;
        case {'rk' 'rkm' 'redkm'}
            if(isempty(varargin{i+1})); continue; end
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}))
                error('matTaup:taupcurve:badInput',...
                    'REDKM must be a scalar number (in km/sec)!');
            end
            nd=nd+2;
            dargs(1,nd-1:nd)={'-redkm' num2str(varargin{i+1})};
            k=true;
        otherwise
            error('matTaup:taupcurve:badOption',...
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
if(d || k); inArgs=[inArgs dargs]; end

% debug
%disp(inArgs);
%arrivals=MatTauP_Curve.run_curve(inArgs);

% attempt run
try
    arrivals=MatTauP_Curve.run_curve(inArgs);
catch
    % oops!
    error('matTaup:taupcurve:runFailed',...
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
    disp('Depth   Phase            Min/Max Travel   Min/Max Ray Param   Min/Max Distance')
    disp(' (km)   Name                Time (s)          p (s/deg)            (deg)      ')
    disp('------------------------------------------------------------------------------')
    
    % initialize dist vs time plot
    fh1=figure('color','k','name','TauP Travel Time Curves');
    pos=get(fh1,'position');
    set(fh1,'position',[pos(1) pos(2)-pos(4)*.75 pos(3) pos(4)*1.75]);
    set(gca,'xcolor','w','ycolor','w','color','k',...
        'yaxislocation','right','position',[0 0.1 0.9 0.8])
    title({['MODEL: ' model '  EVENT DEPTH: ' depth 'km']...
        ['  PHASES: ' phases]},'color','w','interpreter','none')
    xlabel('Distance (deg)')
    ylabel('Time (sec)')
    hold on
    box on
    
    % loop over phases
    colors=hsv(arrivals.length); n=0;
    ph1=nan(1,arrivals.length); ph2=ph1; pn=cell(1,arrivals.length);
    for ii=1:arrivals.length
        % list phase info
        fprintf(' %6.1f   %-10s   %8.2f/%8.2f   %7.3f/%7.3f    %7.2f/%7.2f\n',...
            arrivals(ii).sourceDepth,...
            char(arrivals(ii).phaseName),...
            min(arrivals(ii).time),max(arrivals(ii).time),...
            min(arrivals(ii).rayParam)/R2D,max(arrivals(ii).rayParam)/R2D,...
            min(arrivals(ii).dist),max(arrivals(ii).dist));
        
        if(numel(arrivals(ii).dist)>1)
            n=n+1;
            figure(fh1);
            ph1(n)=plot(arrivals(ii).dist,arrivals(ii).time,...
                'color',colors(ii,:));
            pn{n}=char(arrivals(ii).phaseName);
        end
    end
    
    % force distance range from 0 to 180
    xlim([0 180])
    
    % legend
    lh1=legend(ph1(1:n),pn(1:n),'location','westoutside');
    set(lh1,'color','none','edgecolor','w','textcolor','w',...
        'fontsize',6,'interpreter','none')
    
    % initialize dist vs ray parameter plot
    fh2=figure('color','k','name','TauP Ray Parameter Curves');
    pos=get(fh2,'position');
    set(fh2,'position',[pos(1) pos(2)-pos(4)*.75 pos(3) pos(4)*1.75]);
    set(gca,'xcolor','w','ycolor','w','color','k',...
        'yaxislocation','right','position',[0 0.1 0.9 0.8])
    title({['MODEL: ' model '  EVENT DEPTH: ' depth 'km']...
        ['  PHASES: ' phases]},'color','w','interpreter','none')
    xlabel('Distance (deg)')
    ylabel('Ray Parameter (sec/deg)')
    hold on
    box on
    
    % plot rayparameter vs dist
    n=0;
    for ii=1:arrivals.length
        if(numel(arrivals(ii).dist)>1)
            n=n+1;
            figure(fh2);
            ph2(n)=plot(arrivals(ii).dist,arrivals(ii).rayParam/R2D,...
                'color',colors(ii,:));
        end
    end
    
    % force distance range from 0 to 180
    xlim([0 180])
    
    % legend
    lh2=legend(ph2(1:n),pn(1:n),'location','westoutside');
    set(lh2,'color','none','edgecolor','w','textcolor','w',...
        'fontsize',6,'interpreter','none')
    
    return
end

% struct output
tt(1:arrivals.length)=struct('phase',[],'distance',[],'depth',[],...
    'time',[],'rayparameter',[]);
for ii=1:arrivals.length
    tt(ii).time=arrivals(ii).time;
    tt(ii).distance=arrivals(ii).dist;
    tt(ii).depth=arrivals(ii).sourceDepth;
    tt(ii).phase=char(arrivals(ii).phaseName);
    tt(ii).rayparameter=arrivals(ii).rayParam/R2D;
end

end
