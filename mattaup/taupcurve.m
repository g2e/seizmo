function [varargout]=taupcurve(varargin)
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
%    Description:
%     TAUPCURVE(...) (no outputs, w/ or w/o inputs) displays a formatted
%     list of information on many seismic phases for an earthquake with a
%     depth of 0km.  The 1D Earth model utilized by default is IASP91 (see
%     the MODEL option to adjust this).  The default phase list is
%     'ttbasic' (equivalent to setting phases to BASIC in Brian Kennett's
%     TTIMES program) which lists many common phases by default.  See
%     option PHASES to adjust the phase list.
%
%     TT=TAUPCURVE(...) (1 output, w/ or w/o inputs) returns a structure
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
%     TT=TAUPCURVE(...,'M|MOD|MODEL',MODEL,...) sets the 1D Earth model to
%     MODEL.  Accepts a variety of common models like 'prem', 'iasp91',
%     'ak135'.  See the TauP program/documentation for more.  You may also
%     give a 1DMODEL struct like from CMB_1DMODEL_LIBRARY.  The default
%     model is 'iasp91'.
%
%     TT=TAUPCURVE(...,'H|Z|DEP|EVDP|DEPTH',DEPTH,...) sets the event depth
%     to DEPTH.  DEPTH should be in kilometers!  The default is 0.
%
%     TT=TAUPCURVE(...,'P|PH|PHASES',PHASES,...) sets the phase list to
%     PHASES.  PHASES must be a comma separated list of seismic phases.  
%     See the TauP documentation for conventions and valid phase names.
%     The default phase list is 'ttbasic'.
%
%     TT=TAUPCURVE(...,'RD|RDEG|REDDEG',VELO,...) applies a reduction
%     velocity of VELO degrees per second.  Note that this is in reciprocal
%     units of that of the ray parameter returned.  There is no default
%     value for this option.
%
%     TT=TAUPCURVE(...,'RK|RKM|REDKM',VELO,...) applies a reduction
%     velocity of VELO kilometers per second.  There is no default value
%     for this option.
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
%     % Travel time curves for several P phases from a 300km deep event:
%     taupcurve('mod','prem','dep',300,'ph','ttp+')
%
%     % Apply a reduction velocity of 25km/s
%     taupcurve('mod','prem','dep',300,'ph','ttp+','redkm',25)
%
%    See also: TAUPPATH, TAUPTIME, TAUPPIERCE, TAUPCREATE, PLOT_TAUPCURVE

%     Version History:
%        Sep.  2, 2009 - major revision of script, name change to avoid
%                        breakage due to input/output changes, multiple
%                        calls to PHASES allowed
%        Sep.  5, 2009 - minor doc update
%        Nov. 13, 2009 - dropped some import calls
%        Jan.  6, 2011 - add matTaup.jar to dynamic java classpath if
%                        necessary
%        Feb. 24, 2012 - switch to enhanced taup (needed for output)
%        Feb. 27, 2012 - skip plotting if no arrivals
%        Jan. 26, 2014 - no longer need to update jar filenames
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 26, 2014 at 17:15 GMT

% todo:

% check nargin
if(mod(nargin,2))
    error('TauP:taupcurve:badNumOptions','Unpaired option(s)!');
end

% try adding *.jar if no TauP class exists
if(~exist('edu.sc.seis.TauP.MatTauP_Curve','class'))
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
rdeg=[]; rkm=[];

% check options
if(~iscellstr(varargin(1:2:end)))
    error('TauP:taupcurve:badInput',...
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
        case {'rd' 'rdeg' 'reddeg'}
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}))
                error('TauP:taupcurve:badInput',...
                    'REDDEG must be a scalar number (in deg/sec)!');
            end
            rdeg=varargin{i+1};
        case {'rk' 'rkm' 'redkm'}
            if(~isscalar(varargin{i+1}) || ~isreal(varargin{i+1}))
                error('TauP:taupcurve:badInput',...
                    'REDKM must be a scalar number (in km/sec)!');
            end
            rkm=varargin{i+1};
        otherwise
            error('TauP:taupcurve:badOption',...
                'Unknown Option: %s',varargin{i});
    end
end

% check reduction
if(~isempty(rdeg) && ~isempty(rkm))
    error('TauP:taupcurve:badInput',...
        'REDKM & REDDEG cannot both be specified!');
end

% create time object for velocity model
tcobj=javaObject('edu.sc.seis.TauP.MatTauP_Curve',model);

% calculate curves for the specified depth & phase
tcobj.setSourceDepth(depth);
tcobj.setPhaseNames(phases);
if(~isempty(rdeg))
    tcobj.setReduceTime(true);
    tcobj.setReduceVelDeg(rdeg);
elseif(~isempty(rkm))
    tcobj.setReduceTime(true);
    tcobj.setReduceVelKm(rkm);
end
tcobj.depthCorrect(depth);
tcobj.curvecalculate(0); % the input is meaningless
narr=tcobj.getNumPhases; % note it is NOT getNumArrivals

% conversion
R2D=180/pi;

% struct output
ok=false(narr,1);
tt(1:narr,1)=struct('modelname',[],'depth',[],'phase',[],...
    'puristphase',[],'distance',[],'mindistance',[],'time',[],...
    'rayparameter',[]);
for ii=1:narr
    % get curve info
    arr=tcobj.getMatCurve(ii-1);
    ok(ii)=arr.getNumPoints>0;
    tt(ii).modelname=modelname;
    tt(ii).depth=arr.getSourceDepth;
    tt(ii).phase=char(arr.getPhaseName);
    tt(ii).puristphase=char(arr.getPuristPhaseName);
    
    % Handle Octave/Matlab difference
    if(exist('OCTAVE_VERSION','builtin')==5) % Octave
        % preallocate
        n=arr.getNumPoints;
        [tt(ii).distance,tt(ii).mindistance,...
            tt(ii).time,tt(ii).rayparameter]=deal(nan(n,1));
        
        % loop over points
        for jj=1:n
            tt(ii).distance(jj)=arr.getDistance(jj-1);
            tt(ii).mindistance(jj)=arr.getMinDistance(jj-1);
            tt(ii).time(jj)=arr.getTime(jj-1);
            tt(ii).rayparameter(jj)=arr.getRayParam(jj-1);
        end
        tt(ii).rayparameter=tt(ii).rayparameter/R2D;
    else % Matlab
        tt(ii).distance=arr.getDistances;
        tt(ii).mindistance=arr.getMinDistances;
        tt(ii).time=arr.getTimes;
        tt(ii).rayparameter=arr.getRayParams/R2D;
    end
end
tt(~ok)=[];
narr=sum(ok);

% formatted listing
if(nargout==0)
    % header
    disp(' ')
    disp(['Model: ' modelname])
    disp('Depth   Phase              Min/Max Travel   Min/Max Ray Param   Min/Max Distance')
    disp(' (km)   Name                  Time (s)          p (s/deg)            (deg)      ')
    disp('--------------------------------------------------------------------------------')
    
    % loop over phases
    for ii=1:narr
        % list phase info
        fprintf(' %6.1f   %-10s   %8.2f /%8.2f   %7.3f /%7.3f    %7.2f /%7.2f\n',...
            tt(ii).depth,tt(ii).phase,...
            min(tt(ii).time),max(tt(ii).time),...
            min(tt(ii).rayparameter),max(tt(ii).rayparameter),...
            min(tt(ii).distance),max(tt(ii).distance));
    end
    
    % plot curves
    if(~isempty(tt))
        ax=plot_taupcurve(tt);
        title(ax(1),{['MODEL: ' modelname ...
            '  EVENT DEPTH: ' num2str(depth) 'km'] ...
            ['  PHASES: ' joinwords(phases,',')]},...
            'color','w','interpreter','none');
        title(ax(2),{['MODEL: ' modelname ...
            '  EVENT DEPTH: ' num2str(depth) 'km'] ...
            ['  PHASES: ' joinwords(phases,',')]},...
            'color','w','interpreter','none');
    end
else
    varargout{1}=tt;
end

end
