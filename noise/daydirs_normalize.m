function []=daydirs_normalize(indir,outdir,eqband,pband,tsw,fsw,o)
%DAYDIRS_NORMALIZE    Temporal & spectral normalization of day dir records
%
%    Usage:    daydirs_normalize(indir,outdir)
%              daydirs_normalize(indir,outdir,eqband,pband)
%              daydirs_normalize(indir,outdir,eqband,pband,tsw,fsw)
%              daydirs_normalize(indir,outdir,eqband,pband,tsw,fsw,overwr)
%
%    Description: DAYDIRS_NORMALIZE(INDIR,OUTDIR) performs temporal and
%     spectral normalization on records within day directory layout INDIR.
%     The temporal normalization uses a 2-pass sliding-absolute-mean to
%     de-emphasize glitches and earthquakes.  The spectral normalization
%     applies a sliding-absolute-mean in the frequency domain with a window
%     width of 2mHz.  This is tuned for ambient noise studies of 1sps data.
%
%     DAYDIRS_NORMALIZE(INDIR,OUTDIR,EQBAND,PBAND) adjusts the earthquake
%     band and the passband of the filter for output.  The default EQBAND
%     is [1/100 1/15], which emphasizes teleseismic events.  The default
%     PBAND is [1/150 1/3] and is tuned for 1sps data.
%
%     DAYDIRS_NORMALIZE(INDIR,OUTDIR,EQBAND,PBAND,TSW,FSW) adjusts the
%     sliding window widths.  TSW is the time-domain width and is 75sec by
%     default.  FSW is the frequency domain width and is 2mHz by default.
%
%     DAYDIRS_NORMALIZE(INDIR,OUTDIR,EQBAND,PBAND,TSW,FSW,OVERWR) quietly
%     overwrites pre-existing records in OUTDIR when OVERWRITE is set to
%     TRUE.  By default OVERWRITE is FALSE.
%
%    Notes:
%
%    Header changes: DEP*
%
%    Examples:
%
%    See also: DAYDIRS_MERGECUT_25HRS, DAYDIRS_RESAMPLE, DAYDIRS_RINST,
%              DAYDIRS_CORRELATE, DAYDIRS_ROTCORR, DAYDIRS_STACKCORR,
%              DAYDIRS_MAKE

%     Version History:
%        June 21, 2010 - initial version
%        Aug. 19, 2010 - bit more useful warning messages
%        Sep. 21, 2010 - commented out parallel processing lines
%        Oct.  6, 2010 - catch rotate error (when no records)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct.  6, 2010 at 11:15 GMT

% todo:

% check nargin
error(nargchk(2,7,nargin));

% defaults
if(nargin<3 || isempty(eqband)); eqband=[1/100 1/15]; end
if(nargin<4 || isempty(pband)); pband=[1/150 1/3]; end
if(nargin<5 || isempty(tsw)); tsw=75; end
if(nargin<6 || isempty(fsw)); fsw=0.002; end
if(nargin<7 || isempty(o)); o=false; end
if(numel(eqband)~=2 || ~isreal(eqband) || any(eqband<=0))
    error('seizmo:daydirs_normalize:badInput',...
        'EQBAND must be a positive-real 2-element vector (in Hz)!');
end
if(numel(pband)~=2 || ~isreal(pband) || any(pband<=0))
    error('seizmo:daydirs_normalize:badInput',...
        'PBAND must be a positive-real 2-element vector (in Hz)!');
end
if(~isscalar(tsw) || ~isreal(tsw) || tsw<=0)
    error('seizmo:daydirs_normalize:badInput',...
        'TSW must be a positive-real scalar (in seconds)!');
end
if(~isscalar(fsw) || ~isreal(fsw) || fsw<=0)
    error('seizmo:daydirs_normalize:badInput',...
        'FSW must be a positive-real scalar (in Hz)!');
end
if(~isscalar(o) || ~islogical(o))
    error('seizmo:daydirs_normalize:badInput',...
        'OVERWRITE flag must be a scalar logical!');
end

% check directories
if(~ischar(indir) || ~isvector(indir))
    error('seizmo:daydirs_normalize:fileNotString',...
        'INDIR must be a string!');
end
if(~exist(indir,'dir'))
    error('seizmo:daydirs_normalize:dirConflict',...
        ['Input Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end
if(~ischar(outdir) || ~isvector(outdir))
    error('seizmo:daydirs_normalize:fileNotString',...
        'OUTDIR must be a string!');
end
if(exist(outdir,'file'))
    if(~exist(outdir,'dir'))
        error('seizmo:daydirs_normalize:dirConflict',...
            'Output Directory: %s\nIs a file!',outdir);
    end
    if(~o)
        fprintf('Output Directory: %s\nDirectory Exists!\n',outdir);
        reply=input('Overwrite? Y/N [N]: ','s');
        if(isempty(reply) || ~strncmpi(reply,'y',1))
            disp('Not overwriting!');
            return;
        end
        disp('Overwriting!');
    end
end

% directory separator
fs=filesep;

% parallel processing setup (8 instances)
%matlabpool(8);

% get year directories and day directories
dirs=xdir([indir fs]);
dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % unhidden dirs
years=str2double({dirs.name});
nyears=numel(years);
if(any(isnan(years)))
    error('seizmo:daydirs_normalize:badLayout',...
        'Improper directory layout!');
end
jdays=cell(size(years));
for i=1:nyears
    % get day directories
    dirs=xdir([indir fs num2str(years(i))]);
    dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % unhidden dir
    jdays{i}=str2double({dirs.name});
    if(any(isnan(jdays{i})))
        error('seizmo:daydirs_normalize:badLayout',...
            'Improper directory layout!');
    end
end

% verbosity (turn it off for the loop)
verbose=seizmoverbose(false);
if(verbose); disp('Temporally & spectrally normalizing record(s)'); end

% loop over years
for i=1:nyears
    % working year
    yr=years(i);
    syr=num2str(yr);
    disp(['PROCESSING YEAR ' syr]);
    
    % loop over days
    for j=1:numel(jdays{i})
    %parfor j=1:numel(jdays{i})
        % working julian day
        jday=jdays{i}(j);
        sjday=num2str(jday,'%03d');
        
        % detail message
        disp(['PROCESSING DAY ' syr '.' sjday]);
        
        % attempt normalization
        try
            % read verticals
            skip=false;
            try
                data=readseizmo([indir fs syr fs sjday fs '*LHZ*'],...
                    [indir fs syr fs sjday fs '*BHZ*']);
            catch
                warning('seizmo:noise:missingRecords',...
                    'No vertical records found for this day!');
                skip=true;
            end
            
            if(~skip)
                % setup
                delta=getheader(data(1),'delta');
                nyqhz=1/(2*delta);
                tnorm=[pband(1)/nyqhz (nyqhz-pband(2))/nyqhz];
                tw=ceil(tsw/delta);
                
                % time normalization (target glitches then quakes)
                weights=add(slidingabsmean(data,tw),eps);
                data=dividerecords(data,weights);
                weights=add(slidingabsmean(...
                    iirfilter(data,'bp','b','c',eqband,'o',4,'p',2),...
                    tw),eps);
                data=dividerecords(data,weights);
                
                % spectral normalization (2mhz sliding window)
                data=dft(data,'rlim');
                data=whiten(data,fsw);
                data=taper(data,tnorm,[],'gausswin',10);
                data=idft(data);
                data=taper(data,0.01); % 15min tapers
                
                % write verticals
                writeseizmo(data,'pathchange',{indir outdir});
            end
            
            % read horizontals
            skip=false;
            try
                data=readseizmo([indir fs syr fs sjday fs '*LHE*'],...
                    [indir fs syr fs sjday fs '*LHN*'],...
                    [indir fs syr fs sjday fs '*BHE*'],...
                    [indir fs syr fs sjday fs '*BHN*']);
            catch
                warning('seizmo:noise:missingRecords',...
                    'No horizontal records found for this day!');
                skip=true;
            end
            
            if(~skip)
                % setup
                delta=getheader(data(1),'delta');
                nyqhz=1/(2*delta);
                tnorm=[pband(1)/nyqhz (nyqhz-pband(2))/nyqhz];
                tw=ceil(tsw/delta);
                
                % rotate horz pairs to north and east
                % - this also removes unpaired and cuts pairs to match
                try
                    data=rotate(data,'to',0,'kcmpnm1','N','kcmpnm2','E');
                catch
                    tmp=lasterror;
                    warning(tmp.message);
                    continue;
                end
                
                % skip if none left (this shouldn't happen)
                if(~numel(data)); continue; end
                
                % horz cmp time normalization
                % - divide by max weight
                % - target glitches then quakes
                weights=add(slidingabsmean(data,tw),eps);
                weights=recordfun(...
                    @(x,y)max(x,y),weights(1:2:end),weights(2:2:end));
                data=dividerecords(data,weights([1:end; 1:end]));
                weights=add(slidingabsmean(...
                    iirfilter(data,'bp','b','c',eqband,'o',4,'p',2),...
                    tw),eps);
                weights=recordfun(...
                    @(x,y)max(x,y),weights(1:2:end),weights(2:2:end));
                data=dividerecords(data,weights([1:end; 1:end]));
                
                % spectral normalization
                % - normalize by mean smooth spectra of north & east
                % - roughly 2mhz sliding window
                data=dft(data,'rlim');
                amph=rlim2amph(data);
                amph=slidingmean(amph,ceil(fsw./getheader(data,'delta')));
                amph=seizmofun(amph,@(x)x(:,[1:2:end; 1:2:end])+eps);
                amph=recordfun(@(x,y)(x+y)/2,amph(1:2:end),amph(2:2:end));
                amph=changeheader(amph,'iftype','irlim');
                data=dividerecords(data,amph([1:end; 1:end]));
                data=taper(data,tnorm,[],'gausswin',10);
                data=idft(data);
                data=taper(data,0.01); % 15min tapers
                
                % write horizontals
                writeseizmo(data,'pathchange',{indir outdir});
            end
        catch
            % close pool & fix verbosity
            tmp=lasterror;
            warning(tmp.message);
            %matlabpool close;
            seizmoverbose(verbose);
            
            % ???
            error(lasterror);
        end
    end
end

% parallel processing takedown
%matlabpool close;
seizmoverbose(verbose);

end
