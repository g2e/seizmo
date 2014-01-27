function [varargout]=cmb_1st_pass(phase,indir,varargin)
%CMB_1ST_PASS    Initial (wideband) core-diff rel. arrivals + amplitudes
%
%    Usage:    results=cmb_1st_pass(phase,indir)
%              results=cmb_1st_pass(phase,indir,'option',value,...)
%
%    Description:
%     RESULTS=CMB_1ST_PASS(PHASE,INDIR) provides an interface for the
%     initial alignment of core-diffracted waveforms through a series of
%     menus and plots.  The idea of this is to reduce the time shifts and
%     to figure out the relative polarities prior to narrow-band alignment
%     which has a diminished ability to resolve large time shifts and
%     relative polarities due to their highly oscillatory nature.  It is
%     expected that the data in INDIR was prepared with PREP_CMB_DATA.
%     PHASE is one of the following:
%      'Pdiff', 'SHdiff', 'SVdiff'.
%     If one or both input arguments are omitted, then an interface is
%     presented to graphically select them.
%
%     RESULTS=CMB_1ST_PASS(PHASE,INDIR,'OPTION',VALUE,...) sets specific
%     options pertinent to the processing.  The options directly used by
%     CMB_1ST_PASS are:
%      'odir'       - output directory (current directory by default)
%      'figdir'     - figure output directory (current dir by default)
%      'minradampl' - minimum radiation pattern amplitude (0.2 by default)
%      'delazcut'   - true/false, graphical dataset selection by dist/azi
%     Addtional option/value pairs are assumed to be MULTIBANDALIGN options
%     and so are passed on to it to adjust its parameters.  See
%     MULTIBANDALIGN for more details.
%
%    Notes:
%     - The data will be filtered between 25s & 80s.  This provides a
%       frequency range that is common to core-diffracted waves from 95 to
%       170 degrees from the earthquake which improves the correlation
%       between the waveforms.
%     - All figures and results are saved in the current directory with a
%       prefix of the input directory name.  Figures are saved as .fig
%       files while results are saved as .mat files for quick inspection in
%       Matlab/Octave.
%
%    Examples:
%     % Avoid copy/pasting the directory into Matlab by omitting the INDIR
%     % input and graphically selecting it instead:
%     results=cmb_1st_pass('Pdiff');
%
%     % Nothing is written to disk if the
%     % output directories are set to false:
%     results=cmb_1st_pass('SHdiff','my/in/dir',...
%                          'odir',false,'figdir',false);
%
%    See also: PREP_CMB_DATA, CMB_2ND_PASS, CMB_OUTLIERS, SLOWDECAYPAIRS,
%              SLOWDECAYPROFILES, MAP_CMB_PROFILES, CMB_CLUSTERING

%     Version History:
%        Dec. 12, 2010 - added docs
%        Jan. 15, 2011 - allow gui based selection if arguments not given,
%                        allow data directory selection
%        Jan. 16, 2011 - alterations for new uniform results struct
%        Jan. 18, 2011 - update for improved multibandalign, .time field
%        Jan. 23, 2011 - pre-align on waveform of interest, skip event if
%                        no waveforms
%        Jan. 26, 2011 - synthetics fields added (only reflectivity synth)
%        Jan. 29, 2011 - make earthmodel a string, datetime in front of
%                        output (to avoid overwrite), fix absolute path
%                        bug, fix phase bug
%        Jan. 31, 2011 - allow no output, odir catching
%        Feb.  2, 2011 - update for kuser0-2 containing split model name,
%                        fix odir bug
%        Mar. 10, 2011 - calculate corrections beforehand (to allow better
%                        estimated alignment), use cmt radiation to cut out
%                        near-nodal stations and set polarities, delaz cut
%                        option
%        Mar. 17, 2011 - corrections are not recalculated, fixed breakage
%                        when no waveforms found
%        Mar. 21, 2011 - bugfix for out-of-range stations
%        Mar. 24, 2011 - added detail msg for processing info
%        Mar. 25, 2011 - fix corrections desync from data
%        Apr.  6, 2011 - minor change to minradampl cut
%        Apr. 22, 2011 - update for finalcut field
%        May  20, 2011 - add some code to workaround matlab write bug
%        Mar.  1, 2012 - octave ascii save workaround
%        Mar.  5, 2012 - allow no written output, odir/figdir bugfix
%        Mar. 13, 2012 - use getheader improvements
%        Mar. 15, 2012 - fix for pick functions
%        Feb. 26, 2013 - bugfix: skip event if only 1 waveform
%                        bugfix: no crash when handling no good events
%        Mar. 11, 2013 - skip event if 2 or less waveforms
%        Jan. 27, 2014 - abs path fix & reduced filesep/fullfile calls
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 27, 2014 at 13:35 GMT

% todo:

% check nargin
if(nargin>2 && mod(nargin,2))
    error('seizmo:cmb_1st_pass:badNumInputs',...
        'OPTIONS must be paired with a value!');
end

% select/check phase
valid={'Pdiff' 'SHdiff' 'SVdiff'};
if(nargin<1 || isempty(phase))
    % prompt user for phase choice
    choice=[];
    while(isempty(choice))
        choice=menu('Align which phase?',valid{:});
    end
    phase=valid{choice};
elseif(~isstring(phase) || ~ismember(phase,valid))
    error('seizmo:cmb_1st_pass:badPhase',...
        ['PHASE must be one of the following:\n' ...
        sprintf('''%s'' ',valid{:}) '!']);
end

% directory separator
fs=filesep;

% select/check indir
if(nargin<2 || isempty(indir))
    drawnow;
    indir=uigetdir('.','Choose a directory of data directories:');
end
if(~isstring(indir))
    error('seizmo:cmb_1st_pass:badInput',...
        'INDIR must be a string giving one directory!');
end
if(~isabspath(indir)); indir=[pwd fs indir]; end
if(~isdir(indir))
    error('seizmo:cmb_1st_pass:badInput',...
        'INDIR must be a directory!');
end

% default/extract odir/minradampl/delazcut
odir='.';
minradampl=0.2;
cut_by_delaz=false;
figdir='';
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:cmb_1st_pass:badInput',...
        'OPTION must all be strings!');
end
keep=true(nargin-2,1);
for i=1:2:nargin-2
    switch lower(varargin{i})
        case {'outdir' 'odir'}
            odir=varargin{i+1};
            if(isempty(figdir))
                varargin{i}='figdir';
                figdir=odir;
            else
                keep(i:i+1)=false;
            end
        case 'minradampl'
            minradampl=varargin{i+1};
            keep(i:i+1)=false;
        case 'delazcut'
            cut_by_delaz=varargin{i+1};
            keep(i:i+1)=false;
        case 'figdir'
            figdir=varargin{i+1};
    end
end
varargin=varargin(keep);
if(isempty(figdir)); figdir='.'; end

% fix TRUE directory input
if(islogical(odir) && isscalar(odir) && odir); odir='.'; end
if(islogical(figdir) && isscalar(figdir) && figdir); figdir='.'; end

% check odir/figdir
if(~isstring(odir) && ~(islogical(odir) && isscalar(odir)))
    error('seizmo:cmb_1st_pass:badInput',...
        'ODIR must be a string or TRUE/FALSE!');
elseif(~isstring(figdir) && ~(islogical(figdir) && isscalar(figdir)))
    error('seizmo:cmb_1st_pass:badInput',...
        'FIGDIR must be a string or TRUE/FALSE!');
end
if(islogical(odir)); out=odir; else out=true; end
if(islogical(figdir)); figout=odir; else figout=true; end

% create odir if not there
if(out)
    [ok,msg,msgid]=mkdir(odir);
    if(~ok)
        warning(msgid,msg);
        error('seizmo:cmb_1st_pass:pathBad',...
            'Cannot create directory: %s',odir);
    end
elseif(~out && ~nargout)
    error('seizmo:cmb_1st_pass:badInput',...
        'Output variable must be assigned when no written output!');
end
if(figout)
    [ok,msg,msgid]=mkdir(figdir);
    if(~ok)
        warning(msgid,msg);
        error('seizmo:cmb_1st_pass:pathBad',...
            'Cannot create directory: %s',figdir);
    end
end

% check minimum radiation pattern amplitude
if(~isreal(minradampl) || ~isscalar(minradampl) || minradampl<0 ...
        || minradampl>1)
    error('seizmo:cmb_1st_pass:badInput',...
        'MINRADAMPL must be a scalar within the range 0 to 1 !');
end

% check delazcut
if(~islogical(cut_by_delaz) || ~isscalar(cut_by_delaz))
    error('seizmo:cmb_1st_pass:badInput',...
        'DELAZCUT must be true or false!');
end

% get date directories
dates=dir(indir);
dates(strcmp({dates.name},'.') | strcmp({dates.name},'..'))=[];
dates(~[dates.isdir])=[];
if(isempty(dates))
    % you probably selected a single directory to work on
    % so lets just use that as the date directory
    pathdirs=getwords(indir,fs);
    if(strcmp(indir(1),fs)); pathdirs{1}=[fs pathdirs{1}]; end
    indir=joinwords(pathdirs(1:end-1),fs);
    if(isempty(indir)); indir='.'; end
    dates=dir(indir);
    dates=dates(strcmp({dates.name}.',pathdirs(end)));
    s=1;
else
    % get user selected start date
    datelist=char(strcat({dates.name}.'));
    s=listdlg('PromptString','Select events:',...
              'InitialValue',1:numel(dates),...
              'ListSize',[170 300],...
              'ListString',datelist);
    %s=1:numel(dates);
end

% error if none selected
if(isempty(s))
    error('seizmo:cmb_1st_pass:noDirsSelected',...
        'No directories selected!');
end

% handle output for no good events
varargout={};

% loop over events
for i=1:numel(s)
    % echo directory name to screen
    disp(dates(s(i)).name);
    
    % read in headers
    data=readheader([indir fs dates(s(i)).name]);
    runname=[dates(s(i)).name '_' phase '_1stPass'];
    
    % appropriate components for this phase
    switch phase
        case 'Pdiff'
            % only vertical data
            data=data(vertcmp(data));
            truephase='Pdiff';
            disp(['Found ' num2str(numel(data)) ' Vertical Components']);
        case 'SHdiff'
            % only transverse component
            kcmpnm=getheader(data,'kcmpnm');
            data=data(strcmp(kcmpnm,'BHT') ...
                | strcmp(kcmpnm,'HHT') ...
                | strcmp(kcmpnm,'LHT'));
            truephase='Sdiff';
            disp(['Found ' num2str(numel(data)) ' Transverse Components']);
        case 'SVdiff'
            % only radial component
            kcmpnm=getheader(data,'kcmpnm');
            data=data(strcmp(kcmpnm,'BHR') ...
                | strcmp(kcmpnm,'HHR') ...
                | strcmp(kcmpnm,'LHR'));
            truephase='Sdiff';
            disp(['Found ' num2str(numel(data)) ' Radial Components']);
    end
    
    % skip if 2 or less waveforms
    if(numel(data)<3)
        warning('seizmo:cmb_1st_pass:noWaveforms',...
            'Too few %s waveforms found for event: %s',...
            phase,dates(s(i)).name)
        continue;
    end
    
    % sort by distance
    data=sortbyfield(data,'gcarc');
    
    % remove non-grazing/diffracted recordings
    gcarc=getheader(data,'gcarc');
    data(gcarc<90 | gcarc>170)=[];
    disp(['Found ' num2str(numel(data)) ' components in 90-170deg range']);
    
    % skip if 2 or less waveforms
    if(numel(data)<3)
        warning('seizmo:cmb_1st_pass:noWaveforms',...
            'Too few %s waveforms in 90-170deg range for event: %s',...
            phase,dates(s(i)).name)
        continue;
    end
    
    % get corrections
    corrections=cmb_corrections(phase,data);
    
    % remove stations near nodal planes
    good=abs(corrections.radpatcor)>=minradampl;
    data=data(good);
    corrections=fixcorrstruct(corrections,good);
    disp([num2str(sum(~good)) ' Near-Nodal Stations Removed']);
    
    % skip if 2 or less waveforms
    if(numel(data)<3)
        warning('seizmo:cmb_1st_pass:noWaveforms',...
            'Too few %s waveforms fit criteria for event: %s',...
            phase,dates(s(i)).name)
        continue;
    end
    
    % delazcut if desired
    if(cut_by_delaz)
        [ev,st]=getheader(data,'ev','st');
        [bad,azlim,ddlim,ax]=delazcut(ev(1,1:2),st(:,1:2),[],[],[],...
            'go',{'xaxislocation','middle'});
        data(bad)=[];
        corrections=fixcorrstruct(corrections,~bad);
        disp([num2str(sum(bad)) ...
            ' Stations Outside Distance-Azimuth Limits Removed']);
        if(ishandle(ax))
            if(figout)
                saveas(get(ax,'parent'),[figdir fs datestr(now,30) '_' ...
                    runname '_band_1_delazcut.fig']);
            end
            close(get(ax,'parent'));
        end
    else
        azlim=[0 360];
        ddlim=[0 180];
    end
    
    % skip if 2 or less waveforms
    if(numel(data)<3)
        warning('seizmo:cmb_1st_pass:noWaveforms',...
            'Too few %s waveforms fit criteria for event: %s',...
            phase,dates(s(i)).name)
        continue;
    end
    
    % time shift to phase
    t=findpicks(data,[truephase(1) ',' truephase],1);
    switch phase
        case 'Pdiff'
            t=t+corrections.ellcor ...
                +corrections.crucor.prem ...
                +corrections.mancor.hmsl06p.upswing;
        case {'SHdiff' 'SVdiff'}
            t=t+corrections.ellcor ...
                +corrections.crucor.prem ...
                +corrections.mancor.hmsl06s.upswing;
    end
    data=timeshift(data,-t);
    
    % check if all synthetics (only reflectivity synthetics for now)
    % - reflect2seizmo conventions here
    isynth=unique(getheader(data,'isynth id'));
    if(isscalar(isynth) && strcmpi(isynth,'ireflect'))
        issynth=true;
        [mod1,mod2,mod3]=getheader(data(1),'kuser0','kuser1','kuser2');
        synmodel=strtrim(char(strcat(mod1,mod2,mod3)));
    else
        issynth=false;
        synmodel='DATA';
    end
    
    % read data
    data=readdata(data);
    
    % correct polarity from radiation pattern
    data=multiply(data,sign(corrections.radpatcor));
    
    % set up a single bandpass for multibandalign
    bank=1./[42.5 80 25];
    
    % multibandalign
    tmp=multibandalign(data,...
        'phase',truephase,...
        'bank',bank,...
        'runname',runname,...
        'absxc',false,...
        'estarr',[],...
        'estpol',1,...
        'wgtpow',2,...
        'noisewin',[-125 -15],...
        'signalwin',[-10 70],...
        varargin{:});
    
    % matlab bug workaround
    % - really strange...saving fails to
    %   write all fields if we dont do this
    if(isempty(tmp.usersnr)); tmp.usersnr=[]; end
    
    % add run name, data directory, phase name
    tmp.runname=runname;
    tmp.dirname=[indir fs dates(s(i)).name];
    tmp.phase=phase;
    tmp.synthetics=issynth;
    tmp.earthmodel=synmodel;
    tmp.azlim=azlim;
    tmp.ddlim=ddlim;
    
    % add corrections
    if(~isempty(tmp.useralign))
        good=find(tmp.usersnr.snr>=tmp.usersnr.snrcut);
        good(tmp.userwinnow.cut)=[];
        if(isfield(tmp,'finalcut')); good=good(tmp.finalcut); end
        tmp.corrections=fixcorrstruct(corrections,good);
    else
        tmp.corrections=[];
    end
    
    % add in clustering info (all belong to 1 cluster)
    if(isempty(tmp.useralign))
        nrecs=0;
    else
        nrecs=numel(tmp.useralign.data);
    end
    tmp.usercluster.T=ones(nrecs,1);
    tmp.usercluster.units=zeros(nrecs,1);
    tmp.usercluster.good=true;
    tmp.usercluster.color=[1 0 0]; % default to red
    
    % add in outlier info (no outliers)
    tmp.outliers.bad=false(nrecs,1);
    
    % time of run
    tmp.time=datestr(now);
    timestr=datestr(now,30);
    
    % matlab bug workaround
    % - really strange...saving fails to
    %   write all fields if we dont do this
    if(isempty(tmp.userwinnow)); tmp.userwinnow=[]; end
    
    % save results
    if(out)
        if(isoctave)
            save([odir fs timestr '_' runname '_results.mat'],...
                '-7','-struct','tmp');
        else % matlab
            save([odir fs timestr '_' runname '_results.mat'],...
                '-struct','tmp');
        end
        
        % read in to check
        try
            tmp=load([odir fs timestr '_' runname '_results.mat']);
            error(check_cmb_results(tmp));
        catch
            warning('seizmo:cmb_1st_pass:failedWrite',...
                ['Had trouble writing RESULTS! ' ...
                'Check directory permissions!']);
        end
    end
    
    % export to command line too
    if(nargout); varargout{1}(i)=tmp; end
end

if(nargout && ~numel(varargout))
    error('seizmo:cmb_1st_pass:badIdea',...
        'It appears there was not enough data available for this phase!');
end

end


function [s]=fixcorrstruct(s,good)
fields=fieldnames(s);
for i=1:numel(fields)
    if(isstruct(s.(fields{i})) && isscalar(s.(fields{i})))
        s.(fields{i})=fixcorrstruct(s.(fields{i}),good);
    else
        s.(fields{i})=s.(fields{i})(good);
    end
end
end
