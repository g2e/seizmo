function [results]=cmb_1st_pass(phase,indir,varargin)
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
%     RESULTS=CMB_1ST_PASS(PHASE,INDIR,'OPTION',VALUE,...) passes
%     option/value pairs to MULTIBANDALIGN to adjust its parameters.  See
%     MULTIBANDALIGN for more details.
%
%    Notes:
%     - The data will be filtered between 25s & 80s.  This provides a
%       frequency range that is common to core-diffracted waves from 90 to
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
%    See also: PREP_CMB_DATA, CMB_2ND_PASS, CMB_OUTLIERS, SLOWDECAYPAIRS,
%              SLOWDECAYPROFILES, MAP_CMB_PROFILES

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
%                        output (to avoid overwrite), fix absolute path bug
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 29, 2011 at 13:35 GMT

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

% select/check indir
if(nargin<2 || isempty(indir))
    drawnow;
    indir=uigetdir('.','Choose a directory of data directories:');
elseif(~isstring(indir))
    error('seizmo:cmb_1st_pass:badInput',...
        'INDIR must be a string giving one directory!');
elseif(~isdir(indir))
    error('seizmo:cmb_1st_pass:badInput',...
        'INDIR must be a directory!');
end

% get date directories
dates=dir(indir);
dates(strcmp({dates.name},'.') | strcmp({dates.name},'..'))=[];
dates(~[dates.isdir])=[];
if(isempty(dates))
    % you probably selected a single directory to work on
    % so lets just use that as the date directory
    pathdirs=getwords(indir,filesep);
    if(strcmp(indir(1),filesep)); pathdirs{1}=[filesep pathdirs{1}]; end
    indir=joinwords(pathdirs(1:end-1),filesep);
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
end

% error if none selected
if(isempty(s))
    error('seizmo:cmb_1st_pass:noDirsSelected',...
        'No directories selected!');
end

% loop over events
for i=1:numel(s)
    % echo directory name to screen
    disp(dates(s(i)).name);
    
    % read in headers
    data=readheader([indir filesep dates(s(i)).name]);
    runname=[dates(s(i)).name '_' phase '_1stPass'];
    
    % appropriate components for this phase
    switch phase
        case 'Pdiff'
            % only vertical data
            data=data(vertcmp(data));
            truephase='Pdiff';
        case 'SHdiff'
            % only transverse component
            kcmpnm=getheader(data,'kcmpnm');
            data=data(strcmp(kcmpnm,'BHT') ...
                | strcmp(kcmpnm,'HHT') ...
                | strcmp(kcmpnm,'LHT'));
            truephase='Sdiff';
        case 'SVdiff'
            % only radial component
            kcmpnm=getheader(data,'kcmpnm');
            data=data(strcmp(kcmpnm,'BHR') ...
                | strcmp(kcmpnm,'HHR') ...
                | strcmp(kcmpnm,'LHR'));
            truephase='Sdiff';
    end
    
    % skip if no waveforms
    if(~numel(data))
        warning('seizmo:cmb_1st_pass:noWaveforms',...
            'No %s waveforms found for event: %s',phase,dates(s(i)).name)
        continue;
    end
    
    % time shift to phase
    [t,n]=getarrival(data,truephase);
    data=timeshift(data,-t,strcat('it',num2str(n)));
    
    % check if all synthetics (only reflectivity synthetics for now)
    % - reflect2seizmo conventions here
    isynth=unique(getenumid(data,'isynth'));
    if(isscalar(isynth) && strcmpi(isynth,'ireflect'))
        issynth=true;
        synmodel=char(getheader(data(1),'kuser2'));
    else
        issynth=false;
        synmodel='DATA';
    end
    
    % read data
    data=readdata(data);
    
    % sort by distance
    data=sortbyfield(data,'gcarc');
    
    % set up a singleband for multibandalign
    bank=1./[42.5 80 25];
    
    % multibandalign
    tmp=multibandalign(data,...
        'phase',truephase,...
        'bank',bank,...
        'runname',runname,...
        'absxc',true,...
        'estarr',[],...
        'wgtpow',2,varargin{:});
    
    % matlab bug workaround
    % - really strange...saving fails to
    %   write all fields if we dont do this
    if(isempty(tmp.usersnr)); tmp.usersnr=[]; end
    
    % add run name, data directory, phase name
    tmp.runname=runname;
    tmp.dirname=[indir filesep dates(s(i)).name];
    tmp.phase=phase;
    tmp.synthetics=issynth;
    tmp.earthmodel=synmodel;
    
    % add corrections
    if(~isempty(tmp.useralign))
        tmp.corrections=cmb_corrections(phase,tmp.useralign.data);
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
    
    % save results
    save([timestr '_' runname '_results.mat'],'-struct','tmp');
    
    % read in to check
    try
        tmp=load([timestr '_' runname '_results.mat']);
        error(check_cmb_results(tmp));
    catch
        warning('seizmo:cmb_1st_pass:failedWrite',...
            'Had trouble writing RESULTS! Check directory permissions!');
    end
    
    % export to command line too
    results(i)=tmp;
end

if(~exist('results','var'))
    error('seizmo:cmb_1st_pass:badIdea',...
        'It appears there was no data available for this phase!');
end

end
