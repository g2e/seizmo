function []=noise_process(indir,outdir,steps,varargin)
%NOISE_PROCESS    Performs processing of seismic data for noise analysis
%
%    Usage:    noise_process(indir,outdir)
%              noise_process(indir,outdir,steps)
%              noise_process(indir,outdir,steps,'opt1',val,...,'optN',val)
%
%    Description:
%     NOISE_PROCESS(INDIR,OUTDIR) processes the data under the directory
%     INDIR using noise cross correlation methods.  The resulting data are
%     written to OUTDIR.  The following techniques may be performed on the
%     input dataset (only steps 7, 8, 9, 11, 12, and 14 are performed by
%     default, see the following Usage forms to set an alternate processing
%     step list and to alter how each step is performed):
%      (* 1) remove dead records (no change in recorded value)
%      (* 2) remove short records (spanning less than 70% of time section)
%      (* 3) remove mean & trend
%      (* 4) taper (first/last 1%)
%      (* 5) resample to 1 sample/sec
%      (* 6) remove polezero response (displacement, hp taper: [.004 .008])
%      (  7) reject records based on amplitudes (>Inf nm)
%      (  8) rotate horizontals to North/East (removes unpaired)
%      (  9) t-domain normalize (equal scaling for each time section)
%      ( 10) f-domain normalize (2mHz moving average)
%      ( 11) interpolate/synchronize records to time section
%      ( 12) correlate (keep +/-4000s lagtime)
%      ( 13) correlate the coda of correlations (can be repeated)
%      ( 14) rotate correlations into <RR, RT, TR, TT>
%       * -> Steps 1-6 are done in NOISE_SETUP by default and should be
%            skipped if that function was already run on the dataset.
%
%     NOISE_PROCESS(INDIR,OUTDIR,STEPS) only does the processing steps
%     indicated in STEPS.  STEPS is a vector of numbers corresponding to
%     the valid steps given above.  The default is [7:9 11 12 14] (note
%     that 1:6 are done in NOISE_SETUP by default and should be skipped).
%
%     NOISE_PROCESS(INDIR,OUTDIR,STEPS,'OPT1',VAL,...,'OPTN',VAL) allows
%     changing some of the noise correlation parameters.  The following are
%     configurable:
%      MINIMUMLENGTH - minimum length of records in % of time section [70]
%      TAPERWIDTH    - % width of taper relative to record length [1]
%      TAPERTYPE     - type of taper []
%      TAPEROPTION   - option for taper []
%      SAMPLERATE    - samplerate to synchronize records to in step 5 [1]
%      PZDB          - polezero db to use in step 6 []
%      UNITS         - units of records after polezero removal ['disp']
%      PZTAPERLIMITS - highpass taper to stabilize pz removal [.004 .008]
%      REJECTMETHOD  - data rejection method: 'rms' or 'abs' ['abs']
%      REJECTCUTOFF  - data rejection cutoff: abs value or NxRMS [Inf]
%      REJECTWIDTH   - include 2N neighboring timewindows (+current) when
%                      calculating the rms value for data rejection [0]
%      TDSTYLE       - time domain normalization style:
%                      ['flat'] - single scaling for all in each timewindow
%                       '1bit'  - set amplitudes to +/-1
%                       'clip'  - clip values above some absolute value
%                       'ram'   - normalized using running-absolute mean
%      TDRMS         - use TDCLIP x RMS when TDSTYLE='clip' [true]
%      TDCLIP        - sets clip for TDSTYLE='clip' [1]
%      TDWEIGHTBANDS - frequency bands for TDSTYLE='ram' weights
%                     [1/150 1/3; 1/100 1/15]
%                     (TUNED TO REDUCE GLITCHES & TELESEISMIC EARTHQUAKES)
%      FDSTYLE       - frequency domain normalization style:
%                       '1bit' - set spectral ampl. to 1
%                      ['ram'] - normalized using running absolute mean
%      FDWIDTH       - frequency window width for FDSTYLE='ram' (Hz) [.002]
%      XCMAXLAG      - maximum lag time of correlograms in sec [4000]
%      NORMXC        - normalize the cross correlations [true]
%      COHERENCY     - return coherency correlograms [false]
%      FDOUT         - return correlograms in the frequency-domain [false]
%      NOAUTO        - skip autocorrelations [false]
%      C3VRAYL       - Rayleigh wave velocity for coda window. [2.6km/s]
%      C3WINOFF      - coda window offset [0s]
%      C3WINLEN      - coda window length [1200s]
%      C3STNLIST     - coda station list []
%                      Station list should be a cell array of strings as
%                      {'NET.STN.HOLE.CMP' ...}
%      C3XCCMP       - correlation coda component to analyze.
%                      'causal', 'acausal', 'both', ['symmetric']
%      C3ZTRANS      - stack coda xcorr w/ Fisher's transform? [true]
%      C3MINCODA     - minimum coda stations for coda correlation [1]
%      TIMESTART     - process time sections from this time on []
%      TIMEEND       - process times sections before this time []
%      LATRNG        - include stations in this latitude range []
%      LONRNG        - include stations in this longitude range []
%      NETWORKS      - include records with these network codes []
%      STATIONS      - include records with these station codes []
%      STREAMS       - include records with these stream codes []
%      COMPONENTS    - include records with these component codes []
%      FILENAMES     - limit processing to files with these filenames []
%      QUIETWRITE    - quietly overwrite OUTDIR (default is false)
%
%    Notes:
%     - Good Noise Analysis References:
%        Bensen et al 2007, GJI, doi:10.1111/j.1365-246X.2007.03374.x
%        Yang et al 2007, GJI, doi:10.1111/j.1365-246X.2006.03203.x
%        Lin et al 2008, GJI, doi:10.1111/j.1365-246X.2008.03720.x
%        Harmon et al 2008, GRL, doi:10.1029/2008GL035387
%        Stehly et al 2008, JGR, doi:10.1029/2008JB005693
%        Prieto et al 2009, JGR, doi:10.1029/2008JB006067
%        Ekstrom et al 2009, GRL, doi:10.1029/2009GL039131
%        Weaver 2011, CRG, doi:10.1016/j.crte.2011.07.001
%        Seats et al 2012, GJI, doi:10.1111/j.1365-246X.2011.05263.x
%        Ma & Beroza 2012, GRL, doi:10.1029/2011g1050755
%        Zhang & Yang 2013, JGR, doi:10.1002/jgrb.50186
%     - Some steps for horizontals currently require running seismogram
%       rotation on the same run (but if you forget it is automatically
%       done for you).
%     - While data i/o is .mat files by default, SAC files are okay too.
%       Note that .mat files can be read into Matlab using the LOAD
%       function and written out as SAC files using the function
%       WRITESEIZMO.  You may also convert the inputs & outputs using
%       NOISE_MAT2SAC & NOISE_SAC2MAT.  Be aware that SAC files are ignored
%       when there is a .mat file present.
%
%    Header changes: Varies with steps chosen...
%
%    Examples:
%     % Perform the first 3 steps of noise processing,
%     % writing out the resulting data of each step:
%     noise_process('raw','step1',1)
%     noise_process('step1','step2',2)
%     noise_process('step2','step3',3)
%     % This is great for prototyping and debugging!
%
%     % Skip the normalization steps:
%     noise_process('setup','ncfs',[7 8 11 12 14])
%
%     % Use non-overlapping 15-minute timesections, sampled
%     % at 5Hz to look at noise up to about 1.5Hz:
%     noise_setup('raw','15disp','l',15,'o',0,'sr',5);
%     noise_process('15disp','15xc5',[],'xcmaxlag',500);
%     % We adjusted the lag time b/c 15 minutes is 900 seconds.
%
%     % The default setup output makes 3 hour, non-overlapping timewindows.
%     % Here we use data rejection with the 'rms' method.  Data with a
%     % amplitude exceeding 100 times the rms for the surrounding day are
%     % rejected ('rejectwidth'=6 as (6x2+1) x 3 hours = 25 hours):
%     noise_rms_calc('3hr');
%     noise_process('3hr','ncfs',[],'rm','rms','rc',100,'rw',6);
%
%    See also: NOISE_SETUP, NOISE_STACK, NOISE_OVERVIEW, NOISE_RMS_CALC

%     Version History:
%        Nov. 22, 2011 - initial version (only first 6 steps)
%        Nov. 29, 2011 - more steps, subsetting
%        Dec. 13, 2011 - fix output writing
%        Jan. 12, 2012 - correlogram subsetting, normalization
%        Jan. 24, 2012 - doc update, parameter parsing/checking, several
%                        parameters altered
%        Jan. 25, 2012 - drop 1bit+ram, ram is now multipass capable
%        Jan. 27, 2012 - allow filenames to be a cellstr array, minor
%                        improvement to string option handling, split
%                        reading of header and data for speed
%        Feb.  1, 2012 - all time operations are in utc, turn off checking,
%                        parallelization edits, fdpassband option gone
%        Mar.  8, 2012 - skip dirs outside time limits before reading data,
%                        now checks steps, early h/v split & shortcircuits
%        Mar. 15, 2012 - bad parfor variable bugfix, parallel verbose fixes
%        Apr. 26, 2012 - minor doc fix (lp->hp), step 7 is now amplitude
%                        based rejection
%        May   3, 2012 - fixed bug in amplitude-based rejection
%        May  14, 2012 - set amp rejection to Inf (no reject) by default
%        May  31, 2012 - speed overrides for (divide/add)records
%        June  3, 2012 - immediately change class to doubles
%        June 11, 2012 - fix rms formula
%        June 14, 2012 - make time options names more flexible
%        Aug. 23, 2012 - td tapering after fd norm (removes high freq zero
%                        slowness noise), changed 1st tdnorm filter to be
%                        wider (captures glitches better), remove
%                        correlations against cmp for the same station
%        Sep. 20, 2012 - correlations are now unnormalized
%        Sep. 23, 2012 - horizontals are sorted before correlation now
%        Feb. 27, 2013 - support for new correlate (outputs auto
%                        correlations as well), rotate_correlations rename
%        Mar. 25, 2013 - matio by default, 7:12 by default, fixed mat file
%                        reading/writing, updated docs accordingly
%        Aug.  8, 2013 - rms-based data rejection has been added to step 7
%                        but requires running NOISE_RMS_CALC first, this
%                        adds several new reject* options too
%        Sep. 23, 2013 - update for new rotate_correlations, 4=>3 fix
%        Jan. 26, 2014 - abs path exist fix
%        May  27, 2014 - remove MATIO option: now autodetects sac/mat i/o
%        May  28, 2014 - minor warning fix
%        June  4, 2014 - c3 processing step (step 13)
%        June 13, 2014 - c3 & xc rotate steps switched, bugfix: xc rotate
%                        shortcircuit was too high, added c3 shortcircuits
%        June 25, 2014 - added coherency, fdout & noauto correlate options,
%                        c3 options now passed to noise_c3, flat tdstyle
%                        option (based on Weaver et al 2011)
%        July 14, 2014 - flat td normalization is default, shifted steps to
%                        allow for interpolation of seismograms rather than
%                        correlograms, made step order changing easier to
%                        maintain
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 14, 2014 at 11:15 GMT

% todo:
% - 3cmp support
%   - requires rotate3_correlations

% check nargin
error(nargchk(2,inf,nargin));
if(nargin>=4 && ~mod(nargin,2))
    error('seizmo:noise_process:badInput',...
        'Unpaired option/value pair given!');
end

% directory separator
fs=filesep;

% steps change occasionally so lets make this easier for maintenance
STEP.AR=7;   % amplitude based rejection
STEP.ROT=8;  % seismogram rotation
STEP.TDN=9;  % time-domain normalization
STEP.FDN=10; % frequency-domain normalization
STEP.INT=11; % seismogram interpolation
STEP.XC=12;  % cross correlation
STEP.C3=13;  % correlation coda correlation
STEP.XCR=14; % correlation rotation

% default steps
if(nargin<3 || isempty(steps)); steps=[7:9 11 12 14]; end

% force pre-steps for some
if(any(steps==STEP.TDN | steps==STEP.FDN | steps==STEP.XC))
    % no skipping rotation to North & East
    steps=[steps(:).' STEP.ROT];
end
if(any(steps==STEP.XC))
    % no skipping interpolation & synchronization
    steps=[steps(:).' STEP.INT];
end

% parse/check options
opt=noise_process_parameters(varargin{:});

% check directories
if(~isstring(indir))
    error('seizmo:noise_process:fileNotString',...
        'INDIR must be a string!');
end
if(~isabspath(indir)); indir=[pwd fs indir]; end
if(~exist(indir,'dir'))
    error('seizmo:noise_process:dirConflict',...
        ['Input Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end
if(~isstring(outdir))
    error('seizmo:noise_process:fileNotString',...
        'OUTDIR must be a string!');
end
if(~isabspath(outdir)); outdir=[pwd fs outdir]; end
if(exist(outdir,'file'))
    if(~exist(outdir,'dir'))
        error('seizmo:noise_process:dirConflict',...
            'Output Directory: %s\nIs a file!',outdir);
    end
    if(~opt.QUIETWRITE)
        fprintf('Output Directory: %s\nDirectory Exists!\n',outdir);
        reply=input('Overwrite? Y/N [N]: ','s');
        if(isempty(reply) || ~strncmpi(reply,'y',1))
            disp('Not overwriting!');
            return;
        end
        disp('Overwriting!');
    end
end

% check steps
if(isempty(steps) || ~isnumeric(steps) || ~isreal(steps) ...
        || any(steps<=0) || any(steps~=fix(steps)))
    error('seizmo:noise_process:badInput',...
        'STEPS must be a series of positive integers!');
end

% parallel processing setup
verbose=seizmoverbose;
%matlabpool(4); % PARALLEL

% get year directories and time-section directories
dirs=xdir([indir fs]);
dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % unhidden dirs
yrdir={dirs.name};
nyr=numel(yrdir);
tsdir=cell(size(yrdir));
for i=1:nyr
    dirs=xdir([indir fs yrdir{i}]);
    dirs=dirs([dirs.isdir]' & ~strncmp({dirs.name}','.',1)); % unhidden dir
    tsdir{i}={dirs.name};
end
tsdirs=[tsdir{:}]'; % all timesections
clear dirs yrdir nyr tsdir;

% convert time ranges to numeric arrays
t=char(tsdirs);
tsbgn=str2double([cellstr(t(:,1:4)) cellstr(t(:,6:8)) ...
    cellstr(t(:,10:11)) cellstr(t(:,13:14)) cellstr(t(:,16:17))]);
tsend=str2double([cellstr(t(:,19:22)) cellstr(t(:,24:26)) ...
    cellstr(t(:,28:29)) cellstr(t(:,31:32)) cellstr(t(:,34:35))]);
clear t;

% forget about any outside of user specified time limits
good=true(numel(tsdirs),1);
if(~isempty(opt.TIMESTART))
    good=good & timediff(opt.TIMESTART,tsend)>0;
end
if(~isempty(opt.TIMEEND))
    good=good & timediff(opt.TIMEEND,tsbgn)<0;
end
tsdirs=tsdirs(good); tsbgn=tsbgn(good,:); tsend=tsend(good,:);
clear good;

% rms based rejection setup
if(any(steps==STEP.AR)) % simple amplitude based rejection
    switch lower(opt.REJECTMETHOD)
        case 'rms'
            rmsname=[indir fs 'noise_rms_info.mat'];
            if(exist(rmsname,'file') && ~exist(rmsname,'dir'))
                rmsinfo=load(rmsname);
            else
                error('seizmo:noise_process:badInput',...
                    ['Cannot find noise_rms_info.mat!\n' ...
                    'Make sure to run NOISE_RMS_CALC on ' ...
                    'the input directory before running ' ...
                    'NOISE_PROCESS step ' num2str(STEP.AR) '!']);
            end
            if(any(~ismember({'data' 'tsbgn' 'tsend'},...
                    fieldnames(rmsinfo))))
                error('seizmo:noise_process:badInput',...
                    'noise_rms_info.mat is missing required RMS data!');
            end
            if(opt.REJECTWIDTH>0)
                rmsinfo.data=slidingrms(rmsinfo.data,opt.REJECTWIDTH);
            end
            kname=upper(getheader(rmsinfo.data,'kname'));
            rmsinfo.kname=strcat(kname(:,1),'.',kname(:,2),'.',...
                kname(:,3),'.',kname(:,4));
    end
end

% detail message
if(verbose); disp('PROCESSING SEISMIC DATA FOR NOISE ANALYSIS'); end

% loop over time section directories
%parfor i=1:numel(tsdirs) % PARALLEL
for i=1:numel(tsdirs) % SERIAL
    % force quietness even in parfor (which resets globals)
    seizmoverbose(false);
    
    % read the header
    try
        opt.MATIO=true;
        data=load(strcat(indir,fs,tsdirs{i}(1:4),fs,...
            tsdirs{i},fs,'noise_records.mat'));
        data=data.noise_records;
        if(~isempty(opt.FILENAMES))
            warning('seizmo:noise_stack:unusedOption',...
                'FILENAMES option ignored for MAT input!');
        end
    catch
        opt.MATIO=false;
        try
            data=readheader(strcat(indir,fs,tsdirs{i}(1:4),fs,tsdirs{i},...
                fs,opt.FILENAMES));
        catch
            % no data...
            continue;
        end
    end
    if(isempty(data)); continue; end
    
    % check if records are correlations
    [kuser0,kuser1]=getheader(data,'kuser0','kuser1');
    xc=ismember(kuser0,{'MASTER' 'SLAVE'}) ...
        & ismember(kuser1,{'MASTER' 'SLAVE'});
    if(all(xc)); isxc=true;
    elseif(all(~xc)); isxc=false;
    else
        error('seizmo:noise_process:mixedData',...
            ['Data contains both seismic records and ' ...
            'correlograms! This is NOT allowed.\nYou might try ' ...
            'using the FILENAMES option to limit input to one type.']);
    end
    if(isxc && any(steps<STEP.XC)) % CAREFUL
        error('seizmo:noise_process:invalidProcess4xcdata',...
            'Cannot run earlier processing steps on correlograms!');
    elseif(~isxc && ~any(steps==STEP.XC) && any(steps>STEP.XC))
        error('seizmo:noise_process:invalidProcess4data',...
            'Cannot run later processing steps on non-correlograms!');
    end
    
    % proceed by data type
    if(isxc) % correlograms
        % limit to the stations that the user allows
        % Note: check both fields!
        if(~isempty(opt.LATRNG))
            [stla,evla]=getheader(data,'stla','evla');
            data=data(stla>=min(opt.LATRNG) & stla<=max(opt.LATRNG) ...
                & evla>=min(opt.LATRNG) & evla<=max(opt.LATRNG));
            if(isempty(data)); continue; end
        end
        if(~isempty(opt.LONRNG))
            [stlo,evlo]=getheader(data,'stlo','evlo');
            data=data(stlo>=min(opt.LONRNG) & stlo<=max(opt.LONRNG) ...
                & evlo>=min(opt.LONRNG) & evlo<=max(opt.LONRNG));
            if(isempty(data)); continue; end
        end
        if(~isempty(opt.NETWORKS))
            [knetwk1,knetwk2]=getheader(data,'knetwk','kt0');
            data=data(ismember(lower(knetwk1),opt.NETWORKS) ...
                & ismember(lower(knetwk2),opt.NETWORKS));
            if(isempty(data)); continue; end
        end
        if(~isempty(opt.STATIONS))
            [kstnm1,kstnm2]=getheader(data,'kstnm','kt1');
            data=data(ismember(lower(kstnm1),opt.STATIONS) ...
                & ismember(lower(kstnm2),opt.STATIONS));
            if(isempty(data)); continue; end
        end
        if(~isempty(opt.STREAMS))
            [khole1,khole2]=getheader(data,'khole','kt2');
            data=data(ismember(lower(khole1),opt.STREAMS) ...
                & ismember(lower(khole2),opt.STREAMS));
            if(isempty(data)); continue; end
        end
        if(~isempty(opt.COMPONENTS))
            [kcmpnm1,kcmpnm2]=getheader(data,'kcmpnm','kt3');
            data=data(ismember(lower(kcmpnm1),opt.COMPONENTS) ...
                & ismember(lower(kcmpnm2),opt.COMPONENTS));
            if(isempty(data)); continue; end
        end
    else % seismic records
        % limit to stations user allowed
        if(~isempty(opt.LATRNG))
            stla=getheader(data,'stla');
            data=data(stla>=min(opt.LATRNG) & stla<=max(opt.LATRNG));
            if(isempty(data)); continue; end
        end
        if(~isempty(opt.LONRNG))
            stlo=getheader(data,'stlo');
            data=data(stlo>=min(opt.LONRNG) & stlo<=max(opt.LONRNG));
            if(isempty(data)); continue; end
        end
        if(~isempty(opt.NETWORKS))
            knetwk=getheader(data,'knetwk');
            data=data(ismember(lower(knetwk),opt.NETWORKS));
            if(isempty(data)); continue; end
        end
        if(~isempty(opt.STATIONS))
            kstnm=getheader(data,'kstnm');
            data=data(ismember(lower(kstnm),opt.STATIONS));
            if(isempty(data)); continue; end
        end
        if(~isempty(opt.STREAMS))
            khole=getheader(data,'khole');
            data=data(ismember(lower(khole),opt.STREAMS));
            if(isempty(data)); continue; end
        end
        if(~isempty(opt.COMPONENTS))
            kcmpnm=getheader(data,'kcmpnm');
            data=data(ismember(lower(kcmpnm),opt.COMPONENTS));
            if(isempty(data)); continue; end
        end
    end
    
    % splits data into vertical and horizontal sets
    if(~isxc)
        vdata=data(vertcmp(data));
        hdata=data(horzcmp(data));
        data=[]; % clearing data
        
        % shortcircuits
        if(any(steps==STEP.XC) && numel(vdata)==1); vdata=[]; end
        if(any(steps==STEP.XC) && numel(hdata)==1); hdata=[]; end
        if(any(steps==STEP.C3) && numel(vdata)<3); vdata=[]; end
        if(any(steps==STEP.C3) && numel(hdata)<6); hdata=[]; end
        if(any(steps==STEP.XCR) && numel(hdata)==1); hdata=[]; end
        if(isempty(hdata) && isempty(vdata)); continue; end
    else
        % shortcircuits
        if(any(steps==STEP.C3) && numel(data)<3); continue; end
        if(any(steps==STEP.XCR) && numel(data)<3); continue; end
    end
    
    % read in data
    if(~opt.MATIO)
        if(isxc)
            data=changeclass(readdata(data),'double');
        else
            if(~isempty(vdata))
                vdata=changeclass(readdata(vdata),'double');
            end
            if(~isempty(hdata))
                hdata=changeclass(readdata(hdata),'double');
            end
        end
    end
    
    % detail message
    if(verbose); disp(['PROCESSING: ' tsdirs{i}]); end
    
    % turn off checking
    oldseizmocheckstate=seizmocheck_state(false);
    oldcheckheaderstate=checkheader_state(false);
    
    try
        % process data for noise analysis
        if(any(steps==1)) % remove dead
            if(~isempty(vdata)); vdata=removedeadrecords(vdata); end
            if(~isempty(hdata)); hdata=removedeadrecords(hdata); end
            if(any(steps==STEP.XC) && numel(vdata)==1); vdata=[]; end
            if(any(steps==STEP.XC) && numel(hdata)==1); hdata=[]; end
            if(any(steps==STEP.C3) && numel(vdata)<3); vdata=[]; end
            if(any(steps==STEP.C3) && numel(hdata)<6); hdata=[]; end
            if(any(steps==STEP.XCR) && numel(hdata)==1); hdata=[]; end
            if(isempty(hdata) && isempty(vdata)); continue; end
        end
        if(any(steps==2)) % remove short
            if(~isempty(vdata))
                [b,e]=getheader(vdata,'b','e');
                vdata=vdata(...
                    e-b>opt.MINIMUMLENGTH*timediff(tsbgn(i,:),tsend(i,:)));
            end
            if(~isempty(hdata))
                [b,e]=getheader(hdata,'b','e');
                hdata=hdata(...
                    e-b>opt.MINIMUMLENGTH*timediff(tsbgn(i,:),tsend(i,:)));
            end
            if(any(steps==STEP.XC) && numel(vdata)==1); vdata=[]; end
            if(any(steps==STEP.XC) && numel(hdata)==1); hdata=[]; end
            if(any(steps==STEP.C3) && numel(vdata)<3); vdata=[]; end
            if(any(steps==STEP.C3) && numel(hdata)<6); hdata=[]; end
            if(any(steps==STEP.XCR) && numel(hdata)==1); hdata=[]; end
            if(isempty(hdata) && isempty(vdata)); continue; end
        end
        if(any(steps==3)) % remove trend
            if(~isempty(vdata)); vdata=removetrend(vdata); end
            if(~isempty(hdata)); hdata=removetrend(hdata); end
        end
        if(any(steps==4)) % taper
            if(~isempty(vdata))
                vdata=taper(vdata,opt.TAPERWIDTH,...
                    [],opt.TAPERTYPE,opt.TAPEROPT);
            end
            if(~isempty(hdata))
                hdata=taper(hdata,opt.TAPERWIDTH,...
                    [],opt.TAPERTYPE,opt.TAPEROPT);
            end
        end
        if(any(steps==5)) % resample
            if(~isempty(vdata)); vdata=syncrates(vdata,opt.SAMPLERATE); end
            if(~isempty(hdata)); hdata=syncrates(hdata,opt.SAMPLERATE); end
        end
        if(any(steps==6)) % remove pz
            if(~isempty(opt.PZDB))
                if(~isempty(vdata)); vdata=getsacpz(vdata,opt.PZDB); end
                if(~isempty(hdata)); hdata=getsacpz(hdata,opt.PZDB); end
            end
            if(~isempty(vdata))
                vdata=removesacpz(vdata,...
                    'units',opt.UNITS,'tl',opt.PZTAPERLIMITS);
            end
            if(~isempty(hdata))
                hdata=removesacpz(hdata,...
                    'units',opt.UNITS,'tl',opt.PZTAPERLIMITS);
            end
            if(any(steps==STEP.XC) && numel(vdata)==1); vdata=[]; end
            if(any(steps==STEP.XC) && numel(hdata)==1); hdata=[]; end
            if(any(steps==STEP.C3) && numel(vdata)<3); vdata=[]; end
            if(any(steps==STEP.C3) && numel(hdata)<6); hdata=[]; end
            if(any(steps==STEP.XCR) && numel(hdata)==1); hdata=[]; end
            if(isempty(hdata) && isempty(vdata)); continue; end
        end
        if(any(steps==STEP.AR)) % simple amplitude based rejection
            if(~isempty(vdata))
                [depmin,depmax]=getheader(vdata,'depmin','depmax');
                ampmax=max(abs(depmin),abs(depmax));
                switch lower(opt.REJECTMETHOD)
                    case 'abs'
                        vdata(ampmax>=opt.REJECTCUTOFF)=[];
                    case 'rms'
                        % grab pertinent rms values (the tough part)
                        kname=upper(getheader(vdata,'kname'));
                        kname=strcat(kname(:,1),'.',kname(:,2),'.',...
                            kname(:,3),'.',kname(:,4));
                        [tf,idx]=ismember(kname,rmsinfo.kname);
                        if(~all(tf))
                            error('seizmo:noise_process:badData',...
                                'Missing RMS info for some records!');
                        end
                        [tmp,tidx]=min(abs(...
                            rmsinfo.tsbgn-gregorian2serial(tsbgn(i,:))));
                        rms=getvaluefun(rmsinfo.data(idx),@(x)x(tidx));
                        vdata(ampmax>=opt.REJECTCUTOFF*rms)=[];
                end
            end
            if(~isempty(hdata))
                [depmin,depmax]=getheader(hdata,'depmin','depmax');
                ampmax=max(abs(depmin),abs(depmax));
                switch lower(opt.REJECTMETHOD)
                    case 'abs'
                        hdata(ampmax>=opt.REJECTCUTOFF)=[];
                    case 'rms'
                        % grab pertinent rms values (the tough part)
                        kname=upper(getheader(hdata,'kname'));
                        kname=strcat(kname(:,1),'.',kname(:,2),'.',...
                            kname(:,3),'.',kname(:,4));
                        [tf,idx]=ismember(kname,rmsinfo.kname);
                        if(~all(tf))
                            error('seizmo:noise_process:badData',...
                                'Missing RMS info for some records!');
                        end
                        [tmp,tidx]=min(abs(...
                            rmsinfo.tsbgn-gregorian2serial(tsbgn(i,:))));
                        rms=getvaluefun(rmsinfo.data(idx),@(x)x(tidx));
                        hdata(ampmax>=opt.REJECTCUTOFF*rms)=[];
                end
            end
            if(any(steps==STEP.XC) && numel(vdata)==1); vdata=[]; end
            if(any(steps==STEP.XC) && numel(hdata)==1); hdata=[]; end
            if(any(steps==STEP.C3) && numel(vdata)<3); vdata=[]; end
            if(any(steps==STEP.C3) && numel(hdata)<6); hdata=[]; end
            if(any(steps==STEP.XCR) && numel(hdata)==1); hdata=[]; end
            if(isempty(hdata) && isempty(vdata)); continue; end
        end
        if(any(steps==STEP.ROT)) % rotate horz to NE
            if(~isempty(hdata))
                hdata=rotate(hdata,'to',0,'kcmpnm1','N','kcmpnm2','E');
            end
            if(any(steps==STEP.XC) && numel(hdata)==1); hdata=[]; end
            if(any(steps==STEP.C3) && numel(hdata)<6); hdata=[]; end
            if(any(steps==STEP.XCR) && numel(hdata)==1); hdata=[]; end
            if(isempty(hdata) && isempty(vdata)); continue; end
        end
        if(any(steps==STEP.TDN)) % td norm
            % normalization style
            switch lower(opt.TDSTYLE)
                case 'flat'
                    % use robust rms (better for spikes & bad channels)
                    rms=getvaluefun([vdata; hdata],@(x)median(x.^2));
                    rms=sqrt(median(rms));
                    if(rms==0); rms=1; end
                    if(~isempty(vdata)); vdata=divide(vdata,rms); end
                    if(~isempty(hdata)); hdata=divide(hdata,rms); end
                case '1bit'
                    if(~isempty(vdata))
                        vdata=solofun(vdata,@sign);
                    end
                    if(~isempty(hdata))
                        % orthogonal pair 1bit: x^2+y^2=1
                        weights=solofun(addrecords(...
                            solofun(hdata(1:2:end),@(x)x.^2),...
                            solofun(hdata(2:2:end),@(x)x.^2)),@sqrt);
                        hdata=dividerecords(hdata,...
                            weights([1:end; 1:end]));
                    end
                case 'clip'
                    if(~isempty(vdata))
                        if(opt.TDRMS)
                            % use robust rms (better for spikes)
                            rms=getvaluefun(vdata,...
                                @(x)sqrt(median(x.^2)));
                            tdclip=rms*opt.TDCLIP;
                        else
                            tdclip=opt.TDCLIP;
                        end
                        vdata=clip(vdata,tdclip);
                    end
                    if(~isempty(hdata))
                        % orthogonal pair clipping
                        if(opt.TDRMS)
                            rms=getvaluefun(solofun(addrecords(...
                                solofun(hdata(1:2:end),@(x)x.^2),...
                                solofun(hdata(2:2:end),@(x)x.^2)),...
                                @sqrt),@(x)sqrt(median(x.^2)));
                            tdclip=rms([1:end; 1:end])*opt.TDCLIP;
                        else
                            tdclip=opt.TDCLIP;
                        end
                        hdata=clip(hdata,tdclip);
                    end
                case 'ram'
                    tdwb=opt.TDWEIGHTBAND;
                    if(~isempty(vdata))
                        delta=getheader(vdata(1),'delta');
                        %%% FOR CAMEROON VERIFICATION %%%
                        %weights=add(slidingabsmean(vdata,...
                        %    ceil(1/(2*delta*min(tdwb(:))))),eps);
                        %vdata=dividerecords(vdata,weights);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        weights=add(slidingabsmean(iirfilter(vdata,...
                            'bp','b','c',tdwb(1,:),'o',4,'p',2),...
                            ceil(1/(2*delta*min(tdwb(:))))),eps);
                        for a=2:size(tdwb,1)
                            weights=addrecords(weights,...
                                add(slidingabsmean(iirfilter(vdata,...
                                'bp','b','c',tdwb(a,:),'o',4,'p',2),...
                                ceil(1/(2*delta*min(tdwb(:))))),eps));
                        end
                        vdata=dividerecords(vdata,weights);
                    end
                    if(~isempty(hdata))
                        delta=getheader(hdata(1),'delta');
                        %%% FOR CAMEROON VERIFICATION %%%
                        %weights=add(slidingabsmean(hdata,...
                        %    ceil(1/(2*delta*min(tdwb(:))))),eps);
                        %weights=solofun(addrecords(...
                        %    solofun(weights(1:2:end),@(x)x.^2),...
                        %    solofun(weights(2:2:end),@(x)x.^2)),@sqrt);
                        %hdata=dividerecords(hdata,...
                        %    weights([1:end; 1:end]));
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        weights=add(slidingabsmean(iirfilter(hdata,...
                            'bp','b','c',tdwb(1,:),'o',4,'p',2),...
                            ceil(1/(2*delta*min(tdwb(:))))),eps);
                        for a=2:size(tdwb,1)
                            weights=addrecords(weights,...
                                add(slidingabsmean(iirfilter(hdata,...
                                'bp','b','c',tdwb(a,:),'o',4,'p',2),...
                                ceil(1/(2*delta*min(tdwb(:))))),eps));
                        end
                        weights=solofun(addrecords(...
                            solofun(weights(1:2:end),@(x)x.^2),...
                            solofun(weights(2:2:end),@(x)x.^2)),@sqrt);
                        hdata=dividerecords(hdata,...
                            weights([1:end; 1:end]));
                    end
                otherwise
                    error('seizmo:noise_process:badInput',...
                        'Unknown TDSTYLE: %s',opt.TDSTYLE);
            end
        end
        if(any(steps==STEP.FDN)) % fd norm
            % normalization style
            switch lower(opt.FDSTYLE)
                case '1bit'
                    if(~isempty(vdata))
                        vdata=dft(vdata);
                        vdata=solofun(vdata,@(x)[x(:,1).^0,x(:,2)]);
                        vdata=idft(vdata);
                        vdata=taper(vdata,opt.TAPERWIDTH,...
                            [],opt.TAPERTYPE,opt.TAPEROPT);
                    end
                    if(~isempty(hdata))
                        % orthogonal pair 1bit: x^2+y^2=1
                        hdata=dft(hdata,'rlim');
                        amph=rlim2amph(hdata);
                        amph=solofun(amph,...
                            @(x)x(:,[1:2:end; 1:2:end])+eps);
                        amph=solofun(addrecords(...
                            solofun(amph(1:2:end),@(x)x.^2),...
                            solofun(amph(2:2:end),@(x)x.^2)),@sqrt);
                        amph=changeheader(amph,'iftype','irlim');
                        hdata=dividerecords(hdata,...
                            amph([1:end; 1:end]));
                        hdata=idft(hdata);
                        hdata=taper(hdata,opt.TAPERWIDTH,...
                            [],opt.TAPERTYPE,opt.TAPEROPT);
                    end
                case 'ram'
                    if(~isempty(vdata))
                        %delta=getheader(vdata(1),'delta');         % cam
                        %nyqhz=1/(2*delta);                         % cam
                        %tnorm=[1/150/nyqhz (nyqhz-1/3)/nyqhz];     % cam
                        vdata=dft(vdata,'rlim');
                        vdata=whiten(vdata,opt.FDWIDTH);
                        %vdata=taper(vdata,tnorm,[],'gausswin',10); % cam
                        vdata=idft(vdata);
                        vdata=taper(vdata,opt.TAPERWIDTH,...
                            [],opt.TAPERTYPE,opt.TAPEROPT);
                    end
                    if(~isempty(hdata))
                        % orthogonal pair ram: x^2+y^2=1
                        %delta=getheader(hdata(1),'delta');         % cam
                        %nyqhz=1/(2*delta);                         % cam
                        %tnorm=[1/150/nyqhz (nyqhz-1/3)/nyqhz];     % cam
                        hdata=dft(hdata,'rlim');
                        amph=rlim2amph(hdata);
                        amph=slidingmean(amph,ceil(opt.FDWIDTH...
                            ./getheader(hdata,'delta')));
                        amph=solofun(amph,...
                            @(x)x(:,[1:2:end; 1:2:end])+eps);
                        amph=solofun(addrecords(...
                            solofun(amph(1:2:end),@(x)x.^2),...
                            solofun(amph(2:2:end),@(x)x.^2)),@sqrt);
                        amph=changeheader(amph,'iftype','irlim');
                        hdata=dividerecords(hdata,...
                            amph([1:end; 1:end]));
                        %hdata=taper(hdata,tnorm,[],'gausswin',10); % cam
                        hdata=idft(hdata);
                        hdata=taper(hdata,opt.TAPERWIDTH,...
                            [],opt.TAPERTYPE,opt.TAPEROPT);
                    end
                otherwise
                    error('seizmo:noise_process:badInput',...
                        'Unknown FDSTYLE: %s',opt.TDSTYLE);
            end
        end
        if(any(steps==STEP.INT)) % interpolation
            if(~isempty(vdata))
                [delta,f]=getheader(vdata(1),'delta','f');
                vdata=interpolate(vdata,1/delta,[],0,f,0);
            end
            if(~isempty(hdata))
                [delta,f]=getheader(hdata(1),'delta','f');
                hdata=interpolate(hdata,1/delta,[],0,f,0);
            end
        end
        if(any(steps==STEP.XC)) % xc
            if(numel(vdata)<2 && numel(hdata)<2); continue; end
            if(numel(vdata)>1)
                if(isempty(opt.FDOUT))
                    vdata=correlate(vdata,...
                        opt.NORMXC{:},opt.COHERENCY{:},'mcxc',...
                        opt.NOAUTO{:},opt.XCMAXLAG);
                else
                    vdata=correlate(vdata,...
                        opt.NORMXC{:},opt.COHERENCY{:},opt.FDOUT{:},...
                        opt.NOAUTO{:},'mcxc');
                end
                [vdata.path]=deal([indir fs tsdirs{i}(1:4) fs tsdirs{i}]);
            else
                vdata=vdata([]);
            end
            if(numel(hdata)>1)
                if(isempty(opt.FDOUT))
                    hdata=correlate(hdata,...
                        opt.NORMXC{:},opt.COHERENCY{:},'mcxc',...
                        opt.NOAUTO{:},opt.XCMAXLAG);
                else
                    hdata=correlate(hdata,...
                        opt.NORMXC{:},opt.COHERENCY{:},opt.FDOUT{:},...
                        opt.NOAUTO{:},'mcxc');
                end
                [hdata.path]=deal([indir fs tsdirs{i}(1:4) fs tsdirs{i}]);
            else
                hdata=hdata([]);
            end
        end
        if(any(steps==STEP.C3)) % correlate coda of correlations
            % correlation input or split data from before?
            if(isxc) % correlation input
                for c3=1:sum(steps==STEP.C3)
                    data=noise_c3(data,'vrayl',opt.C3VRAYL,...
                        'winoff',opt.C3WINOFF,'winlen',opt.C3WINLEN,...
                        'xccmp',opt.C3XCCMP,'ztrans',opt.C3ZTRANS,...
                        'mincoda',opt.C3MINCODA,'stnlist',opt.C3STNLIST,...
                        'xcmaxlag',opt.XCMAXLAG,...
                        'normxc',~isempty(opt.NORMXC),...
                        'coherency',~isempty(opt.COHERENCY),...
                        'fdout',~isempty(opt.FDOUT),...
                        'noauto',~isempty(opt.NOAUTO));
                    if(isempty(data)); break; end
                end
                if(isempty(data)); continue; end
                [data.path]=deal([indir fs tsdirs{i}(1:4) fs tsdirs{i}]);
            else % generated correlations
                for c3=1:sum(steps==STEP.C3)
                    if(~isempty(vdata))
                        vdata=noise_c3(vdata,'vrayl',opt.C3VRAYL,...
                            'winoff',opt.C3WINOFF,'winlen',opt.C3WINLEN,...
                            'xccmp',opt.C3XCCMP,'ztrans',opt.C3ZTRANS,...
                            'mincoda',opt.C3MINCODA,...
                            'stnlist',opt.C3STNLIST,...
                            'xcmaxlag',opt.XCMAXLAG,...
                            'normxc',~isempty(opt.NORMXC),...
                            'coherency',~isempty(opt.COHERENCY),...
                            'fdout',~isempty(opt.FDOUT),...
                            'noauto',~isempty(opt.NOAUTO));
                    end
                    if(~isempty(hdata))
                        hdata=noise_c3(hdata,'vrayl',opt.C3VRAYL,...
                            'winoff',opt.C3WINOFF,'winlen',opt.C3WINLEN,...
                            'xccmp',opt.C3XCCMP,'ztrans',opt.C3ZTRANS,...
                            'mincoda',opt.C3MINCODA,...
                            'stnlist',opt.C3STNLIST,...
                            'xcmaxlag',opt.XCMAXLAG,...
                            'normxc',~isempty(opt.NORMXC),...
                            'coherency',~isempty(opt.COHERENCY),...
                            'fdout',~isempty(opt.FDOUT),...
                            'noauto',~isempty(opt.NOAUTO));
                    end
                end
                if(isempty(hdata) && isempty(vdata)); continue; end
                [vdata.path]=deal([indir fs tsdirs{i}(1:4) fs tsdirs{i}]);
                [hdata.path]=deal([indir fs tsdirs{i}(1:4) fs tsdirs{i}]);
            end
        end
        if(any(steps==STEP.XCR)) % rotate xc
            % correlation input or split data from before?
            if(isxc) % correlation input
                % this removes ZZ correlations!
                data=rotate_correlations(data,'rt');
                if(isempty(data)); continue; end
            else % generated correlations
                if(~isempty(hdata))
                    hdata=rotate_correlations(hdata,'rt');
                    if(isempty(hdata) && isempty(vdata)); continue; end
                end
            end
        end
        
        % write the data
        if(opt.MATIO)
            if(isxc)
                if(isempty(data)); continue; end
                noise_records=changepath(data,'path',...
                    [outdir fs tsdirs{i}(1:4) fs tsdirs{i} fs]);
            else
                if(isempty(hdata) && isempty(vdata)); continue; end
                noise_records=changepath([vdata; hdata],'path',...
                    [outdir fs tsdirs{i}(1:4) fs tsdirs{i} fs]);
            end
            if(~exist([outdir fs tsdirs{i}(1:4) fs tsdirs{i}],'dir'))
                mkdir([outdir fs tsdirs{i}(1:4) fs tsdirs{i}]);
            end
            save([outdir fs tsdirs{i}(1:4) fs tsdirs{i} fs ...
                'noise_records.mat'],'noise_records');
            clear noise_records;
        else
            if(isxc)
                if(isempty(data)); continue; end
                writeseizmo(data,'path',...
                    [outdir fs tsdirs{i}(1:4) fs tsdirs{i} fs]);
            else
                if(isempty(hdata) && isempty(vdata)); continue; end
                writeseizmo([vdata; hdata],'path',...
                        [outdir fs tsdirs{i}(1:4) fs tsdirs{i} fs]);
            end
        end
        
        % toggle checking back
        seizmocheck_state(oldseizmocheckstate);
        checkheader_state(oldcheckheaderstate);
    catch
        % toggle checking back
        seizmocheck_state(oldseizmocheckstate);
        checkheader_state(oldcheckheaderstate);
        
        % parallel processing takedown & fix verbosity
        %matlabpool close; % PARALLEL
        seizmoverbose(verbose);
        
        % rethrow error
        error(lasterror);
    end
end

% parallel processing takedown & fix verbosity
%matlabpool close; % PARALLEL
seizmoverbose(verbose);

end


function [opt]=noise_process_parameters(varargin)
% parses/checks noise_process pv pairs

% defaults
varargin=[{'minlen' 70 'tw' 1 'tt' [] 'topt' [] 'sr' 1 'rw' 0 ...
    'pzdb' [] 'units' 'disp' 'pztl' [.004 .008] 'rc' Inf 'rm' 'abs' ...
    'tds' 'flat' 'tdrms' true 'tdclip' 1 'tdfb' [1/150 1/3; 1/100 1/15] ...
    'fds' 'ram' 'fdw' .002 'lag' 4000 'nxc' true 'coh' false 'c3mc' 1 ...
    'c3v' 2.6 'c3wo' 0 'c3wl' 4e3 'c3x' 'sym' 'c3z' true 'c3s' {} ...
    'ts' [] 'te' [] 'lat' [] 'lon' [] 'fdout' false 'noauto' false ...
    'net' [] 'sta' [] 'str' [] 'cmp' [] 'file' [] 'q' false} varargin];

% get user input
for i=1:2:numel(varargin)
    switch lower(varargin{i})
        case {'minimumlength' 'ml' 'minlen' 'mlen' 'minl' 'minlength'}
            if(isempty(varargin{i+1})); continue; end
            opt.MINIMUMLENGTH=varargin{i+1};
        case {'taperwidth' 'tw' 'taperw' 'tapw' 'twidth'}
            if(isempty(varargin{i+1})); continue; end
            opt.TAPERWIDTH=varargin{i+1};
        case {'tapertype' 'tt' 'tapert' 'tapt' 'ttype'}
            opt.TAPERTYPE=varargin{i+1};
        case {'taperopt' 'topt' 'tapero' 'tapo' 'tapopt' 'taperoption'}
            opt.TAPEROPT=varargin{i+1};
        case {'samplerate' 'sr' 'srate'}
            if(isempty(varargin{i+1})); continue; end
            opt.SAMPLERATE=varargin{i+1};
        case {'pzdb' 'db' 'pz'}
            opt.PZDB=varargin{i+1};
        case {'units' 'to' 'u' 'unit' 'un'}
            if(isempty(varargin{i+1})); continue; end
            opt.UNITS=varargin{i+1};
        case {'pztaperlimits' 'pztl' 'tl' 'taperlim' 'tlim' 'taplim'}
            if(isempty(varargin{i+1})); continue; end
            opt.PZTAPERLIMITS=varargin{i+1};
        case {'rejectcutoff' 'cutoff' 'rejcut' 'cut' 'rc' 'rco'}
            if(isempty(varargin{i+1})); continue; end
            opt.REJECTCUTOFF=varargin{i+1};
        case {'rejectwidth' 'width' 'rejwid' 'wid' 'rw'}
            if(isempty(varargin{i+1})); continue; end
            opt.REJECTWIDTH=varargin{i+1};
        case {'rejectmethod' 'method' 'rejmeth' 'rm' 'meth'}
            if(isempty(varargin{i+1})); continue; end
            opt.REJECTMETHOD=varargin{i+1};
        case {'tdstyle' 'td' 'tds'}
            if(isempty(varargin{i+1})); continue; end
            opt.TDSTYLE=varargin{i+1};
        case {'tdrms' 'rms'}
            if(isempty(varargin{i+1})); continue; end
            opt.TDRMS=varargin{i+1};
        case {'tdclip' 'clip'}
            if(isempty(varargin{i+1})); continue; end
            opt.TDCLIP=varargin{i+1};
        case {'tdweightband' 'tdfreqband' 'tdwb' 'tdfb'}
            if(isempty(varargin{i+1})); continue; end
            opt.TDWEIGHTBAND=varargin{i+1};
        case {'fdstyle' 'fd' 'fds'}
            if(isempty(varargin{i+1})); continue; end
            opt.FDSTYLE=varargin{i+1};
        case {'fdwidth' 'fdw'}
            if(isempty(varargin{i+1})); continue; end
            opt.FDWIDTH=varargin{i+1};
        case {'xcmaxlag' 'xcmax' 'xclag' 'maxlag' 'lag'}
            if(isempty(varargin{i+1})); continue; end
            opt.XCMAXLAG=varargin{i+1};
        case {'normxc' 'nxc'}
            if(isempty(varargin{i+1})); continue; end
            opt.NORMXC=varargin{i+1};
        case {'coherency' 'cohere' 'co' 'coh'}
            if(isempty(varargin{i+1})); continue; end
            opt.COHERENCY=varargin{i+1};
        case {'fdout' 'fdo' 'fo'}
            if(isempty(varargin{i+1})); continue; end
            opt.FDOUT=varargin{i+1};
        case {'noauto' 'noa'}
            if(isempty(varargin{i+1})); continue; end
            opt.NOAUTO=varargin{i+1};
        case {'c3v' 'c3vr' 'c3vrayl'}
            if(isempty(varargin{i+1})); continue; end
            opt.C3VRAYL=varargin{i+1};
        case {'c3wo' 'c3wino' 'c3winoff'}
            if(isempty(varargin{i+1})); continue; end
            opt.C3WINOFF=varargin{i+1};
        case {'c3wl' 'c3winl' 'c3winlen'}
            if(isempty(varargin{i+1})); continue; end
            opt.C3WINLEN=varargin{i+1};
        case {'c3x' 'c3xc' 'c3xccmp'}
            if(isempty(varargin{i+1})); continue; end
            opt.C3XCCMP=varargin{i+1};
        case {'c3z' 'c3zt' 'c3ztrans' 'c3ztransform'}
            if(isempty(varargin{i+1})); continue; end
            opt.C3ZTRANS=varargin{i+1};
        case {'c3s' 'c3sl' 'c3stn' 'c3stnlist'}
            opt.C3STNLIST=varargin{i+1};
        case {'c3m' 'c3mc' 'c3min' 'c3minc' 'c3mincoda'}
            if(isempty(varargin{i+1})); continue; end
            opt.C3MINCODA=varargin{i+1};
        case {'ts' 'tstart' 'timestart' 'startt' 'stime' 'starttime'}
            opt.TIMESTART=varargin{i+1};
        case {'te' 'tend' 'timeend' 'et' 'endtime' 'etime' 'endt'}
            opt.TIMEEND=varargin{i+1};
        case {'lat' 'la' 'lar' 'latr' 'larng' 'latitude' 'latrng'}
            opt.LATRNG=varargin{i+1};
        case {'lon' 'lo' 'lor' 'lonr' 'lorng' 'longitude' 'lonrng'}
            opt.LONRNG=varargin{i+1};
        case {'knetwk' 'n' 'net' 'netwk' 'network' 'nets' 'networks'}
            opt.NETWORKS=varargin{i+1};
        case {'kstnm' 'st' 'sta' 'stn' 'stns' 'stations' 'station'}
            if(~isempty(varargin{i+1}) && isnumeric(varargin{i+1}))
                % timestart/starttime catch
                warning('seizmo:noise_process:badInput',...
                    'TIMESTART/STATION mixup!  Assuming TIMESTART!');
                opt.TIMESTART=varargin{i+1};
            else
                opt.STATIONS=varargin{i+1};
            end
        case {'khole' 'hole' 'holes' 'str' 'strs' 'stream' 'streams'}
            opt.STREAMS=varargin{i+1};
        case {'kcmpnm' 'cmpnm' 'cmp' 'cmps' 'component' 'components'}
            opt.COMPONENTS=varargin{i+1};
        case {'f' 'file' 'filename' 'files' 'filenames'}
            opt.FILENAMES=varargin{i+1};
        case {'q' 'qw' 'quiet' 'qwrite' 'quietwrite'}
            if(isempty(varargin{i+1})); continue; end
            opt.QUIETWRITE=varargin{i+1};
        otherwise
            error('seizmo:noise_process:badInput',...
                'Unknown Option: %s !',varargin{i});
    end
end

% fix string options to be cellstr vectors
if(ischar(opt.C3STNLIST)); opt.C3STNLIST=cellstr(opt.C3STNLIST); end
if(ischar(opt.NETWORKS)); opt.NETWORKS=cellstr(opt.NETWORKS); end
if(ischar(opt.STATIONS)); opt.STATIONS=cellstr(opt.STATIONS); end
if(ischar(opt.STREAMS)); opt.STREAMS=cellstr(opt.STREAMS); end
if(ischar(opt.COMPONENTS)); opt.COMPONENTS=cellstr(opt.COMPONENTS); end
if(ischar(opt.FILENAMES)); opt.FILENAMES=cellstr(opt.FILENAMES); end
if(iscellstr(opt.C3STNLIST))
    opt.C3STNLIST=unique(lower(opt.C3STNLIST(:)));
end
if(iscellstr(opt.NETWORKS))
    opt.NETWORKS=unique(lower(opt.NETWORKS(:)));
end
if(iscellstr(opt.STATIONS))
    opt.STATIONS=unique(lower(opt.STATIONS(:)));
end
if(iscellstr(opt.STREAMS))
    opt.STREAMS=unique(lower(opt.STREAMS(:)));
end
if(iscellstr(opt.COMPONENTS))
    opt.COMPONENTS=unique(lower(opt.COMPONENTS(:)));
end
if(iscellstr(opt.FILENAMES))
    opt.FILENAMES=unique(opt.FILENAMES(:));
end

% check options
szs=size(opt.TIMESTART);
sze=size(opt.TIMEEND);
szp=size(opt.PZTAPERLIMITS);
if(~isscalar(opt.MINIMUMLENGTH) || ~isreal(opt.MINIMUMLENGTH) ...
        || opt.MINIMUMLENGTH<0 || opt.MINIMUMLENGTH>100)
    error('seizmo:noise_process:badInput',...
        'MINIMUMLENGTH must be a scalar within 0 & 100 (%%)!');
elseif(~isscalar(opt.TAPERWIDTH) || ~isreal(opt.TAPERWIDTH) ...
        || opt.TAPERWIDTH<0 || opt.TAPERWIDTH>100)
    error('seizmo:noise_process:badInput',...
        'TAPERWIDTH must be a scalar within 0 & 100 (%%)!');
elseif(~isempty(opt.TAPERTYPE) && (~ischar(opt.TAPERTYPE) ...
        || ~isvector(opt.TAPERTYPE) || size(opt.TAPERTYPE,1)~=1))
    error('seizmo:noise_process:badInput',...
        'TAPERTYPE must be a string indicating a valid taper!');
elseif(~isempty(opt.TAPEROPT) && (~isscalar(opt.TAPEROPT) ...
        || ~isreal(opt.TAPEROPT)))
    error('seizmo:noise_process:badInput',...
        'TAPEROPT must be a real-valued scalar!');
elseif(~isscalar(opt.SAMPLERATE) || ~isreal(opt.SAMPLERATE) ...
        || opt.SAMPLERATE<=0)
    error('seizmo:noise_process:badInput',...
        'SAMPLERATE must be a positive real-valued scalar!');
elseif(~isempty(opt.PZDB) && ~isstruct(opt.PZDB) && (~ischar(opt.PZDB) ...
        || ~isvector(opt.PZDB) || size(opt.PZDB,1)~=1))
    error('seizmo:noise_process:badInput',...
        'PZDB must be either a struct or a string!');
elseif(~ischar(opt.UNITS) || ~isvector(opt.UNITS) || size(opt.UNITS,1)~=1)
    error('seizmo:noise_process:badInput',...
        'UNITS must be a string!');
elseif(~isempty(opt.PZTAPERLIMITS) && (numel(szp)>2 || szp(1)~=1 ...
        || all(szp(2)~=[2 4]) || ~isnumeric(opt.PZTAPERLIMITS) ...
        || ~isreal(opt.PZTAPERLIMITS) || any(opt.PZTAPERLIMITS<0)))
    error('seizmo:noise_process:badInput',...
        'PZTAPERLIMITS must be a [LOWSTOP LOWPASS] or [LS LP HP HS]!');
elseif(~ischar(opt.REJECTMETHOD) || ~isvector(opt.REJECTMETHOD) ...
        || size(opt.REJECTMETHOD,1)~=1 ...
        || ~ismember(lower(opt.REJECTMETHOD),{'abs' 'rms'}))
    error('seizmo:noise_process:badInput',...
        'REJECTMETHOD must be either ''ABS'' or ''RMS''!');
elseif(~isscalar(opt.REJECTWIDTH) || ~isreal(opt.REJECTWIDTH) ...
        || fix(opt.REJECTWIDTH)~=opt.REJECTWIDTH || opt.REJECTWIDTH<0)
    error('seizmo:noise_process:badInput',...
        'REJECTWIDTH must be a positive integer!');
elseif(~isscalar(opt.REJECTCUTOFF) || ~isreal(opt.REJECTCUTOFF) ...
        || opt.REJECTCUTOFF<=0)
    error('seizmo:noise_process:badInput',...
        'REJECTCUTOFF must be a positive scalar!');
elseif(~ischar(opt.TDSTYLE) || ~isvector(opt.TDSTYLE) ...
        || size(opt.TDSTYLE,1)~=1)
    error('seizmo:noise_process:badInput',...
        'TDSTYLE must be a string!');
elseif(~isscalar(opt.TDRMS) || ~islogical(opt.TDRMS))
    error('seizmo:noise_process:badInput',...
        'TDRMS must be true or false!');
elseif(~isscalar(opt.TDCLIP) || ~isreal(opt.TDCLIP))
    error('seizmo:noise_process:badInput',...
        'TDCLIP must be a real-valued scalar!');
elseif(~isnumeric(opt.TDWEIGHTBAND) || ~isreal(opt.TDWEIGHTBAND) ...
        || size(opt.TDWEIGHTBAND,2)~=2 ...
        || numel(size(opt.TDWEIGHTBAND))~=2 || any(opt.TDWEIGHTBAND(:)<0))
    error('seizmo:noise_process:badInput',...
        'TDWEIGHTBAND must be [LOW HIGH] in Hz!');
elseif(~ischar(opt.FDSTYLE) || ~isvector(opt.FDSTYLE) ...
        || size(opt.FDSTYLE,1)~=1)
    error('seizmo:noise_process:badInput',...
        'FDSTYLE must be a string!');
elseif(~isscalar(opt.FDWIDTH) || ~isreal(opt.FDWIDTH) || opt.FDWIDTH<=0)
    error('seizmo:noise_process:badInput',...
        'FDWIDTH must be a positive real-valued scalar!');
elseif(~isscalar(opt.XCMAXLAG) || ~isreal(opt.XCMAXLAG) ...
        || opt.XCMAXLAG<=0)
    error('seizmo:noise_process:badInput',...
        'XCMAXLAG must be a positive real-valued scalar in seconds!');
elseif(~isscalar(opt.NORMXC) || ~islogical(opt.NORMXC))
    error('seizmo:noise_process:badInput',...
        'NORMXC must be true or false!');
elseif(~isscalar(opt.COHERENCY) || ~islogical(opt.COHERENCY))
    error('seizmo:noise_process:badInput',...
        'COHERENCY must be true or false!');
elseif(~isscalar(opt.FDOUT) || ~islogical(opt.FDOUT))
    error('seizmo:noise_process:badInput',...
        'FDOUT must be true or false!');
elseif(~isscalar(opt.NOAUTO) || ~islogical(opt.NOAUTO))
    error('seizmo:noise_process:badInput',...
        'NOAUTO must be true or false!');
elseif(~isnumeric(opt.C3VRAYL) || ~isscalar(opt.C3VRAYL) ...
        || ~isreal(opt.C3VRAYL) || opt.C3VRAYL<=0)
    error('seizmo:noise_process:badInput',...
        'C3VRAYL must be a positive real scalar!');
elseif(~isnumeric(opt.C3WINOFF) || ~isscalar(opt.C3WINOFF) ...
        || ~isreal(opt.C3WINOFF))
    error('seizmo:noise_process:badInput',...
        'C3WINOFF must be a real-valued scalar!');
elseif(~isnumeric(opt.C3WINLEN) || ~isscalar(opt.C3WINLEN) ...
        || ~isreal(opt.C3WINLEN) || opt.C3WINLEN<=0)
    error('seizmo:noise_process:badInput',...
        'C3WINLEN must be a positive real scalar!');
elseif(~ischar(opt.C3XCCMP) || ~isvector(opt.C3XCCMP) ...
        || size(opt.C3XCCMP,1)~=1 || ~ismember(lower(opt.C3XCCMP),...
        {'causal' 'acausal' 'both' 'sym'}))
    error('seizmo:noise_process:badInput',...
        'C3XCCMP must be ''CAUSAL'' ''ACAUSAL'' ''BOTH'' or ''SYM''!');
elseif(~isscalar(opt.C3ZTRANS) || ~islogical(opt.C3ZTRANS))
    error('seizmo:noise_process:badInput',...
        'C3ZTRANS must be true or false!');
elseif(~isempty(opt.C3STNLIST) && (~iscellstr(opt.C3STNLIST)))
    error('seizmo:noise_process:badInput',...
        'C3STNLIST must be a string list of KNAME codes!');
elseif(~isnumeric(opt.C3MINCODA) || ~isscalar(opt.C3MINCODA) ...
        || ~isreal(opt.C3MINCODA) || opt.C3MINCODA<=0 ...
        || opt.C3MINCODA~=fix(opt.C3MINCODA))
    error('seizmo:noise_process:badInput',...
        'C3MINCODA must be a positive integer>0!');
elseif(~isscalar(opt.QUIETWRITE) || ~islogical(opt.QUIETWRITE))
    error('seizmo:noise_process:badInput',...
        'QUIETWRITE flag must be a scalar logical!');
elseif(~isempty(opt.TIMESTART) && (numel(szs)>2 || szs(1)~=1 ...
        || all(szs(2)~=[2 3 5 6]) || ~isnumeric(opt.TIMESTART) ...
        || ~isreal(opt.TIMESTART)))
    error('seizmo:noise_process:badInput',...
        'TIMESTART must be a recognized date-time vector!');
elseif(~isempty(opt.TIMEEND) && (numel(sze)>2 || sze(1)~=1 ...
        || all(sze(2)~=[2 3 5 6]) || ~isnumeric(opt.TIMEEND) ...
        || ~isreal(opt.TIMEEND)))
    error('seizmo:noise_process:badInput',...
        'TIMEEND must be a recognized date-time vector!');
elseif(~isempty(opt.LATRNG) && (~isnumeric(opt.LATRNG) ...
        || ~isreal(opt.LATRNG) || numel(opt.LATRNG)~=2 ...
        || size(opt.LATRNG,2)~=2 || numel(size(opt.LATRNG))~=2))
    error('seizmo:noise_process:badInput',...
        'LATRNG must be a 2 element numeric vector as [LOW HIGH]!');
elseif(~isempty(opt.LONRNG) && (~isnumeric(opt.LONRNG) ...
        || ~isreal(opt.LONRNG) || numel(opt.LONRNG)~=2 ...
        || size(opt.LONRNG,2)~=2 || numel(size(opt.LONRNG))~=2))
    error('seizmo:noise_process:badInput',...
        'LONRNG must be a 2 element numeric vector as [LOW HIGH]!');
elseif(~isempty(opt.NETWORKS) && (~iscellstr(opt.NETWORKS)))
    error('seizmo:noise_process:badInput',...
        'NETWORKS must be a string list of allowed network codes!');
elseif(~isempty(opt.STATIONS) && (~iscellstr(opt.STATIONS)))
    error('seizmo:noise_process:badInput',...
        'STATIONS must be a string list of allowed station codes!');
elseif(~isempty(opt.STREAMS) && (~iscellstr(opt.STREAMS)))
    error('seizmo:noise_process:badInput',...
        'STREAMS must be a string list of allowed stream codes!');
elseif(~isempty(opt.COMPONENTS) && (~iscellstr(opt.COMPONENTS)))
    error('seizmo:noise_process:badInput',...
        'COMPONENTS must be a string list of allowed component codes!');
elseif(~isempty(opt.FILENAMES) && (~iscellstr(opt.FILENAMES)))
    error('seizmo:noise_process:badInput',...
        'FILENAMES must be a string list of allowed files!');
end

% percent to fraction
opt.MINIMUMLENGTH=opt.MINIMUMLENGTH/100;
opt.TAPERWIDTH=opt.TAPERWIDTH/100;

% xc logical to option
if(opt.NORMXC); opt.NORMXC={'normxc'}; else opt.NORMXC={}; end
if(opt.COHERENCY); opt.COHERENCY={'coherency'}; else opt.COHERENCY={}; end
if(opt.FDOUT); opt.FDOUT={'fdout'}; else opt.FDOUT={}; end
if(opt.NOAUTO); opt.NOAUTO={'noauto'}; else opt.NOAUTO={}; end

end


function [d1]=dividerecords(d1,d2,varargin)
% simple hack for speed (also updates dep* fields)
nrecs=numel(d1);
[depmin,depmax,depmen]=deal(nan(nrecs,1));
for i=1:nrecs
    d1(i).dep=d1(i).dep./d2(i).dep;
    depmin(i)=min(d1(i).dep(:));
    depmax(i)=max(d1(i).dep(:));
    depmen(i)=nanmean(d1(i).dep(:));
end
d1=changeheader(d1,'depmin',depmin,'depmax',depmax,'depmen',depmen);
end


function [d1]=addrecords(d1,d2,varargin)
% simple hack for speed (no dep* update)
for i=1:numel(d1); d1(i).dep=d1(i).dep+d2(i).dep; end
end

