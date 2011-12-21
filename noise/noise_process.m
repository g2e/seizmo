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
%     written to OUTDIR.  The following techniques are done on the input
%     dataset:
%      ( 1) remove flat records (no change in recorded value)
%      ( 2) remove short records (any less that 70% of time section)
%      ( 3) remove mean & trend
%      ( 4) taper (first/last 1%)
%      ( 5) resample to 1 sample/sec
%      ( 6) remove polezero response (to: disp, taper: [1/250 1/150])
%      ( 7) NOT IMPLEMENTED
%      ( 8) rotate horizontals to North/East (also removes unpaired, cuts)
%      ( 9) t-domain normalize (2pass moving average: unfilt,15-100s)
%      (10) f-domain normalize (2mHz moving average)
%      (11) correlate (keep +/-4000s lagtime)
%      (12) rotate correlations
%     See below for details on how to alter or skip some of these steps.
%
%     NOISE_PROCESS(INDIR,OUTDIR,STEPS) only does the processing steps
%     indicated in STEPS.  STEPS should be a vector of numbers
%     corresponding to valid steps given above.  The default is [] which
%     does all of the above steps.
%
%     NOISE_PROCESS(INDIR,OUTDIR,STEPS,'OPT1',VAL,...,'OPTN',VAL) allows
%     changing some of the noise correlation parameters.  The following
%     options are configurable:
%      MINIMUMLENGTH - minimum length of records in % of time section [70]
%      TAPERWIDTH - [0.01]
%      TAPERTYPE - []
%      TAPEROPTION - []
%      SAMPLERATE - [1]
%      PZDB - polezero db to use: []
%      UNITS - ['disp']
%      PZTAPERLIMITS - [1/250 1/150]
%      MWTHRESH - minimum CMT mw for reject [inf]
%      TDSTYLE - '1bit'
%                'rmsclip' - clip values above ?xRMS
%                'rmsreject' - reject values (+nearby) above ?xRMS
%                'ram' - normalized using running-absolute mean
%                ['1bit+ram'] - 1bit & running-absolute mean
%      TDRMSDETECT - ? x RMS [1]
%      TDRMSREJECT - minutes to reject around point [30]
%      TDFREQBAND - freq-band for running-abs mean [1/100 1/15]
%      FDSTYLE - '1bit' - set all amplitudes to 1 (horizontals are special)
%                'ram' - normalized using running absolute mean
%      FDWIDTH - width of window for running-abs mean [2mHz]
%      FDPASSBAND - taper outside this band [1/100 1/3]
%      XCMAXLAG - maximum lag time of output correlograms in sec [4000]
%      TIMESTART - process time sections from this time on []
%      TIMEEND - process times sections before this time []
%      LATRNG - include stations in this latitude range []
%      LONRNG - include stations in this longitude range []
%      NETWORKS - include records from these networks []
%      STATIONS - include records from these stations []
%      STREAMS - include records with these stream codes []
%      COMPONENTS - include records with these components []
%      FILENAMES - limit processing to files with these filenames []
%
%    Notes:
%     - Good Noise Analysis References:
%        Bensen et al 2007, GJI, doi:10.1111/j.1365-246X.2007.03374.x
%        Yang et al 2007, GJI, doi:10.1111/j.1365-246X.2006.03203.x
%        Lin et al 2008, GJI, doi:10.1111/j.1365-246X.2008.03720.x
%        Harmon et al 2008, GRL, doi:10.1029/2008GL035387
%        Prieto et al 2009, JGR, doi:10.1029/2008JB006067
%        Ekstrom et al 2009, GRL, doi:10.1029/2009GL039131
%
%    Examples:
%     % Perform the first 3 steps of noise processing,
%     % writing out the resulting data of each step:
%     noise_process('raw','step1',1)
%     noise_process('step1','step2',2)
%     noise_process('step2','step3',3)
%
%     % Skip the normalization step:
%     noise_process('raw','xc',[1:8 11:12])
%
%    See also: NOISE_SETUP, NOISE_STACK, NOISE_WORKFLOW

%     Version History:
%        Nov. 22, 2011 - initial version (only first 6 steps)
%        Nov. 29, 2011 - more steps, subsetting
%        Dec. 13, 2011 - fix output writing, correlogram subsetting
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec. 13, 2011 at 11:15 GMT

% todo:
% - normalization steps
% - correlogram input needs to respect subsetting

% check nargin
error(nargchk(2,inf,nargin));
if(nargin>=4 && ~mod(nargin,2))
    error('seizmo:noise_process:badInput',...
        'Unpaired option/value pair given!');
end

% default steps to [] (all)
if(nargin<3); steps=[]; end

% parse/check options
opt=noise_process_parameters(varargin{:});

% check directories
if(~ischar(indir) || ~isvector(indir))
    error('seizmo:noise_process:fileNotString',...
        'INDIR must be a string!');
end
if(~exist(indir,'dir'))
    error('seizmo:noise_process:dirConflict',...
        ['Input Directory: %s\n' ...
        'Does not exist (or is not a directory)!'],indir);
end
if(~ischar(outdir) || ~isvector(outdir))
    error('seizmo:noise_process:fileNotString',...
        'OUTDIR must be a string!');
end
if(exist(outdir,'file'))
    if(~exist(outdir,'dir'))
        error('seizmo:noise_process:dirConflict',...
            'Output Directory: %s\nIs a file!',outdir);
    end
    if(~qw)
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

% parallel processing setup (up to 8 instances)
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

% verbosity  (turn it off for the loop)
verbose=seizmoverbose(false);
if(verbose); disp('Processing seismic data for noise analysis'); end

% loop over years
for i=1:nyr
    % skip if outside user-defined time range
    if((~isempty(opt.TIMESTART) ...
            && str2double(yrdir{i})<opt.TIMESTART(1)) ...
            || (~isempty(opt.TIMEEND) ...
            && str2double(yrdir{i})>opt.TIMEEND(1)))
        continue;
    end
    
    % loop over time section directories
    for j=1:numel(tsdir{i})
    %parfor j=1:numel(tsdir{i})
        % read the data
        try
            data=readseizmo(...
                [indir fs yrdir{i} fs tsdir{i}{j} fs opt.FILENAMES]);
            hvsplit=false;
        catch
            % no data...
            continue;
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
        if(isxc && any(steps<12))
            error('seizmo:noise_process:invalidProcess4xcdata',...
                'Cannot run earlier processing steps on correlograms!');
        end
        
        % proceed by data type
        if(isxc) % correlograms
            %%% REQUIRES ATTENTION!!!!! %%%
            % this could get complicated...
        else % seismic records
            % find/check time section limits
            [tsbgn,tsend]=getheader(data,'a utc','f utc');
            tsbgn=unique(tsbgn,'rows');
            tsend=unique(tsend,'rows');
            if(size(tsbgn,1)>1 || size(tsend,1)>1)
                error('seizmo:noise_process:inconsistentSetup',...
                    'Time window limits are inconsistent!');
            end
            
            % skip if outside user-defined time range
            if((~isempty(opt.TIMESTART) ...
                    && timediff(opt.TIMESTART,tsend)<=0) ...
                    || (~isempty(opt.TIMEEND) ...
                    && timediff(opt.TIMEEND,tsbgn)>=0))
                continue;
            end
            
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
                data=data(ismember(knetwk,opt.NETWORKS));
                if(isempty(data)); continue; end
            end
            if(~isempty(opt.STATIONS))
                kstnm=getheader(data,'kstnm');
                data=data(ismember(kstnm,opt.STATIONS));
                if(isempty(data)); continue; end
            end
            if(~isempty(opt.STREAMS))
                khole=getheader(data,'khole');
                data=data(ismember(khole,opt.STREAMS));
                if(isempty(data)); continue; end
            end
            if(~isempty(opt.COMPONENTS))
                kcmpnm=getheader(data,'kcmpnm');
                data=data(ismember(kcmpnm,opt.COMPONENTS));
                if(isempty(data)); continue; end
            end
        end
        
        % detail message
        if(verbose); disp(['PROCESSING: ' tsdir{i}{j}]); end
        
        % process data for noise analysis
        if(any(steps==1)) % remove dead
            data=removedeadrecords(data);
            if(isempty(data)); continue; end
        end
        if(any(steps==2)) % remove short
            [b,e]=getheader(data,'b','e');
            data=data(e-b>opt.MINIMUMLENGTH*timediff(tsbgn,tsend));
            if(isempty(data)); continue; end
        end
        if(any(steps==3)) % remove trend
            data=removetrend(data);
        end
        if(any(steps==4)) % taper
            data=taper(data,opt.TAPERWIDTH,[],opt.TAPERTYPE,opt.TAPEROPT);
        end
        if(any(steps==5)) % resample
            data=syncrates(data,opt.SAMPLERATE);
        end
        if(any(steps==6)) % remove pz
            if(~isempty(opt.PZDB)); data=getsacpz(data,opt.PZDB); end
            data=removesacpz(data,...
                'units',opt.UNITS,'tl',opt.PZTAPERLIMITS);
            if(isempty(data)); continue; end
        end
        if(any(steps==7)) % cmt reject?
            % PLACEHOLDER -- CURRENTLY UNIMPLEMENTED
            % Search just time section or some time before too?
            % How do we define that time? maximum...
            % Find events meeting mw & distance criteria?
            % Rhie & Romanowicz 2006 is a good paper for this method.
            % This seems better figured out in an a posteriori way...
        end
        
        % HIDDEN STEP required for 8+
        % split data into vertical and horizontal sets
        if(~isxc && any(steps>7))
            vdata=data(vertcmp(data));
            hdata=data(horzcmp(data));
            hvsplit=true;
            clear data;
            if(isempty(hdata) && isempty(vdata)); continue; end
        end
        
        % continue processing data for noise analysis
        if(any(steps==8)) % rotate horz to NE
            hdata=rotate(hdata,'to',0,'kcmpnm1','N','kcmpnm2','E');
        end
        if(any(steps==9)) % td norm
            % dont forget EE/EN/NE/NN
            % we have several ways
            % - 1bit
            % - rms clip
            % - rms reject
            % - running abs mean
            % - 1bit + running abs mean (degault)
            switch opt.TDSTYLE
                case '1bit'
                    if(~isempty(vdata))
                        
                    end
                    if(~isempty(hdata))
                        
                    end
                case 'rmsclip'
                case 'rmsreject'
                case 'ram'
                case '1bit+ram'
            end
        end
        if(any(steps==10)) % fd norm
            % dont forget EE/EN/NE/NN
            
        end
        if(any(steps==11)) % xc
            % need to separate ZZ and EE/EN/NE/NN
            if(numel(vdata)<2 && numel(hdata)<2); continue; end
            if(numel(vdata)>1)
                delta=getheader(vdata(1),'delta');
                vdata=interpolate(correlate(...
                    cut(vdata,'a','f','fill',true),...
                    'lags',(opt.XCMAXLAG+4*delta).*[-1 1]),...
                    1/delta,[],-opt.XCMAXLAG,opt.XCMAXLAG);
                [vdata.path]=deal([indir fs yrdir{i} fs tsdir{i}{j}]);
            else
                vdata=vdata([]);
            end
            if(numel(hdata)>1)
                delta=getheader(hdata(1),'delta');
                hdata=interpolate(correlate(...
                    cut(hdata,'a','f','fill',true),...
                    'lags',(opt.XCMAXLAG+4*delta).*[-1 1]),...
                    1/delta,[],-opt.XCMAXLAG,opt.XCMAXLAG);
                [hdata.path]=deal([indir fs yrdir{i} fs tsdir{i}{j}]);
            else
                hdata=hdata([]);
            end
        end
        if(any(steps==12)) % rotate xc
            % this removes ZZ correlations!
            if(isxc)
                data=rotate_correlations(data);
            else
                hdata=rotate_correlations(hdata);
            end
        end
        
        % write the data
        if(isxc || ~hvsplit)
            if(isempty(data)); continue; end
            writeseizmo(data,'path',...
                [outdir fs yrdir{i} fs tsdir{i}{j} fs]);
        else
            if(~isempty(vdata))
                writeseizmo(vdata,'path',...
                    [outdir fs yrdir{i} fs tsdir{i}{j} fs]);
            end
            if(~isempty(hdata))
                writeseizmo(hdata,'path',...
                    [outdir fs yrdir{i} fs tsdir{i}{j} fs]);
            end
        end
    end
end

% parallel processing takedown & fix verbosity
%matlabpool close; % PARALLEL
seizmoverbose(verbose);

end
