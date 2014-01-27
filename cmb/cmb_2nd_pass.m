function [varargout]=cmb_2nd_pass(results,sr,varargin)
%CMB_2ND_PASS    Narrow-band core-diff relative arrivals + amplitudes
%
%    Usage:    results=cmb_2nd_pass(results)
%              results=cmb_2nd_pass(results,sr)
%              results=cmb_2nd_pass(results,sr,'option',value,...)
%
%    Description:
%     RESULTS=CMB_2ND_PASS(RESULTS) aligns the data in RESULTS (returned
%     from CMB_1ST_PASS, CMB_CLUSTERING, or CMB_OUTLIERS) for a series of
%     25 frequency bands from 80s to 8s.  See the last usage form to adjust
%     the frequency bands and filters.
%
%     RESULTS=CMB_2ND_PASS(RESULTS,SR) resamples records to a sample rate
%     of SR.  SR must be in Hz (ie SR==5 is 5Hz sample rate).  The default
%     is no resampling.
%
%     RESULTS=CMB_2ND_PASS(RESULTS,SR,'OPTION',VALUE,...) passes
%     option/value pairs.  Any options allowed by MULTIBANDALIGN may be
%     set.  See MULTIBANDALIGN for a comprehensive list of those.  Some
%     useful options:
%      'odir'    - directory to save all output (current dir is default)
%      'figdir'  - directory to save figures (overrides odir option)
%      'distrng' - specify a distance window as [DISTMIN DISTMAX]
%      'azrng'   - specify an azimuth window as [AZMIN AZMAX]
%
%    Notes:
%     - If you have corrected the ground units of the underlying data in
%       the .dirname directory using .usercluster.units then you should
%       also set .usercluster.units to all zeros!
%
%    Examples:
%     % This is the typical usage case for me:
%     sr=5; % 5sps
%     results2=cmb_2nd_pass(results,sr);
%
%     % Nothing is written to disk if the
%     % output directories are set to false:
%     results=cmb_2nd_pass(results,5,'odir',false,'figdir',false);
%
%    See also: PREP_CMB_DATA, CMB_1ST_PASS, CMB_CLUSTERING, CMB_OUTLIERS,
%              SLOWDECAYPAIRS, SLOWDECAYPROFILES, MAP_CMB_PROFILES

%     Version History:
%        Dec. 12, 2010 - added docs
%        Jan. 14, 2011 - corrections are properly edited to match output,
%                        support for .adjustclusters.units, fixed outlier
%                        bug (forgot to remove them)
%        Jan. 18, 2011 - update for results struct and multibandalign
%                        updates, edit output names & runname to keep
%                        further output informative if we cluster/outlier a
%                        narrow band result
%        Jan. 26, 2011 - .synthetics & .earthmodel fields, 2-digit cluster
%        Jan. 29, 2011 - output now has creation datetime string prepended
%        Jan. 31, 2011 - fix bug on too few high snr, allow no output, odir
%                        catching
%        Mar. 18, 2011 - handle raypaths in correction info
%        Mar. 25, 2011 - handle new radiation pattern corrections
%        Apr.  8, 2011 - more informative message if you forgot to input sr
%        Apr. 11, 2011 - improve docs, add gcrng/azrng options
%        Apr. 17, 2011 - optimizations for multibandalign, clear gc/azrng
%        Apr. 22, 2011 - update for finalcut field
%        May  20, 2011 - add some code to workaround matlab write bug
%        Mar.  1, 2012 - octave ascii save workaround
%        Mar.  5, 2012 - allow no written output, odir/figdir bugfix
%        Feb. 14, 2013 - bugfix figdir handling
%        Mar. 11, 2013 - directory input (reads indir/*.mat), selection
%                        list, advanced clustering commented out
%        Jan. 27, 2014 - abs path fix & reduced filesep/fullfile calls
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 27, 2014 at 13:35 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));

% directory separator
fs=filesep;

% handle directory input
if(isstring(results))
    if(~isabspath(results)); results=[pwd fs results]; end
    if(isdir(results))
        files=xdir([results fs '*.mat']);
        clear results;
        for i=1:numel(files)
            results(i)=load([files(i).path files(i).name]);
        end
    end
end

% check results
error(check_cmb_results(results));

% check sample rate
if(nargin<2); sr=[]; end
if(~isempty(sr) && (~isnumeric(sr) || ~isscalar(sr) || sr<=0))
    error('seizmo:cmb_2nd_pass:badInput',...
        'SR must be the new sample rate (in samples/second)!');
end

% require parameter/value pairs after sample rate
if(nargin>2 && mod(nargin,2))
    error('seizmo:cmb_2nd_pass:badNumInputs',...
        'OPTIONS must be paired with a value!');
end

% default/extract odir
odir='.';
figdir=[];
gcrng=[0 180];
azrng=[0 360];
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:cmb_2nd_pass:badInput',...
        'OPTION must all be strings!');
end
keep=true(nargin-2,1);
for i=1:2:nargin-2
    switch lower(varargin{i})
        case {'outdir' 'odir'}
            if(isempty(varargin{i+1}))
                odir='.';
                keep(i:i+1)=false;
            else
                odir=varargin{i+1};
                if(isempty(figdir))
                    varargin{i}='figdir';
                    figdir=odir;
                else
                    keep(i:i+1)=false;
                end
            end
        case {'gcrng' 'distrng'}
            if(~isnumeric(varargin{i+1}) || numel(varargin{i+1})~=2 ...
                    || varargin{i+1}(1)>varargin{i+1}(2) ...
                    || any(varargin{i+1})<0)
                error('seizmo:cmb_2nd_pass:badInput',...
                    'DISTRNG must be given as [DISTMIN DISTMAX]!');
            end
            gcrng=varargin{i+1};
            keep(i:i+1)=false;
        case 'azrng'
            if(~isnumeric(varargin{i+1}) || numel(varargin{i+1})~=2 ...
                    || varargin{i+1}(1)>varargin{i+1}(2))
                error('seizmo:cmb_2nd_pass:badInput',...
                    'AZRNG must be given as [AZMIN AZMAX]!');
            end
            azrng=varargin{i+1};
            keep(i:i+1)=false;
        case 'figdir'
            figdir=varargin{i+1};
    end
end
varargin=varargin(keep);
if(isempty(figdir)); figdir=odir; end

% fix TRUE directory input
if(islogical(odir) && isscalar(odir) && odir); odir='.'; end
if(islogical(figdir) && isscalar(figdir) && figdir); figdir=odir; end

% check odir
if(~isstring(odir) && ~(islogical(odir) && isscalar(odir)))
    error('seizmo:cmb_2nd_pass:badInput',...
        'ODIR must be a string or TRUE/FALSE!');
elseif(~isstring(figdir) && ~(islogical(figdir) && isscalar(figdir)))
    error('seizmo:cmb_2nd_pass:badInput',...
        'FIGDIR must be a string or TRUE/FALSE!');
end
if(islogical(odir)); out=odir; else out=true; end
%if(islogical(figdir)); figout=odir; else figout=true; end
% UNCOMMENT THE ABOVE IF FIGURE HANDLING IS EVER INSERTED BELOW HERE

% create odir if not there
if(out)
    [ok,msg,msgid]=mkdir(odir);
    if(~ok)
        warning(msgid,msg);
        error('seizmo:cmb_2nd_pass:pathBad',...
            'Cannot create directory: %s',odir);
    end
elseif(~out && ~nargout)
    error('seizmo:cmb_2nd_pass:badInput',...
        'Output variable must be assigned when no written output!');
end

% select events
datelist=char({results.runname}.');
s=listdlg('PromptString','Select events:',...
          'InitialValue',1:numel(results),...
          'ListSize',[170 300],...
          'ListString',datelist);

% error if none selected
if(isempty(s))
    error('seizmo:cmb_2nd_pass:noDirsSelected',...
        'No earthquakes selected!');
end

% loop over each event
cnt=0;
for i=1:numel(s)
    % run name
    runname=regexprep(results(s(i)).runname,'1stPass','2ndPass');
    disp(runname);
    
    % skip if no useralign
    if(isempty(results(s(i)).useralign)); continue; end
    
    % read in data
    data=readseizmo(strcat(results(s(i)).dirname,fs,...
        {results(s(i)).useralign.data.name}'));
    
    % get some header info
    [gcarc,az]=getheader(data,'gcarc','az');
    
    % adjust ground units
    %if(any(results(s(i)).usercluster.units~=0))
    %    units=results(s(i)).usercluster.units;
    %    if(any(units==-2))
    %        data(units==-2)=differentiate(differentiate(data(units==-2)));
    %    end
    %    if(any(units==-1))
    %        data(units==-1)=differentiate(data(units==-1));
    %    end
    %    if(any(units==1))
    %        data(units==1)=integrate(data(units==1));
    %    end
    %    if(any(units==2))
    %        data(units==2)=integrate(integrate(data(units==2)));
    %    end
    %end
    
    % resample
    if(~isempty(sr)); data=syncrates(data,sr); end
    
    % align data using 1stPass results
    arr=results(s(i)).useralign.solution.arr;
    pol=results(s(i)).useralign.solution.pol;
    data=multiply(data,pol);
    data=timeshift(data,-getheader(data,'o')-arr);
    
    % loop over good clusters
    for j=find(results(s(i)).usercluster.good(:)')
        % get cluster info
        sj=num2str(j,'%02d');
        disp(['Aligning cluster ' sj]);
        good=find(results(s(i)).usercluster.T==j ...
            & ~(results(s(i)).outliers.bad) ...
            & gcarc>=gcrng(1) & gcarc<=gcrng(2) ...
            & ((az>=azrng(1) & az<=azrng(2)) ...
            | (az>=azrng(1)-360 & az<=azrng(2)-360) ...
            | (az>=azrng(1)+360 & az<=azrng(2)+360)));
        
        % census
        pop=numel(good);
        if(pop<3)
            warning('seizmo:cmb_2nd_pass:tooFewGood',...
                ['Cluster ' sj ' has <3 good members. Skipping!']);
            continue;
        end
        
        % extract appropriate corrections
        correct=results(s(i)).corrections;
        correct=fixcorrstruct(correct,good);
        
        % apply sign corrections for radiation pattern
        data(good)=multiply(data(good),sign(correct.radpatcor));
        
        % multiband alignment
        tmp=multibandalign(data(good),...
            'runname',[runname '_cluster_' sj],'lags',1/2,'minlag',.01,...
            'npeaks',1,'estarr',0,'estpol',1,'wgtpow',2,...
            'method','reweight',varargin{:});

        % loop over each band in the result to add more info
        for k=1:numel(tmp)
            % matlab bug workaround
            % - really strange...saving fails to
            %   write all fields if we dont do this
            if(isempty(tmp(k).usersnr)); tmp(k).usersnr=[]; end

            % add run name, quake name, data directory name, syn stuff
            tmp(k).phase=results(s(i)).phase;
            tmp(k).runname=[runname '_cluster_' sj ...
                '_band_' num2str(k,'%02d')];
            tmp(k).dirname=results(s(i)).dirname;
            tmp(k).synthetics=results(s(i)).synthetics;
            tmp(k).earthmodel=results(s(i)).earthmodel;

            % fix corrections
            if(~isempty(tmp(k).usersnr) && ~isempty(tmp(k).userwinnow))
                good=find(tmp(k).usersnr.snr>=tmp(k).usersnr.snrcut);
                good(tmp(k).userwinnow.cut)=[];
                if(isfield(tmp(k),'finalcut'))
                    good=good(tmp(k).finalcut);
                end
                tmp(k).corrections=fixcorrstruct(correct,good);
            else
                tmp(k).corrections=fixcorrstruct(correct,[]);
            end
            
            % add in clustering info (all belong to 1 cluster)
            if(isempty(tmp(k).useralign))
                nrecs=0;
            else
                nrecs=numel(tmp(k).useralign.data);
            end
            tmp(k).usercluster.T=ones(nrecs,1);
            tmp(k).usercluster.units=zeros(nrecs,1);
            tmp(k).usercluster.good=true;
            tmp(k).usercluster.color=[1 0 0]; % default to red
            
            % add in outlier info (no outliers)
            tmp(k).outliers.bad=false(nrecs,1);
            
            % time of run
            tmp(k).time=datestr(now);
            
            % matlab bug workaround
            % - really strange...saving fails to
            %   write all fields if we dont do this
            if(isempty(tmp(k).userwinnow)); tmp(k).userwinnow=[]; end
        end

        % save results (all bands together)
        if(out)
            if(isoctave)
                save([odir fs datestr(now,30) '_' runname '_cluster_' ...
                    sj '_allband_results.mat'],'-7','tmp');
            else % matlab
                save([odir fs datestr(now,30) '_' runname '_cluster_' ...
                    sj '_allband_results.mat'],'tmp');
            end
        end

        % export to command line too
        if(nargout)
            cnt=cnt+1;
            varargout{1}(1:k,cnt)=tmp;
        end
    end
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

