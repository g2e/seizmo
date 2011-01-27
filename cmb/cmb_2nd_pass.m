function [results2]=cmb_2nd_pass(results,sr,varargin)
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
%     option/value pairs to MULTIBANDALIGN to adjust its parameters.  See
%     MULTIBANDALIGN for more details.
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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 26, 2011 at 13:35 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));
if(nargin>2 && mod(nargin,2))
    error('seizmo:cmb_2nd_pass:badNumInputs',...
        'OPTIONS must be paired with a value!');
end

% check results
error(check_cmb_results(results));

% check sample rate
if(nargin<2); sr=[]; end
if(~isempty(sr) && (~isscalar(sr) || sr<=0))
    error('seizmo:cmb_2nd_pass:badInput',...
        'SR must be the new sample rate (in samples/second)!');
end

% loop over each event
cnt=0;
for i=1:numel(results)
    % run name
    runname=regexprep(results(i).runname,'1stPass','2ndPass');
    disp(runname);
    
    % skip if no useralign
    if(isempty(results(i).useralign)); continue; end
    
    % read in data
    data=readseizmo(strcat(results(i).dirname,filesep,...
        {results(i).useralign.data.name}'));
    
    % adjust ground units
    if(any(results(i).usercluster.units~=0))
        units=results0(i).usercluster.units;
        if(any(units==-2))
            data(units==-2)=differentiate(differentiate(data(units==-2)));
        end
        if(any(units==-1))
            data(units==-1)=differentiate(data(units==-1));
        end
        if(any(units==1))
            data(units==1)=integrate(data(units==1));
        end
        if(any(units==2))
            data(units==2)=integrate(integrate(data(units==2)));
        end
    end
    
    % resample
    if(~isempty(sr)); data=syncrates(data,sr); end
    
    % align data using 1stPass results
    arr=results(i).useralign.solution.arr;
    pol=results(i).useralign.solution.pol;
    data=multiply(data,pol);
    data=timeshift(data,-getheader(data,'o')-arr);
    
    % loop over good clusters
    for j=find(results(i).usercluster.good(:)')
        % get cluster info
        sj=num2str(j,'%02d');
        disp(['Aligning cluster ' sj]);
        good=find(results(i).usercluster.T==j ...
            & ~(results(i).outliers.bad));
        
        % census
        pop=numel(good);
        if(pop<3)
            warning('seizmo:cmb_2nd_pass:tooFewGood',...
                ['Cluster ' sj ' has <3 good members. Skipping!']);
            continue;
        end
        
        % extract appropriate corrections
        correct=results(i).corrections;
        correct=fixcorrstruct(correct,good);
        
        % multiband alignment
        tmp=multibandalign(data(good),...
            'runname',[runname '_cluster_' sj],...
            'absxc',false,'estarr',0,'estpol',1,'wgtpow',2,...
            varargin{:});

        % loop over each band in the result to add more info
        for k=1:numel(tmp)
            % matlab bug workaround
            % - see cmb_1st_pass
            if(isempty(tmp(k).usersnr))
                tmp(k).usersnr=[];
            end

            % add run name, quake name, data directory name, syn stuff
            tmp(k).phase=results(i).phase;
            tmp(k).runname=[runname '_cluster_' sj ...
                '_band_' num2str(k,'%02d')];
            tmp(k).dirname=results(i).dirname;
            tmp(k).synthetics=results(i).synthetics;
            tmp(k).earthmodel=results(i).earthmodel;

            % fix corrections
            if(~isempty(tmp(k).usersnr))
                good=find(tmp(k).usersnr.snr>=tmp(k).usersnr.snrcut);
                good(tmp(k).userwinnow.cut)=[];
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
        end

        % save results (all bands together)
        save([runname '_cluster_' sj '_allband_results.mat'],'tmp');

        % export to command line too
        cnt=cnt+1;
        results2(1:k,cnt)=tmp;
    end
end

end


function [s]=fixcorrstruct(s,good)
fields=fieldnames(s);
for i=1:numel(fields)
    if(isstruct(s.(fields{i})))
        s.(fields{i})=fixcorrstruct(s.(fields{i}),good);
    else
        s.(fields{i})=s.(fields{i})(good);
    end
end
end

