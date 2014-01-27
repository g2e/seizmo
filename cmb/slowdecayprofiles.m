function [pf]=slowdecayprofiles(results,azrng,gcrng,odir)
%SLOWDECAYPROFILES    Returns multi-station profile measurements
%
%    Usage:    pf=slowdecayprofiles(results,azrng,gcrng)
%              pf=slowdecayprofiles(results,azrng,gcrng,odir)
%
%    Description:
%     PF=SLOWDECAYPROFILES(RESULTS,AZRNG,GCRNG) takes the relative arrival
%     time and amplitude measurements contained in RESULTS produced by
%     CMB_1ST_PASS, CMB_CLUSTERING, CMB_OUTLIERS, or CMB_2ND_PASS and
%     calculates the slowness and decay rate for a profile of stations
%     within the criteria set by azimuthal range AZRNG and distance range
%     GCRNG.  Note that AZRNG & GCRNG are absolute ranges, meaning an AZRNG
%     of [0 360] will not exclude any stations by azimuth.  AZRNG & GCRNG
%     must both be 2-element vectors of [AZMIN AZMAX] & [GCMIN GCMAX].
%     They are by default [0 360] & [0 180] (ie no exclusion) and are
%     optional.  The output PF is a struct with as many elements as there
%     are profiles found (depends on number of clusters and how many
%     elements the input RESULTS struct had).  The format of the PF struct
%     is described in the Notes section below.
%
%     PF=SLOWDECAYPROFILES(RESULTS,AZRNG,GCRNG,ODIR) sets the output
%     directory where the PF struct is saved.  By default ODIR is '.' (the
%     current directory.  You may set ODIR to FALSE for no written output.
%
%    Notes:
%     - The PF struct is also written to disk as:
%           TIMESTAMP_EVENTDIR_EARTHMODEL_nstn_profiles.mat
%       where TIMESTAMP is the time when the file is written and uses
%       format 30 from the DATESTR function, EVENTDIR is the directory
%       where the waveforms reside, & EARTHMODEL is derived from the
%       results.earthmodel field.  If you give a results struct of 2+
%       elements with differing event directories or earthmodels then they
%       are set to 'misc'.
%     - The PF struct has the following fields:
%       .gcdist         - degree distance difference between stations
%       .azwidth        - azimuthal difference between stations
%       .slow           - horizontal slowness (s/deg)
%       .slowerr        - horizontal slowness standard error
%       .decay          - decay rate
%       .decayerr       - decay rate standard error
%       .cslow          - corrected horizontal slowness***
%       .cslowerr       - corrected horizontal slowness standard error
%       .cdecay         - corrected decay rate
%       .cdecayerr      - corrected decay rate standard error
%       .cluster        - cluster id
%       .kname          - {net stn stream cmp}
%       .st             - [lat lon elev(m) depth(m)]
%       .ev             - [lat lon elev(m) depth(m)]
%       .delaz          - [degdist az baz kmdist]
%       .corrcoef       - max correlation coefficient between waveforms
%       .synthetics     - TRUE if synthetic data (only reflect synthetics)
%       .earthmodel     - model used to make synthetics or 'DATA'
%       .freq           - filter corners of bandpass
%       .phase          - core-diffracted wave type
%       .runname        - name of this run, used for naming output
%       .dirname        - directory containing the waveforms
%       .time           - date string of time of this struct's creation
%
%      *** Correction is different between data and synthetics.  For data
%          the .cslow value is found by subtracting out the corrections
%          (and hence attempts to go from 3D to 1D by removing the lateral
%          heterogeneity).  For synthetics the .cslow value is essentially
%          the opposite (it is corrected to 3D).  So basically:
%                     +----------+---------------+
%                     |   DATA   |   SYNTHETICS  |
%            +--------+----------+---------------+
%            |  .slow |    3D    |       1D      |
%            +--------+----------+---------------+
%            | .cslow |    1D    |       3D      |
%            +--------+----------+---------------+
%
%          To compare the data & sythetics you should compare 3D values of
%          data to 3D values of synthetics or 1D values of data to 1D
%          values of synthetics.  Drawing conclusions from comparison of 3D
%          to 1D is not recommended (except to see the affect corrections
%          have on data or synthetics).
%
%    Examples:
%     % Return station profiles with an azimuth
%     % of 200-220deg and a distance of 90-160deg:
%     pf=slowdecayprofiles(results,[200 220],[90 160])
%
%    See also: SLOWDECAYPAIRS, CMB_2ND_PASS, CMB_OUTLIERS, CMB_1ST_PASS,
%              CMB_CLUSTERING, PREP_CMB_DATA, PLOT_CMB_MEASUREMENTS,
%              MAP_CMB_PROFILES, PLOT_CMB_PDF

%     Version History:
%        Dec. 12, 2010 - initial version
%        Jan. 18, 2011 - update for results struct standardization, added
%                        corrections & correlation coefficients to output,
%                        time is now a string, require common event
%        Jan. 23, 2011 - fix indexing bug
%        Jan. 26, 2011 - pass on new .synthetics & .earthmodel fields,
%                        .cslow depends on .synthetics, added Notes
%                        about PF struct format
%        Jan. 29, 2011 - save output, fix corrections bug
%        Jan. 31, 2011 - allow no output, odir input, better checks
%        Feb.  5, 2011 - fix bug when no output specified
%        Feb. 12, 2011 - include snr-based arrival time error
%        Feb. 17, 2011 - fixed decay constant error
%        Mar.  1, 2011 - combined write rather than individually, added
%                        notes about output
%        Mar.  2, 2011 - earthmodel in output name
%        Mar. 18, 2011 - handle raypaths in correction info
%        Mar. 30, 2011 - doc update
%        Apr. 22, 2011 - update for finalcut field
%        Mar.  2, 2012 - handle unset earthmodel bugfix, octave save ascii
%                        workaround, .mat output name includes eventdir
%        Mar.  5, 2012 - allow no written output
%        Mar.  6, 2012 - basic .weights field support
%        Oct. 11, 2012 - drop corrections field (huge)
%        July 25, 2013 - directory input (reads indir/*.mat)
%        Jan. 27, 2014 - put unique output file detection in while loop,
%                        abs path fix & reduced filesep/fullfile calls
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 27, 2014 at 13:35 GMT

% todo:

% check nargin
error(nargchk(1,4,nargin));

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

% check results struct
error(check_cmb_results(results));

% default azrng, gcrng, odir
if(nargin<2 || isempty(azrng)); azrng=[0 360]; end
if(nargin<3 || isempty(gcrng)); gcrng=[0 180]; end
if(nargin<4 || isempty(odir)); odir='.'; end
if(islogical(odir) && isscalar(odir) && odir); odir='.'; end

% check azrng & gcrng
if(~isreal(azrng) || numel(azrng)~=2)
    error('seizmo:slowdecayprofiles:badInput',...
        'AZRNG must be [AZMIN AZMAX]!');
elseif(any(abs(azrng)>540))
    error('seizmo:slowdecayprofiles:badInput',...
        'Keep AZRNG within +/-540deg!');
elseif(~isreal(gcrng) || numel(gcrng)~=2)
    error('seizmo:slowdecayprofiles:badInput',...
        'GCRNG must be [GCMIN GCMAX]!');
elseif(~isstring(odir) && ~(islogical(odir) && isscalar(odir)))
    error('seizmo:slowdecayprofiles:badInput',...
        'ODIR must be a string or TRUE/FALSE!');
end
if(islogical(odir)); out=odir; else out=true; end

% make sure odir exists (create it if it does not)
if(out)
    [ok,msg,msgid]=mkdir(odir);
    if(~ok)
        warning(msgid,msg);
        error('seizmo:slowdecayprofiles:pathBad',...
            'Cannot create directory: %s',odir);
    end
elseif(~out && ~nargout)
    error('seizmo:slowdecayprofiles:badInput',...
        'Output variable must be assigned when no written output!');
end

% verbosity
%verbose=seizmoverbose;

% loop over every result
for a=1:numel(results)
    % skip if results.useralign is empty
    if(isempty(results(a).useralign)); continue; end
    
    % number of records
    %nrecs=numel(results(a).useralign.data);
    
    % extract header details
    [st,ev,delaz,kname]=getheader(results(a).useralign.data,...
        'st','ev','delaz','kname');
    
    % check event info matches
    ev=unique(ev,'rows');
    if(size(ev,1)>1)
        error('seizmo:slowdecaypairs:badInput',...
            'EVENT location varies between records!');
    end
    
    % corrected relative arrival times and amplitudes
    rtime=results(a).useralign.solution.arr;
    if(results(a).synthetics)
        % we add corrections here to go from 1D to 3D
        switch results(a).phase
            case 'Pdiff'
                crtime=results(a).useralign.solution.arr...
                    +results(a).corrections.ellcor...
                    +results(a).corrections.crucor.prem...
                    +results(a).corrections.mancor.hmsl06p.upswing;
            case {'SHdiff' 'SVdiff'}
                crtime=results(a).useralign.solution.arr...
                    +results(a).corrections.ellcor...
                    +results(a).corrections.crucor.prem...
                    +results(a).corrections.mancor.hmsl06s.upswing;
        end
    else % data
        % we subtract corrections here to go from 3D to 1D
        switch results(a).phase
            case 'Pdiff'
                crtime=results(a).useralign.solution.arr...
                    -results(a).corrections.ellcor...
                    -results(a).corrections.crucor.prem...
                    -results(a).corrections.mancor.hmsl06p.upswing;
            case {'SHdiff' 'SVdiff'}
                crtime=results(a).useralign.solution.arr...
                    -results(a).corrections.ellcor...
                    -results(a).corrections.crucor.prem...
                    -results(a).corrections.mancor.hmsl06s.upswing;
        end
    end
    snr=results(a).usersnr.snr;
    snr=snr(snr>=results(a).usersnr.snrcut);
    snr(results(a).userwinnow.cut)=[];
    if(isfield(results(a),'finalcut')); snr=snr(results(a).finalcut); end
    rtimeerr=sqrt((results(a).useralign.solution.arrerr).^2 ...
        +(max(1./results(a).filter.corners)...
        ./(2*pi).*snr2phaseerror(snr)).^2);
    rampl=results(a).useralign.solution.amp;
    crampl=results(a).useralign.solution.amp...
        ./results(a).corrections.geomsprcor;
    ramplerr=results(a).useralign.solution.amperr;
    
    % get cluster indexing
    cidx=results(a).usercluster.T;
    good=results(a).usercluster.good';
    
    % get outliers
    outliers=results(a).outliers.bad;
    
    % loop over "good" clusters
    cnt=0;
    for b=find(good)
        % indices of members
        m=cidx==b & ~outliers;
        
        % get stations in range (need the indices)
        delaz(:,2)=delaz(:,2)-360*ceil((delaz(:,2)-azrng(2))/360);
        idx=find(m & delaz(:,1)>=gcrng(1) & delaz(:,1)<=gcrng(2) ...
            & delaz(:,2)>=azrng(1) & delaz(:,2)<=azrng(2));
        nsta=numel(idx);
        
        % skip if none/one
        if(nsta<2); continue; end
        
        % initialize struct
        cnt=cnt+1;
        tmp(cnt)=struct('gcdist',[],'azwidth',[],...
            'slow',[],'slowerr',[],'decay',[],'decayerr',[],...
            'cslow',[],'cslowerr',[],'cdecay',[],'cdecayerr',[],...
            'cluster',b,'kname',[],'st',[],'ev',[],'delaz',[],...
            ...%'corrections',[],...
            'corrcoef',[],...
            'synthetics',results(a).synthetics,...
            'earthmodel',results(a).earthmodel,...
            'freq',results(a).filter.corners,'phase',results(a).phase,...
            'runname',results(a).runname,'dirname',results(a).dirname,...
            'time',datestr(now));
        
        % insert known info
        tmp(cnt).kname=kname(idx,:);
        tmp(cnt).st=st(idx,:);
        tmp(cnt).ev=ev;
        tmp(cnt).delaz=delaz(idx,:);
        
        % great circle distance and width
        tmp(cnt).gcdist=max(delaz(idx,1))-min(delaz(idx,1));
        tmp(cnt).azwidth=max(delaz(idx,2))-min(delaz(idx,2));
        
        % corrections
        %tmp(cnt).corrections=fixcorrstruct(results(a).corrections,idx);
        
        % correlation coefficients
        tmp(cnt).corrcoef=...
            submat(ndsquareform(results(a).useralign.xc.cg),1:2,idx,3,1);
        
        % find slowness & decay rate
        if(isfield(results(a),'weights'))
            try
                [m,covm]=wlinem(delaz(idx,1),rtime(idx),1,...
                    diag(rtimeerr(idx).^2),diag(results(a).weights(idx)));
                tmp(cnt).slow=m(2);
                tmp(cnt).slowerr=sqrt(covm(2,2));
                [m,covm]=wlinem(delaz(idx,1),crtime(idx),1,...
                    diag(rtimeerr(idx).^2),diag(results(a).weights(idx)));
                tmp(cnt).cslow=m(2);
                tmp(cnt).cslowerr=sqrt(covm(2,2));
                [m,covm]=wlinem(delaz(idx,1),log(rampl(idx)),1,...
                    diag((log(rampl(idx)+ramplerr(idx))...
                    -log(rampl(idx))).^2),diag(results(a).weights(idx)));
                tmp(cnt).decay=m(2);
                tmp(cnt).decayerr=sqrt(covm(2,2));
                [m,covm]=wlinem(delaz(idx,1),log(crampl(idx)),1,...
                    diag((log(crampl(idx)+ramplerr(idx))...
                    -log(crampl(idx))).^2),diag(results(a).weights(idx)));
                tmp(cnt).cdecay=m(2);
                tmp(cnt).cdecayerr=sqrt(covm(2,2));
            catch
                error('seizmo:slowdecayprofiles:corruptResults',...
                    'Something is wrong with the .weights field!');
            end
        else
            [m,covm]=wlinem(delaz(idx,1),rtime(idx),1,...
                diag(rtimeerr(idx).^2));
            tmp(cnt).slow=m(2);
            tmp(cnt).slowerr=sqrt(covm(2,2));
            [m,covm]=wlinem(delaz(idx,1),crtime(idx),1,...
                diag(rtimeerr(idx).^2));
            tmp(cnt).cslow=m(2);
            tmp(cnt).cslowerr=sqrt(covm(2,2));
            [m,covm]=wlinem(delaz(idx,1),log(rampl(idx)),1,...
                diag((log(rampl(idx)+ramplerr(idx))-log(rampl(idx))).^2));
            tmp(cnt).decay=m(2);
            tmp(cnt).decayerr=sqrt(covm(2,2));
            [m,covm]=wlinem(delaz(idx,1),log(crampl(idx)),1,...
                diag((log(crampl(idx)+ramplerr(idx))...
                -log(crampl(idx))).^2));
            tmp(cnt).cdecay=m(2);
            tmp(cnt).cdecayerr=sqrt(covm(2,2));
        end
    end
    
    % skip if none
    if(~cnt); continue; end
    
    % output
    if(~exist('pf','var'))
        pf=tmp;
    else
        pf=[pf; tmp];
    end
end

% save profiles
if(out && exist('pf','var'))
    % strings for saving
    wfdirs={results.dirname}';
    wfdirs(cellfun('isempty',wfdirs))=[];
    wfdir=unique(wfdirs);
    if(~isscalar(wfdir))
        wfdir='misc';
    else
        [a,b,c]=fileparts(char(wfdir));
        wfdir=[b c];
    end
    emods={results.earthmodel}';
    emods(cellfun('isempty',emods))=[];
    emod=unique(emods);
    if(~isscalar(emod)); emod='misc'; else emod=char(emod); end
    
    % avoid clobber by waiting until unique time
    while(exist([odir fs datestr(now,30) '_' wfdir '_' emod ...
            '_nstn_profiles.mat'],'file'))
        pause(1);
    end
    
    % saving
    if(isoctave)
        save([odir fs datestr(now,30) '_' wfdir '_' emod ...
            '_nstn_profiles.mat'],'-7','pf');
    else % matlab
        save([odir fs datestr(now,30) '_' wfdir '_' emod ...
            '_nstn_profiles.mat'],'pf');
    end
end

% check for output
if(nargout && ~exist('pf','var'))
    error('seizmo:slowdecayprofiles:noPairs',...
        'Not enough stations meet the specified profile criteria!');
elseif(~nargout && exist('pf','var'))
    clear pf;
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
