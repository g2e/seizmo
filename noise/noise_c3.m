function [sdata,c3]=noise_c3(data,varargin)
%NOISE_C3    Correlates the Coda of Correlations
%
%    Usage:    ncfs=noise_c3(ncfs)
%              ncfs=noise_c3(ncfs,'opt1',val,...,'optN',val)
%
%    Description:
%     NCFS=NOISE_C3(NCFS) windows the Rayleigh & Love waves codas from the
%     seismic noise correlations (NCFs) in the SEIZMO dataset NCFS and
%     cross correlates those codas to recapture the empirical response of
%     the Earth between seismic stations (aka empirical Green's functions).
%     This operation is useful because NCFs are often biased by the poor
%     distribution of microseismic noise (e.g., non-ambience) while codas
%     are comprised of the scattered wavefield that tends to be more
%     directionally diffuse and therefore less biased.  Another interesting
%     possibility with this function is the recovery of Green's functions
%     from asynchronous observations (e.g., Ma & Beroza 2012).
%
%     NCFS=NOISE_C3(NCFS,'OPT1',VAL,...,'OPTN',VAL) allows modifying
%     several of the parameters used in the C3 method.  The following are
%     configurable:
%      XCMAXLAG  - maximum lag time of correlograms in sec [4000]
%      NORMXC    - normalize the cross correlations [true]
%      COHERENCY - compute coherency correlograms [false]
%      FDOUT     - output correlations in the frequency-domain [false]
%      NOAUTO    - skip autocorrelations [false]
%      VRAYL     - Rayleigh wave velocity.  Default is 2.6km/s.
%      WINOFF    - Coda offset from TRAYL in seconds.  Default is 300s.
%      WINLEN    - Coda length in seconds.  Default is 1200s.
%      XCCMP     - Correlation component to analyze.  Default: 'sym'.
%                  Options: 'causal', 'acausal', 'both' & 'sym'
%      ZTRANS    - Fisher's transform for stacking?  Default is true.
%      MINCODA   - minimum coda stations for coda correlation [1]
%      STNLIST   - coda station list []
%                  Station list should be a cell array of strings as:
%                   {'NET.STN.HOLE.CMP' ...}
%
%    Notes:
%     - C3 References:
%        Stehly et al 2008, JGR, doi:10.1029/2008JB005693
%        de Ridder et al 2009, SEG, pp. 1622-1626
%        Garnier and Papanicolaou 2009, SIAM, pp. 396-437
%        Froment et al 2011, CRG, pp. 623-632
%        Ma & Beroza 2012, GRL, doi:10.1029/2011g1050755
%        Zhang & Yang 2013, JGR, doi:10.1002/jgrb.50186
%     - Redundant correlation pairs in NCFS will generate an error!  That
%       means you should only input one set of correlations and that set
%       must not include reversed correlations (these are redundant too!).
%       Using the function NO_REDUNDANT_CORRELATIONS can help to resolve
%       this issue (but no guarantees!).
%
%    Header changes: SCALE, DEP*
%     A, F, T0-1, T3-4, REFTIME may change (if they vary amongst records)
%     VRAYL  => RESP0
%     WINOFF => RESP1
%     WINLEN => RESP2
%     NCODA  => RESP3
%     XCCMP  => KUSER2
%
%    Examples:
%     % Say the variable "ncfs" contains the noise correlation stacks for
%     % an array.  To check how the C3 method improves the apparent noise
%     % distribution in the period range of 10=20s use the FSS function:
%     fss(ncfs,50,101,[1/20 1/10]);
%     fss(noise_c3(ncfs),50,101,[1/20 1/10]);
%
%    See also: NOISE_PROCESS, CORRELATE, NOISE_STACK, STACK2STACK,
%              NOISE_BESSEL_FIT

%     Version History:
%        June 25, 2014 - initial version
%        July 11, 2014 - fd is converted to complex so Fisher transform
%                        works properly (FISHER was updated)
%        July 17, 2014 - bugfix: solofun needs func handles not strings,
%                        bugfix: fix option naming, bugfix: fix creation of
%                        option defaults, output is column vector
%        July 21, 2014 - bugfix: verbosity restored to previous state
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 21, 2014 at 11:15 GMT

% todo:
% - result does not look impressive so we need to do a step by step check
%   - where is the coda window (make a nice plot)
%     - are symmetric correct?
%   - all correlations for a pair
%   - stacks of each xccmp type

% check nargin
error(nargchk(1,inf,nargin));
if(~mod(nargin,2))
    error('seizmo:noise_c3:badInput',...
        'Unpaired option/value pair given!');
end

% check structure
error(seizmocheck(data,'dep'));

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% turn off verbosity
verbose=seizmoverbose(false);

% detail message
if(verbose); disp('CORRELATING THE CODA OF CORRELATIONS'); end

% attempt header check
try
    % check headers
    data=checkheader(data,...
        'MULCMP_DEP','ERROR',...
        'FALSE_LEVEN','ERROR',...
        'XYZ_IFTYPE','ERROR',...
        'UNSET_ST_LATLON','ERROR',...
        'UNSET_EV_LATLON','ERROR');
    
    % turn off header checking
    oldcheckheaderstate=checkheader_state(false);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % toggle verbosity back
    seizmoverbose(verbose);
    
    % rethrow error
    error(lasterror);
end

% attempt c3 processing
try
    % necessary header info
    [kuser,iftype,sdelta]=getheader(data,'kuser','iftype id','sdelta');
    
    % require all correlograms
    if(~all(strcmp(kuser(:,1),'MASTER') & strcmp(kuser(:,2),'SLAVE')))
        error('seizmo:noise_c3:badInput',...
            'Some records appear not to be correlations!');
    end
    
    % convert spectral filetypes to time domain
    rlim=strcmpi(iftype,'irlim');
    amph=strcmpi(iftype,'iamph');
    if(any(rlim | amph))
        data(rlim | amph)=dephase_correlations(data(rlim | amph),false);
        data(rlim | amph)=idft(data(rlim | amph));
        data(rlim | amph)=multiply(data(rlim | amph),sdelta(rlim | amph));
        data(rlim | amph)=unwrap_correlations(data(rlim | amph));
    end
    
    % parse/check parameters
    opt=noise_c3_parameters(varargin{:});
    
    % require no redundant correlations
    if(numel(data)~=numel(no_redundant_correlations(data)))
        error('seizmo:noise_c3:badInput',...
            'Redundant correlations are not allowed!');
    end
    
    % briefly split off auto correlations
    [adata,data]=split_auto_correlations(data);
    
    % require cross correlations
    if(isempty(data))
        error('seizmo:noise_c3:badInput',...
            'NCFS must contain cross-correlations!');
    end
    
    % now get reversed correlations to complete the cross matrix
    data=[adata; data; reverse_correlations(data)];
    
    % number of records
    nrecs=numel(data);
    
    % necessary header info
    [gcarc,mnm,snm,scale,t0utc,t1utc,delta,b,e]=getheader(data,...
        'gcarc','kt','kname','scale','t0 utc','t1 utc','delta','b','e');
    scale(isnan(scale))=1;
    t0utc=cell2mat(t0utc);
    t1utc=cell2mat(t1utc);
    
    % delta values must match for all records
    if(numel(unique(delta))~=1)
        error('seizmo:noise_c3:deltaMismatch',...
            'DELTA (sample spacing) must match for all records!');
    end
    delta=delta(1);
    
    % extract coda
    tRAYL=gcarc*6371*pi/180/opt.C3VRAYL;
    ncodapts=ceil(opt.C3WINLEN/delta)+1;
    switch opt.C3XCCMP
        case 'causal'
            coda{1}=cut(data,tRAYL+opt.C3WINOFF,'n',ncodapts,'pad',true);
        case 'acausal'
            coda{1}=cut(changeheader(reverse(data),'b',-e,'e',-b),...
                tRAYL+opt.C3WINOFF,'n',ncodapts,'pad',true);
        case {'both' 'xboth'}
            coda{1}=cut(data,tRAYL+opt.C3WINOFF,'n',ncodapts,'pad',true);
            coda{2}=cut(changeheader(reverse(data),'b',-e,'e',-b),...
                tRAYL+opt.C3WINOFF,'n',ncodapts,'pad',true);
        case 'sym'
            coda{1}=cut(symmetric_correlations(data),tRAYL+opt.C3WINOFF,...
                'n',ncodapts,'pad',true);
    end
    nsets=numel(coda);
    
    % station list
    mnm=lower(strcat(mnm(:,1),'.',mnm(:,2),'.',mnm(:,3),'.',mnm(:,4)));
    snm=lower(strcat(snm(:,1),'.',snm(:,2),'.',snm(:,3),'.',snm(:,4)));
    [stns,idx1,idx2]=unique([mnm;snm]);
    
    % number of stations
    nstns=numel(stns);
    
    % max number of output correlations
    npairs=nstns*(nstns+1-2*opt.NOAUTO)/2;

    % detail message
    if(verbose); print_time_left(0,npairs); thispair=0; end
    
    % double over indexing for easy master/slave work
    idx2=reshape(idx2,nrecs,2);
    
    % master station list
    if(isempty(opt.C3STNLIST))
        midx=(1:nstns).';
    else % subset
        [midx,midx]=intersect(stns,opt.C3STNLIST);
        if(isempty(midx))
            error('seizmo:noise_c3:noCodaStns',...
                'No coda match station list in STNLIST!');
        end
    end
    
    % loop over "master" slave
    cnt=0; % c3 counter
    for s1=1:nstns-opt.NOAUTO
        % find coda with s1 as slave & master in master list
        in1=idx2(:,2)==s1 & ismember(idx2(:,1),midx);
        if(~any(in1)); continue; end
        
        % loop over "slave" slave
        for s2=s1+opt.NOAUTO:nstns
            % find coda with s2 as slave & master in master list
            in2=idx2(:,2)==s2 & ismember(idx2(:,1),midx);
            if(~any(in2)); continue; end
            
            % only allow slave sets that share a master
            midx12=intersect(idx2(in1,1),idx2(in2,1));
            if(isempty(midx12)); continue; end
            
            % organize coda sets by master index
            inb1=find(idx2(:,2)==s1 & ismember(idx2(:,1),midx12));
            inb2=find(idx2(:,2)==s2 & ismember(idx2(:,1),midx12));
            [cidx1,cidx1]=sort(idx2(inb1,1));
            [cidx2,cidx2]=sort(idx2(inb2,1));
            inb1=inb1(cidx1);
            inb2=inb2(cidx2);
            nin=numel(inb1);
            
            % enough coda
            if(nin<opt.MINCODA); continue; end
            
            % loop over coda sets
            c3=cell(nsets,1);
            for cmp=1:nsets
                % frequency- or time-domain?
                if(isempty(opt.FDOUT)) % time
                    % correlate & interpolate
                    c3{cmp}=interpolate(correlate(coda{cmp}(inb1),...
                        coda{cmp}(inb2),'reltime',opt.COHERENCY{:},...
                        opt.NORMXC{:},(opt.XCMAXLAG+4*delta).*[-1 1]),...
                        1/delta,[],-opt.XCMAXLAG,opt.XCMAXLAG);
                else % frequency
                    % correlate
                    c3{cmp}=correlate(coda{cmp}(inb1),...
                        coda{cmp}(inb2),'reltime',opt.COHERENCY{:},...
                        opt.NORMXC{:},opt.FDOUT{:});
                    
                    % convert fd to cplx
                    c3{cmp}=solofun(c3{cmp},@(x)complex(x(:,1),x(:,2)));
                end
                
                % apply Fisher's transform
                if(opt.C3ZTRANS); c3{cmp}=solofun(c3{cmp},@fisher); end
                
                % multiply by scale
                c3{cmp}=multiply(c3{cmp},scale(inb1).*scale(inb2));
            end
            
            % cross acausal/causal coda correlation
            switch opt.C3XCCMP
                case 'xboth'
                    % frequency- or time-domain?
                    if(isempty(opt.FDOUT)) % time
                        % correlate between causal & acausal & interpolate
                        c3{3}=interpolate(correlate(...
                            coda{1}(inb1),coda{2}(inb2),'reltime',...
                            opt.NORMXC{:},opt.COHERENCY{:},...
                            (opt.XCMAXLAG+4*delta).*[-1 1]),...
                            1/delta,[],-opt.XCMAXLAG,opt.XCMAXLAG);
                        c3{4}=interpolate(correlate(...
                            coda{2}(inb1),coda{1}(inb2),'reltime',...
                            opt.NORMXC{:},opt.COHERENCY{:},...
                            (opt.XCMAXLAG+4*delta).*[-1 1]),...
                            1/delta,[],-opt.XCMAXLAG,opt.XCMAXLAG);
                    else % frequency
                        % correlate between causal & acausal
                        c3{3}=correlate(coda{1}(inb1),coda{2}(inb2),...
                            'reltime',opt.COHERENCY{:},...
                            opt.NORMXC{:},opt.FDOUT{:});
                        c3{4}=correlate(coda{2}(inb1),coda{1}(inb2),...
                            'reltime',opt.COHERENCY{:},...
                            opt.NORMXC{:},opt.FDOUT{:});
                        
                        % convert fd to cplx
                        c3{3}=solofun(c3{3},@(x)complex(x(:,1),x(:,2)));
                        c3{4}=solofun(c3{4},@(x)complex(x(:,1),x(:,2)));
                    end
                    
                    % apply Fisher's transform
                    if(opt.ZTRANS); c3{3}=solofun(c3{3},@fisher); end
                    if(opt.ZTRANS); c3{4}=solofun(c3{4},@fisher); end
                    
                    % multiply by scale
                    c3{3}=multiply(c3{3},scale(inb1).*scale(inb2));
                    c3{4}=multiply(c3{4},scale(inb1).*scale(inb2));
                    
                    % twice as many coda sets
                    nsets=4;
            end
            
            % stack
            cnt=cnt+1; % increment c3 counter
            sdata(cnt,1)=addrecords(c3{:});
            sscale(cnt)=sum(nsets*scale(inb1).*scale(inb2));
            ncoda(cnt)=nin;
            
            % time range of input data
            [mint0,mint0]=min(timediff(t0utc(inb1(1),:),t0utc(inb1,:),'utc'));
            [maxt1,maxt1]=max(timediff(t1utc(inb1(1),:),t1utc(inb1,:),'utc'));
            sautc{cnt}=t0utc(inb1(mint0),:);
            sfutc{cnt}=t1utc(inb1(maxt1),:);
            [mint0,mint0]=min(timediff(t0utc(inb2(1),:),t0utc(inb2,:),'utc'));
            [maxt1,maxt1]=max(timediff(t1utc(inb2(1),:),t1utc(inb2,:),'utc'));
            st0utc{cnt}=t0utc(inb2(mint0),:);
            st1utc{cnt}=t1utc(inb2(maxt1),:);
            
            % for debugging
            sdata(cnt).misc.c3master=strcat({coda{1}(inb1).path}',...
                {coda{1}(inb1).name}');
            sdata(cnt).misc.c3slave=strcat({coda{1}(inb2).path}',...
                {coda{1}(inb2).name}');

            % detail message
            if(verbose)
                thispair=(s1-1)*nstns-sum(1:s1-1)+s2-s1*opt.NOAUTO;
                print_time_left(thispair,npairs);
            end
            %if(cnt==3); break; end
        end
        %if(cnt==3); break; end
    end
    % detail message
    if(verbose && thispair~=npairs); print_time_left(npairs,npairs); end
    
    % any coda correlations?
    if(cnt==0); sdata=[]; return; end
    
    % divide by scale to get back to an average
    % - updates dep* stats skipped by addrecords hack
    sdata=divide(sdata,sscale);
    
    % unapply Fisher's transform
    if(opt.C3ZTRANS); sdata=solofun(sdata,@ifisher); end
    
    % convert cplx to fd
    % - updates dep* to not be complex
    if(~isempty(opt.FDOUT))
        sdata=solofun(sdata,@(x)[real(x),imag(x)]);
    end
    
    % update headers
    sdata=changeheader(sdata,'scale',sscale,'resp0',opt.C3VRAYL,...
        'resp1',opt.C3WINOFF,'resp2',opt.C3WINLEN,'resp3',ncoda,...
        'kuser2',opt.C3XCCMP,'a utc',sautc,'f utc',sfutc,...
        't0 utc',st0utc,'t1 utc',st1utc,'t3 utc',sautc,'t4 utc',sfutc,...
        'iztype','iunkn');

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % toggle verbosity back
    seizmoverbose(verbose);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
    
    % toggle verbosity back
    seizmoverbose(verbose);
    
    % rethrow error
    error(lasterror);
end

end


function [opt]=noise_c3_parameters(varargin)
% parses/checks noise_c3 pv pairs

% defaults
varargin=[{'c3v' 2.6 'c3wo' 300 'c3wl' 1200 'c3x' 'sym' 'co' false ...
    'c3z' true 'c3s' [] 'maxlag' 4000 'normxc' true 'fdo' false ...
    'noa' false 'minc' 1} varargin];

% get user input
for i=1:2:numel(varargin)
    switch lower(varargin{i})
        case {'v' 'vr' 'vrayl' 'c3v' 'c3vr' 'c3vrayl'}
            if(isempty(varargin{i+1})); continue; end
            opt.C3VRAYL=varargin{i+1};
        case {'wo' 'wino' 'winoff' 'c3wo' 'c3wino' 'c3winoff'}
            if(isempty(varargin{i+1})); continue; end
            opt.C3WINOFF=varargin{i+1};
        case {'wl' 'winl' 'winlen' 'c3wl' 'c3winl' 'c3winlen'}
            if(isempty(varargin{i+1})); continue; end
            opt.C3WINLEN=varargin{i+1};
        case {'x' 'xc' 'xccmp' 'c3x' 'c3xc' 'c3xccmp'}
            if(isempty(varargin{i+1})); continue; end
            opt.C3XCCMP=varargin{i+1};
        case {'z' 'zt' 'ztrans' 'ztransform' ...
                'c3z' 'c3zt' 'c3ztrans' 'c3ztransform'}
            if(isempty(varargin{i+1})); continue; end
            opt.C3ZTRANS=varargin{i+1};
        case {'s' 'sl' 'stn' 'stnlist' ...
                'c3s' 'c3sl' 'c3stn' 'c3stnlist'}
            opt.C3STNLIST=varargin{i+1};
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
        case {'mincoda' 'minco' 'minc' 'minc3'}
            if(isempty(varargin{i+1})); continue; end
            opt.MINCODA=varargin{i+1};
        otherwise
            error('seizmo:noise_c3:badInput',...
                'Unknown Option: %s !',varargin{i});
    end
end

% fix string options to be cellstr vectors
if(ischar(opt.C3STNLIST)); opt.C3STNLIST=cellstr(opt.C3STNLIST); end
if(iscellstr(opt.C3STNLIST))
    opt.C3STNLIST=unique(lower(opt.C3STNLIST(:)));
end

% check options
if(~isnumeric(opt.C3VRAYL) || ~isscalar(opt.C3VRAYL) ...
        || ~isreal(opt.C3VRAYL) || opt.C3VRAYL<=0)
    error('seizmo:noise_c3:badInput',...
        'VRAYL must be a positive real scalar!');
elseif(~isnumeric(opt.C3WINOFF) || ~isscalar(opt.C3WINOFF) ...
        || ~isreal(opt.C3WINOFF))
    error('seizmo:noise_c3:badInput',...
        'WINOFF must be a real-valued scalar!');
elseif(~isnumeric(opt.C3WINLEN) || ~isscalar(opt.C3WINLEN) ...
        || ~isreal(opt.C3WINLEN) || opt.C3WINLEN<=0)
    error('seizmo:noise_c3:badInput',...
        'WINLEN must be a positive real scalar!');
elseif(~ischar(opt.C3XCCMP) || ~isvector(opt.C3XCCMP) ...
        || size(opt.C3XCCMP,1)~=1 || ~ismember(lower(opt.C3XCCMP),...
        {'causal' 'acausal' 'both' 'xboth' 'sym'}))
    error('seizmo:noise_c3:badInput',...
        'XCCMP must be ''CAUSAL'' ''ACAUSAL'' ''BOTH'' or ''SYM''!');
elseif(~isscalar(opt.C3ZTRANS) || ~islogical(opt.C3ZTRANS))
    error('seizmo:noise_c3:badInput',...
        'ZTRANS must be true or false!');
elseif(~isempty(opt.C3STNLIST) && (~iscellstr(opt.C3STNLIST)))
    error('seizmo:noise_c3:badInput',...
        'STNLIST must be a string list of KNAME codes!');
elseif(~isscalar(opt.XCMAXLAG) || ~isreal(opt.XCMAXLAG) ...
        || opt.XCMAXLAG<=0)
    error('seizmo:noise_c3:badInput',...
        'XCMAXLAG must be a positive real-valued scalar in seconds!');
elseif(~isscalar(opt.NORMXC) || ~islogical(opt.NORMXC))
    error('seizmo:noise_c3:badInput',...
        'NORMXC must be true or false!');
elseif(~isscalar(opt.COHERENCY) || ~islogical(opt.COHERENCY))
    error('seizmo:noise_c3:badInput',...
        'COHERENCY must be true or false!');
elseif(~isscalar(opt.FDOUT) || ~islogical(opt.FDOUT))
    error('seizmo:noise_c3:badInput',...
        'FDOUT must be true or false!');
elseif(~isscalar(opt.NOAUTO) || ~islogical(opt.NOAUTO))
    error('seizmo:noise_c3:badInput',...
        'NOAUTO must be true or false!');
elseif(~isnumeric(opt.MINCODA) || ~isscalar(opt.MINCODA) ...
        || ~isreal(opt.MINCODA) || opt.MINCODA<=0 ...
        || opt.MINCODA~=fix(opt.MINCODA))
    error('seizmo:noise_c3:badInput',...
        'MINCODA must be a positive integer>0!');
elseif(~opt.NORMXC && opt.ZTRANS && ~opt.COHERENCY)
    warning('seizmo:noise_c3:badInput',...
        'ZTRANS=TRUE not allowed when NORMXC=FALSE & COHERENCY=FALSE!');
    opt.ZTRANS=false;
end

% normxc logical to option
if(opt.NORMXC); opt.NORMXC={'normxc'}; else opt.NORMXC={}; end
if(opt.COHERENCY); opt.COHERENCY={'coherency'}; else opt.COHERENCY={}; end
if(opt.FDOUT); opt.FDOUT={'fdout'}; else opt.FDOUT={}; end

end


function [d1]=addrecords(varargin)
% simple hack for speed (no dep* update)
try
    for i=1:nargin
        varargin{i}(1).dep=mean(cat(3,varargin{i}.dep),3);
        varargin{i}(2:end)=[];
    end
    d1=varargin{1};
    d1.dep=mean(cat(3,varargin{:}.dep),3);
catch
    error('seizmo:noise_c3:badNCFs',...
        'NCFs differ in number of points! Cannot stack!');
end
end

