function [data]=correlate(master,varargin)
%CORRELATE    Compute cross correlograms of SEIZMO data records
%
%    Usage:    correlograms=correlate(master,slave)
%              correlograms=correlate(master)
%              correlograms=correlate(...,'fdout',...)
%              correlograms=correlate(...,'normxc',...)
%              correlograms=correlate(...,'coherency',...)
%              correlograms=correlate(...,'mcxc',...)
%              correlograms=correlate(...,'noauto',...)
%              correlograms=correlate(...,'reltime',...)
%              correlograms=correlate(...,[lagmin lagmax],...)
%              correlograms=correlate(...,'absxc',...)
%              peaks=correlate(...,'peaks',{'param1',val1,...},...)
%
%    Description:
%     CORRELOGRAMS=CORRELATE(MASTER,SLAVE) cross correlates the pairs given
%     by the MASTER and SLAVE datasets on a index by index basis.  That
%     means that the record MASTER(3) is correlated with SLAVE(3).  If
%     either MASTER or SLAVE is scalar it is applied to all the records of
%     the other.  Please note that the correlation operation is done in the
%     frequency-domain but the results are returned in the time domain
%     (unless FDOUT is set to true -- see below).  The output CORRELOGRAMS
%     will match the size of the non-scalar input (if both are scalar then
%     CORRELOGRAMS is scalar).
%
%     CORRELOGRAMS=CORRELATE(MASTER) returns the autocorrelations for each
%     record in MASTER.  The output CORRELOGRAMS is the same size as the
%     input.
%
%     CORRELOGRAMS=CORRELATE(...,'FDOUT',...) outputs the correlograms in
%     the frequency domain.  The records are in 'RLIM' format (see the
%     function DFT for more info).  Note that the spectral amplitudes do
%     NOT include the DELTA factor included by DFT so converting back to
%     the time-domain using IDFT requires multiplication by the time-domain
%     DELTA as well as a few other adjustments - see the example in the
%     function DEPHASE_CORRELATIONS for details!  Please note that regular
%     correlograms and coherency correlograms have a different number of
%     spectral points as coherency correlograms need more resolution to
%     avoid aliasing issues.  The 'FDOUT' option is ignored if the 'PEAKS'
%     option is set!  This option also causes the 'ABSXC' & [LAGMIN LAGMAX]
%     options to be ignored unless the 'PEAKS' option is set!
%
%     CORRELOGRAMS=CORRELATE(...,'NORMXC',...) outputs correlograms that
%     are normalized by the zero lag value of the autocorrelations.  The
%     values for a normalized correlogram are in the range from -1 to 1,
%     with 1 being a perfect correlation between the two records and -1
%     a perfect anticorrelation.  This requires slightly more computation
%     time as the zero-lag values of the autocorrelations must be computed.
%
%     CORRELOGRAMS=CORRELATE(...,'COHERENCY',...) outputs coherency
%     correlograms that are normalized frequency by frequency using the
%     autocorrelation frequency values.  The result is frequency values in
%     the range from -1 to 1, with 1 being a perfectly coherent correlation
%     between the two records and -1 a perfectly coherent anticorrelation.
%     Please be aware that converting the coherency correlograms to the
%     time-domain can affect the frequency values due to windowing of the
%     data.  This can be avoided by setting FDOUT to TRUE.  The 'COHERENCY'
%     option requires more computation time as autocorrelations must be
%     computed.  The NORMXC option is ignored when this option is set.
%
%     CORRELOGRAMS=CORRELATE(...,'MCXC',...) cross correlates all possible
%     pairings between the MASTER & SLAVE datasets (or if only the MASTER
%     set is given this computes all the unique cross correlations as well
%     as autocorrelations).  This helps to reduce computation time and
%     memory by making use of redundancies in the computations.
%
%     CORRELOGRAMS=CORRELATE(...,'NOAUTO',...) skips autocorrelations for
%     the MASTER only 'MCXC' case.
%
%     CORRELOGRAMS=CORRELATE(...,'RELTIME',...) makes the lags based on the
%     relative timing of the records.  That is, using this flag will cause
%     CORRELATE to ignore any differences in the reference timing of the
%     input records.  By default the lags are based on the absolute timing
%     of the correlated records.
%
%     CORRELOGRAMS=CORRELATE(...,[LAGMIN LAGMAX],...) limits the output
%     correlograms to the lag points within the specified range.  Note that
%     if the range extends past the default range (all possible non-zero
%     points), then the correlograms are padded with zeros as needed.  The
%     default (all possible non-zero points) may be specified with an empty
%     matrix (ie. []).  Also you may specify a single value to get a
%     symmetric range.  This option is ignored if FDOUT is set to true and
%     PEAKS is not set!
%
%     CORRELOGRAMS=CORRELATE(...,'ABSXC',...) takes the absolute value of
%     the correlograms.  This is only useful for peak picking (see the next
%     option) which forces the peak picker to look at troughs in the
%     negative correlation range.  Note that troughs in the positive
%     correlation range will not become peaks by this action.  This option
%     is ignored if FDOUT is set to true and PEAKS is not set!
%
%     PEAKS=CORRELATE(...,'PEAKS',{'PARAM1',VAL1,...},...) passes the
%     correlograms to a peak picker.  The output struct PEAKS contains the
%     fields .cg, .lg, & .pg of size NPAIRSx1xNPEAKSxNADJACENT where NPEAKS
%     & NADJACENT are defined by the peaks options in the cell array input.
%     There are two additional fields: .m & .s which give the indices of
%     the correlated records for each pair as NPAIRSx1 vectors.  See the
%     function GETPEAKS for more info on paramter/value pairs that may be
%     given.  The .pg field indicates the polarity of the peak (always 1
%     unless the 'ABSXC' parameter is passed in which case the elements are
%     either 1 or -1).  This option overrides the FDOUT option!
%
%    Notes:
%     - Correlated records are required to have a common sample rate (DELTA
%       field must be equal), be evenly sampled (LEVEN field must be TRUE),
%       and single component (NCMP field must be 1).  All records are
%       also passed through CHECKHEADER so sane settings for all header
%       fields are enforced.
%     - Frequency domain records are okay for input but require a few extra
%       constraints that can be tricky!  The first is that the frequency
%       domain sample rate (header field DELTA) is equal for these records.
%       Not so bad, right?  Well then NPTS (the number of frequency domain
%       points) must also be equal and satisfies the constraint
%       NPTS>=2^NEXTPOW2(2*MAX(NSPTS)-1) where NSPTS is the number of time
%       domain points.  This condition is necessary to get the correlation
%       values and is difficult to satify.  Typically, if you are lucky and
%       all your records have equal number of time domain points, this
%       requires setting POW2PAD=1 or, if the number of points varies
%       between records, setting POW2PAD=1i*2^NEXTPOW2(2*MAX(NPTS)-1) in
%       DFT when making the frequency domain data will be enough.  For
%       coherency correlograms it is best to triple this requirement.
%     - The correlograms are given filenames using the following format:
%       CORR_-_MASTER_-_REC<idx>_-_<kname>_-_SLAVE_-_REC<idx>_-_<kname>
%       where <idx> is the index of the record in MASTER/SLAVE and <kname>
%       is the fields knetwk, kstnm, khole, kcmpnm of record <idx> joined
%       with periods ('.') in between.  The path is set to the current
%       directory ('.'), while byte-order uses that which is native to the
%       current system.  Filetype is SAC v6 binary file.  See the Header
%       changes section for details on info retained in the header.
%     - The frequency-domain and time-domain outputs are typically slightly
%       different due to numerical noise.  This is exacerbated by coherency
%       correlograms which involve a deconvolution via spectral division.
%       In that case, values outside the maximum lag range are large and
%       their removal during time-domain conversion causes departure from
%       the frequency-domain representation.  For coherency correlograms it
%       is best to stay in the frequency-domain to avoid this issue.
%
%    Header Changes:
%     DEPMEN, DEPMIN, DEPMAX, NPTS
%     Z is the reference time of the master record.
%     A, F is the absolute time limits (B, E) of the master record.
%     T0, T1 is the absolute time limits (B, E) of the slave record.
%     T2=1 if lags are in absolute time, T2=0 if lags are in relative time.
%     T3, T4 is the relative time limits (B, E) of the slave record.
%     NXSIZE, NYSIZE are the original NPTS of the master, slave record.
%     B, E give the lag range.
%     SB, SDELTA, NSPTS are from the frequency-domain multiplication.
%     USER0 is the index of master record & KUSER0 is 'MASTER'.
%     USER1 is the index of slave record & KUSER1 is 'SLAVE'.
%
%     The following info is also retained:
%      SLAVE RECORD FIELD   CORRELOGRAM FIELD
%       STLA                 STLA
%       STLO                 STLO
%       STEL                 STEL
%       STDP                 STDP
%       KNETWK               KNETWK
%       KSTNM                KSTNM
%       KHOLE                KHOLE
%       KCMPNM               KCMPNM
%       CMPINC               CMPINC
%       CMPAZ                CMPAZ
%      MASTER RECORD FIELD  CORRELOGRAM FIELD
%       STLA                 EVLA
%       STLO                 EVLO
%       STEL                 EVEL
%       STDP                 EVDP
%       KNETWK               KT0
%       KSTNM                KT1
%       KHOLE                KT2
%       KCMPNM               KT3
%       CMPINC               USER2
%       CMPAZ                USER3
%
%    Examples:
%     % Roughly equivalent to 'correlate' in SAC:
%     xc=correlate(data(1),data);
%
%     % Perform a multi-channel cross correlation of several records where
%     % correlations are normalized, no autocorrelations are returned, and
%     % the lags are based on the relative timing of each record:
%     xc=correlate(data,'mcxc','normxc','noauto','reltime');
%
%     % The same as above but only returning info about the 3 highest peaks
%     % (in an absolute value sense) for every correlogram:
%     xc=correlate(data,'m','nor','noa','r','a','p',{'n',3});
%
%    See also: CONVOLVE, MCXC, GETPEAKS, REVERSE_CORRELATIONS, ISXC,
%              ROTATE_CORRELATIONS, SPLIT_AUTO_CORRELATIONS,
%              NAME_CORRELATIONS, NO_REDUNDANT_CORRELATIONS,
%              HORZ_CORRELATIONS_SETS, IS_FULL_MATRIX_OF_CORRELATIONS,
%              DEPHASE_CORRELATIONS, WRAP_CORRELATIONS, UNWRAP_CORRELATIONS

%     Version History:
%        June 27, 2009 - first fully functional version
%        Oct.  7, 2009 - slave records NZ* info passed on, fixed record
%                        names in 2 dataset correlogram case, fixed ST/EV
%                        info in 1 dataset correlogram case
%        Dec.  5, 2009 - calculates delaz stuff
%        Dec. 13, 2009 - minor doc update
%        Jan. 27, 2010 - proper SEIZMO handling, fixed bug where delaz info
%                        does not get calculated, seizmoverbose support
%        Jan. 29, 2010 - cleaned up cross_check_data subfunction
%        Feb. 11, 2011 - mass seizmocheck fix
%        Dec.  1, 2011 - doc update, save CMPINC/CMPAZ info, warn if
%                        reftimes vary, remove cross check via better
%                        checkheader usage
%        Jan. 24, 2012 - checkheader post-correlation bugfix
%        Sep. 28, 2012 - force delta & leven of output
%        Oct. 19, 2012 - complete rewrite
%        Oct. 21, 2012 - more fixes related to testing
%        Jan. 28, 2013 - doc update
%        Jan. 30, 2013 - peaks output dimension 2 is forced scalar for
%                        compatibility with peaks output of old version,
%                        allow empty matrix [] to specify lagrng default
%                        and single element lagrng for a symmetric range
%        Feb. 27, 2013 - MAJOR BUGFIX: lagrng usage offset lags by 1 sample
%        Mar. 28, 2013 - checkheader post-correlation bugfix (forgot to add
%                        this to update the delaz fields)
%        Aug. 16, 2013 - bugfixes for no peak detection
%        Sep. 20, 2013 - updated See also section
%        June  4, 2014 - bugfix: preserve checkheader state
%        June 10, 2014 - bugfix: dataless class handling no longer errors,
%                        added coherency option (ignores normxc option),
%                        added fdout option (peaks option overrides this,
%                        lagrng & absxc are overridden otherwise), npts1/2
%                        are now output as nxsize/nysize, sb/sdelta,nspts
%                        also now output
%        June 11, 2014 - fd input allowed now, bugfix: seizmocheck state
%                        restored, bugfix: delta allowed to vary per
%                        pairing for non-mcxc case
%        June 12, 2014 - allow npts of spectral records to not be a power
%                        of 2 (although this is not allowed by checkheader)
%        June 16, 2014 - bugfix: nspts of spectral records set to # of
%                        spectra points so IDFT does not truncate data
%        June 23, 2014 - updated See also list
%        June 24, 2014 - fd output is now dephased by default
%        June 25, 2014 - bugfix: t2 header field now set, t3/t4 are b/e
%                        slave field relative times (from slave reftime)
%        July  8, 2014 - turned on 2x length multiplier for coherency
%        July 11, 2014 - changed 2x to 3x multiplier for coherency
%        July 21, 2014 - bugfix: lagrng improperly grabbed window offset by
%                        -1 sample (timing was valid)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 21, 2014 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));

% check data structure
error(seizmocheck(master,'dep'));

% verbosity
verbose=seizmoverbose;

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt header check
try
    % 1 or 2 datasets
    nrecs1=numel(master);
    if(nargin>1 && isseizmo(varargin{1},'dep'))
        slave=varargin{1};
        varargin(1)=[];
        nrecs2=numel(slave); 
        twodata=true;
        
        % checks
        master=checkheader(master,...
            'FALSE_LEVEN','ERROR',...
            'XYZ_IFTYPE','ERROR',...
            'MULCMP_DEP','ERROR');
        slave=checkheader(slave,...
            'FALSE_LEVEN','ERROR',...
            'XYZ_IFTYPE','ERROR',...
            'MULCMP_DEP','ERROR');
    else
        twodata=false;
        master=checkheader(master,...
            'FALSE_LEVEN','ERROR',...
            'XYZ_IFTYPE','ERROR',...
            'MULCMP_DEP','ERROR');
    end
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

% get options
peaksopt={};
cp=cellfun('isclass',varargin,'cell');
if(any(cp)); peaksopt=varargin{find(cp,1,'last')}; end
[pn,ps,pa,ph,pf]=check_getpeaks_options(peaksopt{:});
varargin(cp)=[];
lagrng=[];
np=cellfun('isclass',varargin,'double');
if(any(np));
    lagrng=varargin{find(np,1,'last')};
    if(isempty(lagrng))
        % do nothing
    elseif(isscalar(lagrng))
        if(~isnumeric(lagrng) || ~isreal(lagrng))
            error('seizmo:correlate:badInput',...
                'LAGRNG must be real-valued!');
        end
        lagrng=abs(lagrng)*[-1 1];
    elseif(~isnumeric(lagrng) || ~isreal(lagrng) ...
            || numel(lagrng)~=2 || diff(lagrng)<0)
        error('seizmo:correlate:badInput',...
            'LAGRNG must be a real-valued vector as [LAGMIN LAGMAX]!');
    end
end
varargin(np)=[];
if(~iscellstr(varargin))
    error('seizmo:correlate:badInput',...
        'Parameter(s) of wrong type given!');
end
[mcxc,normxc,noauto,reltime,absxc,peaks,coherency,fdout]=deal(false);
valid=['mcxc     ';'normxc   ';'noauto   ';'reltime  ';...
    'absxc    ';'peaks    ';'coherency';'fdout    '];
for i=1:numel(varargin)
    switch strmatch(lower(varargin{i}),valid)
        case 1 % mcxc
            mcxc=true;
        case 2 % normxc
            normxc=true;
        case 3 % noauto
            noauto=true;
        case 4 % reltime
            reltime=true;
        case 5 % absxc
            absxc=true;
        case 6 % peaks
            peaks=true;
        case 7 % coherency
            coherency=true;
        case 8 % fdout
            fdout=true;
        otherwise
            error('seizmo:correlate:badInput',...
                'Unknown option: %s !',varargin{i});
    end
end

% check for option conflicts
if(fdout && peaks)
    warning('seizmo:correlate:inputConflict',...
        'FDOUT option ignored because PEAKS is set!');
elseif(fdout && ~isempty(lagrng))
    warning('seizmo:correlate:inputConflict',...
        '[LAGMIN LAGMAX] option ignored because FDOUT=TRUE!');
elseif(fdout && absxc)
    warning('seizmo:correlate:inputConflict',...
        'ABSXC option ignored because FDOUT=TRUE!');
elseif(coherency && normxc)
    % turn off normxc for efficiency
    normxc=false;
end

% coherency needs to be LENMULxcorrelogram length in the
% time-domain to recover the time-domain values accurately
% - testing found that LENMUL=2 avoided wrapping significant energy
%   but is still different from LENMUL=3+ so we use 3
if(coherency); lenmul=3; else lenmul=1; end

% initial checks & adjustments
if(twodata)
    % necessary header info
    [iftype1,sdelta1,nspts1,delta1,npts1,z1,b01,e01,b1,e1,st1,kname1,...
        cmp1,sb01,sb1]=getheader(master,'iftype id','sdelta','nspts',...
        'delta','npts','z','b','e','b utc','e utc','st','kname','cmp',...
        'sb','sb utc');
    [iftype2,sdelta2,nspts2,delta2,npts2,b02,e02,b2,e2,st2,kname2,cmp2,...
        sb02,sb2]=getheader(slave,'iftype id','sdelta','nspts','delta',...
        'npts','b','e','b utc','e utc','st','kname','cmp','sb','sb utc');
    
    % who are the spectral records
    rlim1=strcmpi(iftype1,'irlim');
    amph1=strcmpi(iftype1,'iamph');
    rlim2=strcmpi(iftype2,'irlim');
    amph2=strcmpi(iftype2,'iamph');
    sp1=rlim1 | amph1;
    sp2=rlim2 | amph2;
    issp1=any(sp1);
    issp2=any(sp2);
    
    % switch npts/delta/b/butc/e for spectral records
    if(issp1)
        [sdelta1(sp1),delta1(sp1)]=deal(delta1(sp1),sdelta1(sp1));
        [nspts1(sp1),npts1(sp1)]=deal(npts1(sp1),nspts1(sp1));
        [b01(sp1),b1(sp1)]=deal(sb01(sp1),sb1(sp1));
        e01(sp1)=b01(sp1)+delta1(sp1).*(npts(sp1)-1);
        e1(sp1)=mat2cell(fixtimes(cell2mat(b1(sp1))...
            +[zeros(sum(sp1),4) (npts1(sp1)-1).*delta1(sp1)],'utc'),...
            ones(sum(sp1),1));
    end
    if(issp2)
        [sdelta2(sp2),delta2(sp2)]=deal(delta2(sp2),sdelta2(sp2));
        [nspts2(sp2),npts2(sp2)]=deal(npts2(sp2),nspts2(sp2));
        [b02(sp2),b2(sp2)]=deal(sb02(sp2),sb2(sp2));
        e02(sp2)=b02(sp2)+delta2(sp2).*(npts2(sp2)-1);
        e2(sp2)=mat2cell(fixtimes(cell2mat(b2(sp2))...
            +[zeros(sum(sp2),4) (npts2(sp2)-1).*delta2(sp2)],'utc'),...
            ones(sum(sp2),1));
    end
    
    % checking delta/npts
    if(mcxc)
        % all delta must be equal
        if(numel(unique([delta1; delta2]))~=1)
            error('seizmo:correlate:deltaNotEqual',...
                'Time sample spacing must be equal for all records!');
        end
        
        % minimum number of spectral points as a power of 2
        % - this may change to a higher number for spectral record input
        nspts=2^nextpow2(lenmul*(max(npts1)+max(npts2)-1));
        
        % check nspts/sdelta are equal across spectral records
        if(issp1 || issp2)
            % check nspts/sdelta are equal across spectral records
            if(numel(unique([sdelta1(sp1); sdelta2(sp2)]))~=1)
                error('seizmo:correlate:sdeltaNotEqual',...
                    'DELTA must be equal for all spectral records!');
            elseif(numel(unique([nspts1(sp1); nspts2(sp2)]))~=1)
                error('seizmo:correlate:nsptsNotEqual',...
                    'NPTS must be equal for all spectral records!');
            end
            
            % is nspts1/2 >= nspts?
            nspts_in=unique([nspts1(sp1); nspts2(sp2)]);
            if(max(npts1)+max(npts2)-1>nspts_in)
                error('seizmo:correlate:notEnoughFrequencies',...
                    ['NPTS of spectral records not high enough to ' ...
                    'correlate!  See the Notes section of CORRELATE.']);
            end
            
            % adjust nspts to match nspts1/2
            nspts=nspts_in;
        end
    else
        % all delta must be equal at pair level
        if(~all(delta1==delta2))
            error('seizmo:correlate:deltaNotEqual',...
                'Time sample spacing must be equal for paired records!');
        end
        
        % minimum number of spectral points as a power of 2
        % - this may change to a higher number for spectral record input
        nspts=2.^nextpow2n(lenmul*(npts1+npts2-1));
        nspts_min=npts1+npts2-1;
        
        % check nspts/sdelta are equal across spectral record pairs
        if(issp1 || issp2)
            % handle scalar input
            if(nrecs1==1)
                sdelta1=sdelta1(ones(nrecs2,1));
                nspts1=nspts1(ones(nrecs2,1));
                sp1=sp1(ones(nrecs2,1));
            end
            if(nrecs2==1)
                sdelta2=sdelta2(ones(nrecs1,1));
                nspts2=nspts2(ones(nrecs1,1));
                sp2=sp2(ones(nrecs1,1));
            end
            
            % check nspts/sdelta are equal across spectral record pairs
            bothsp=sp1 & sp2;
            if(~isequal(sdelta1(bothsp),sdelta2(bothsp)))
                error('seizmo:correlate:sdeltaNotEqual',...
                    'DELTA must be equal for paired spectral records!');
            elseif(~isequal(nspts1(bothsp),nspts2(bothsp)))
                error('seizmo:correlate:nsptsNotEqual',...
                    'NPTS must be equal for paired spectral records!');
            end
            
            % is nspts1/2 >= nspts?
            if(any(nspts_min(sp1)>nspts1(sp1)) ...
                    || any(nspts_min(sp2)>nspts2(sp2)))
                error('seizmo:correlate:notEnoughFrequencies',...
                    ['NPTS of spectral records not high enough to ' ...
                    'correlate!  See the Notes section of CORRELATE.']);
            end
            
            % adjust nspts to match nspts1/2
            nspts(sp1)=nspts1(sp1);
            nspts(sp2)=nspts2(sp2);
        end
    end
    
    % convert spectral data to complex
    % - also fixes time domain delta factor used for SAC compatibility
    for i=1:nrecs1
        if(rlim1(i))
            master(i).dep=...
                complex(master(i).dep(:,1),master(i).dep(:,2))/delta1(i);
        elseif(amph1(i))
            master(i).dep=...
                master(i).dep(:,1).*exp(1j*master(i).dep(:,2))/delta1(i);
        end
    end
    for i=1:nrecs2
        if(rlim2(i))
            slave(i).dep=...
                complex(slave(i).dep(:,1),slave(i).dep(:,2))/delta2(i);
        elseif(amph2(i))
            slave(i).dep=...
                slave(i).dep(:,1).*exp(1j*slave(i).dep(:,2))/delta2(i);
        end
    end
else % single dataset
    % necessary header info
    [iftype1,sdelta1,nspts1,delta1,npts1,z1,b01,e01,b1,e1,st1,kname1,...
        cmp1,sb01,sb1]=getheader(master,'iftype id','sdelta','nspts',...
        'delta','npts','z','b','e','b utc','e utc','st','kname','cmp',...
        'sb','sb utc');
    
    % who are the spectral records
    rlim1=strcmpi(iftype1,'irlim');
    amph1=strcmpi(iftype1,'iamph');
    sp1=rlim1 | amph1;
    issp1=any(sp1);
    
    % switch npts/delta for spectral records
    if(issp1)
        [sdelta1(sp1),delta1(sp1)]=deal(delta1(sp1),sdelta1(sp1));
        [nspts1(sp1),npts1(sp1)]=deal(npts1(sp1),nspts1(sp1));
        [b01(sp1),b1(sp1)]=deal(sb01(sp1),sb1(sp1));
        e01(sp1)=b01(sp1)+delta1(sp1).*(npts(sp1)-1);
        e1(sp1)=mat2cell(fixtimes(cell2mat(b1(sp1))...
            +[zeros(sum(sp1),4) (npts1(sp1)-1).*delta1(sp1)],'utc'),...
            ones(sum(sp1),1));
    end
    
    % checking delta/npts
    if(mcxc)
        % all delta must be equal
        if(numel(unique(delta1))~=1)
            error('seizmo:correlate:deltaNotEqual',...
                'Time sample spacing must be equal for all records!');
        end
        
        % minimum number of spectral points as a power of 2
        % - this may change to a higher number for spectral record input
        nspts=2^nextpow2(lenmul*(2*max(npts1)-1));
        
        % check nspts/sdelta are equal across spectral records
        if(issp1)
            % check nspts/sdelta are equal across spectral records
            if(numel(unique(sdelta1(sp1)))~=1)
                error('seizmo:correlate:sdeltaNotEqual',...
                    'DELTA must be equal for all spectral records!');
            elseif(numel(unique(nspts1(sp1)))~=1)
                error('seizmo:correlate:nsptsNotEqual',...
                    'NPTS must be equal for all spectral records!');
            end
            
            % is nspts1/2 >= nspts?
            nspts_in=unique(nspts1(sp1));
            if(2*max(npts1)-1>nspts_in)
                error('seizmo:correlate:notEnoughFrequencies',...
                    ['NPTS of spectral records not high enough to ' ...
                    'correlate!  See the Notes section of CORRELATE.']);
            end
            
            % adjust nspts to match nspts1/2
            nspts=nspts_in;
        end
    else
        % minimum number of spectral points as a power of 2
        % - this may change to a higher number for spectral record input
        nspts=2.^nextpow2n(lenmul*(2*npts1-1));
        
        % check nspts/sdelta are equal across spectral record pairs
        if(issp1)
            % is nspts1/2 >= nspts?
            if(any(2*npts(sp1)-1>nspts1(sp1)))
                error('seizmo:correlate:notEnoughFrequencies',...
                    ['NPTS of spectral records not high enough to ' ...
                    'correlate!  See the Notes section of CORRELATE.']);
            end
            
            % adjust nspts to match nspts1/2
            nspts(sp1)=nspts1(sp1);
        end
    end
    
    % convert spectral data to complex
    % - also fixes time domain delta factor used for SAC compatibility
    for i=1:nrecs1
        if(rlim1(i))
            master(i).dep=...
                complex(master(i).dep(:,1),master(i).dep(:,2))/delta1(i);
        elseif(amph1(i))
            master(i).dep=...
                master(i).dep(:,1).*exp(1j*master(i).dep(:,2))/delta1(i);
        end
    end
end

% getting correlation indices
if(twodata)
    if(mcxc) % all
        [s,m]=find(true(nrecs2,nrecs1));
        s=s(:);
        m=m(:);
    else % only pairs
        m=(1:nrecs1)';
        s=(1:nrecs2)';
        if(nrecs1==1); m=ones(nrecs2,1); nrecs1=nrecs2; end
        if(nrecs2==1); s=ones(nrecs1,1); nrecs2=nrecs1; end
        if(nrecs1~=nrecs2)
            error('seizmo:correlate:badInput',...
                'MASTER & SLAVE datasets must be the same size!');
        end
    end
else
    % correlate against each other or just self?
    if(mcxc) % each other
        if(noauto) % only lower
            [s,m]=find(tril(true(nrecs1),-1));
        else % lower & diag
            [s,m]=find(tril(true(nrecs1)));
        end
        s=s(:);
        m=m(:);
    else % self
        [m,s]=deal((1:nrecs1)');
    end
end

% number of correlations
npairs=numel(m);

% expand delta/nspts as they may vary with each pair
delta=delta1(m,1);
if(mcxc); nspts=nspts(ones(npairs,1),1); end

% get 1st lag time & npts of correlogram
if(twodata)
    b0=-delta.*(npts2(s,1)-1);
    if(reltime); b0=b02(s,1)-b01(m,1)+b0;
    else b0=timediff(cell2mat(b1(m,1)),cell2mat(b2(s,1)),'utc')+b0;
    end
    npts0=npts1(m,1)+npts2(s,1)-1;
else
    b0=-delta.*(npts1(s,1)-1);
    if(reltime); b0=b01(s,1)-b01(m,1)+b0;
    else b0=timediff(cell2mat(b1(m,1)),cell2mat(b1(s,1)),'utc')+b0;
    end
    npts0=npts1(m,1)+npts1(s,1)-1;
end

% .name
if(~peaks)
    maxm=num2str(fix(log10(max(m))+1));
    maxs=num2str(fix(log10(max(s))+1));
    if(twodata)
        name=strcat('CORR_-_MASTER_-_REC',...
            num2str(m,['%0' maxm 'd']),'_-_',...
            kname1(m,1),'.',kname1(m,2),'.',...
            kname1(m,3),'.',kname1(m,4),...
            '_-_SLAVE_-_REC',...
            num2str(s,['%0' maxs 'd']),'_-_',...
            kname2(s,1),'.',kname2(s,2),'.',...
            kname2(s,3),'.',kname2(s,4));
    else
        name=strcat('CORR_-_MASTER_-_REC',...
            num2str(m,['%0' maxm 'd']),'_-_',...
            kname1(m,1),'.',kname1(m,2),'.',...
            kname1(m,3),'.',kname1(m,4),...
            '_-_SLAVE_-_REC',...
            num2str(s,['%0' maxs 'd']),'_-_',...
            kname1(s,1),'.',kname1(s,2),'.',...
            kname1(s,3),'.',kname1(s,4));
    end
end

% quietly allocate output
if(peaks)
    p.cg=[]; p.lg=[]; p.pg=[];
    p=p(ones(npairs,1),1);
else
    seizmoverbose(false);
    data=bseizmo([]);
    seizmoverbose(verbose);
    data=data(ones(npairs,1));
end

% preallocate
[depmin,depmen,depmax]=deal(nan(npairs,1));
if(twodata)
    [zlac1]=deal(nan(nrecs1,1));
    [zlac2]=deal(nan(nrecs2,1));
    cmplx1=false(nrecs1,1);
    cmplx2=false(nrecs2,1);
    oclass1=cell(nrecs1,1);
else
    [zlac1]=deal(nan(nrecs1,1));
    cmplx1=false(nrecs1,1);
    oclass1=cell(nrecs1,1);
end

% detail message
if(verbose)
    disp('Correlating Records');
    print_time_left(0,npairs);
end

% loop over pairs
for i=1:npairs
    % skip dataless
    if((twodata && (npts1(m(i))==0 || npts2(s(i))==0)) ...
            || (~twodata && (npts1(m(i))==0 || npts1(s(i))==0)))
        if(peaks)
            [p(i).cg,p(i).lg]=getpeaks(zeros(0,1),peaksopt{:});
            p(i).pg=nan(size(lg));
        else
            oclass1{m(i)}=class(master(m(i)).dep);
            oclass=str2func(oclass1{m(i)});
            data(i).dep=oclass(zeros(0,1));
        end
        
        % shortcut
        if(verbose); print_time_left(i,npairs); end
        continue;
    end
    
    % convert to fd cmplx, get zlac if normalizing
    if(~cmplx1(m(i)))
        oclass1{m(i)}=class(master(m(i)).dep);
        master(m(i)).dep=double(master(m(i)).dep);
        if(normxc)
            if(sp1(m(i))) % fd zlac
                % based on Parseval's Theorem
                zlac1(m(i))=sqrt(sum(master(m(i)).dep...
                    .*conj(master(m(i)).dep))/npts1(m(i)));
            else % td zlac
                zlac1(m(i))=sqrt(sum(master(m(i)).dep.^2));
            end
            if(zlac1(m(i))==0); zlac1(m(i))=eps; end
        end
        if(~sp1(m(i)))
            master(m(i)).dep=fft(master(m(i)).dep,nspts(i),1);
        end
        cmplx1(m(i))=true;
    end
    if(twodata)
        if(~cmplx2(s(i)))
            slave(s(i)).dep=double(slave(s(i)).dep);
            if(normxc)
                if(sp2(s(i))) % fd zlac
                    % from Parseval's Theorem
                    zlac2(s(i))=sqrt(sum(slave(s(i)).dep...
                        .*conj(slave(s(i)).dep))/npts2(s(i)));
                else % td zlac
                    zlac2(s(i))=sqrt(sum(slave(s(i)).dep.^2));
                end
                if(zlac2(s(i))==0); zlac2(s(i))=eps; end
            end
            if(~sp2(s(i)))
                slave(s(i)).dep=fft(slave(s(i)).dep,nspts(i),1);
            end
            cmplx2(s(i))=true;
        end
    else % single dataset
        if(~cmplx1(s(i)))
            oclass1{s(i)}=class(master(s(i)).dep);
            master(s(i)).dep=double(master(s(i)).dep);
            if(normxc)
                if(sp1(s(i))) % fd zlac
                    % from Parseval's Theorem
                    zlac1(s(i))=sqrt(sum(master(s(i)).dep...
                        .*conj(master(s(i)).dep))/npts1(s(i)));
                else % td zlac
                    zlac1(s(i))=sqrt(sum(master(s(i)).dep.^2));
                end
                if(zlac1(s(i))==0); zlac1(s(i))=eps; end
            end
            if(~sp1(s(i)))
                master(s(i)).dep=fft(master(s(i)).dep,nspts(i),1);
            end
            cmplx1(s(i))=true;
        end
    end
    
    % correlate
    if(twodata)
        if(coherency)
            tmp=conj(master(m(i)).dep).*slave(s(i)).dep...
                ./sqrt(conj(master(m(i)).dep).*master(m(i)).dep...
                .*conj(slave(s(i)).dep).*slave(s(i)).dep);
        else
            tmp=conj(master(m(i)).dep).*slave(s(i)).dep;
            if(normxc); tmp=tmp./(zlac1(m(i))*zlac2(s(i))); end
        end
        if(peaks || ~fdout)
            tmp=ifft(tmp,[],1,'symmetric');
            tmp=tmp([nspts(i)-npts2(s(i))+2:nspts(i) 1:npts1(m(i))],1);
        end
    else % single dataset
        if(coherency)
            tmp=conj(master(m(i)).dep).*master(s(i)).dep...
                ./sqrt(conj(master(m(i)).dep).*master(m(i)).dep...
                .*conj(master(s(i)).dep).*master(s(i)).dep);
        else
            tmp=conj(master(m(i)).dep).*master(s(i)).dep;
            if(normxc); tmp=tmp./(zlac1(m(i))*zlac1(s(i))); end
        end
        if(peaks || ~fdout)
            tmp=ifft(tmp,[],1,'symmetric');
            tmp=tmp([nspts(i)-npts1(s(i))+2:nspts(i) 1:npts1(m(i))],1);
        end
    end
    
    % convert back to original class
    oclass=str2func(oclass1{m(i)});
    tmp=oclass(tmp);
    
    % lagrng
    if((peaks || ~fdout) && ~isempty(lagrng))
        % ceil/floor keep points to within lag range
        lidx(1)=ceil((lagrng(1)-b0(i))/delta(i))+1;
        lidx(2)=floor((lagrng(2)-b0(i))/delta(i))+1;
        lidx=lidx(1):lidx(2);
        tmp=[zeros(sum(lidx<=0),1); ...
            tmp(lidx(lidx>0 & lidx<=npts0(i))); ...
            zeros(sum(lidx>npts0(i)),1)];
        npts0(i)=numel(tmp);
        b0(i)=b0(i)+(min(lidx)-1)*delta(i);
    end
    
    % peaks
    if(peaks) % peak output
        if(absxc)
            [p(i).cg,p(i).lg]=getpeaks(abs(tmp),pn,ps,pa,ph,pf);
            p(i).pg=nan(size(p(i).lg));
            in=p(i).lg>0 | p(i).lg<=npts0(i);
            p(i).pg(in)=sign(abs(tmp(p(i).lg(in))).*tmp(p(i).lg(in)));
        else
            [p(i).cg,p(i).lg]=getpeaks(tmp,pn,ps,pa,ph,pf);
            p(i).pg=nan(size(p(i).lg));
            p(i).pg(p(i).lg>0 | p(i).lg<=npts0(i))=1;
        end
        p(i).lg=b0(i)+delta(i)*(p(i).lg-1);
    else % correlogram output
        % assign tmp to output data struct
        if(fdout)
            data(i).dep=[real(tmp) imag(tmp)];
        else
            data(i).dep=tmp;
            
            % absolute value only applied to time-domain output
            if(absxc); data(i).dep=abs(data(i).dep); end
        end
        
        % dep*
        depmen(i)=nanmean(data(i).dep(:));
        depmin(i)=min(data(i).dep(:));
        depmax(i)=max(data(i).dep(:));
        
        % .name
        data(i).name=name{i};
    end
    
    % detail message
    if(verbose); print_time_left(i,npairs); end
end

% output for peaks
if(peaks)
    data=[];
    data.cg=permute(cat(3,p.cg),[3 4 1 2]);
    data.lg=permute(cat(3,p.lg),[3 4 1 2]);
    data.pg=permute(cat(3,p.pg),[3 4 1 2]);
    data.m=m;
    data.s=s;
    return;
end

% set header info
if(fdout)
    if(twodata)
        data=changeheader(data,'z',z1(m),...
            'sdelta',delta,'nspts',nspts,'sb',b0,'b',0,...
            'delta',1./(delta.*nspts),'npts',nspts,'e',1./(2*delta),...
            'kuser0','MASTER','user0',m,'kuser1','SLAVE','user1',s,...
            'depmin',depmin,'depmen',depmen,'depmax',depmax,...
            'a utc',b1(m),'f utc',e1(m),'t0 utc',b2(s),'t1 utc',e2(s),...
            't2',~reltime,'t3',b02(s),'t4',e02(s),...
            'ev',st1(m,:),'st',st2(s,:),'kt0',kname1(m,1),...
            'kt1',kname1(m,2),'kt2',kname1(m,3),'kt3',kname1(m,4),...
            'kname',kname2(s,:),'user2',cmp1(m,1),'user3',cmp1(m,2),...
            'cmp',cmp2(s,:),'iftype','irlim','nxsize',npts1(m),...
            'nysize',npts2(s));
    else
        data=changeheader(data,'z',z1(m),...
            'sdelta',delta,'nspts',nspts,'sb',b0,'b',0,...
            'delta',1./(delta.*nspts),'npts',nspts,'e',1./(2*delta),...
            'kuser0','MASTER','user0',m,'kuser1','SLAVE','user1',s,...
            'depmin',depmin,'depmen',depmen,'depmax',depmax,...
            'a utc',b1(m),'f utc',e1(m),'t0 utc',b1(s),'t1 utc',e1(s),...
            't2',~reltime,'t3',b01(s),'t4',e01(s),...
            'ev',st1(m,:),'st',st1(s,:),'kt0',kname1(m,1),...
            'kt1',kname1(m,2),'kt2',kname1(m,3),'kt3',kname1(m,4),...
            'kname',kname1(s,:),'user2',cmp1(m,1),'user3',cmp1(m,2),...
            'cmp',cmp1(s,:),'iftype','irlim','nxsize',npts1(m),...
            'nysize',npts1(s));
    end
else
    if(twodata)
        data=changeheader(data,'z',z1(m),...
            'delta',delta,'npts',npts0,'b',b0,'e',b0+delta.*(npts0-1),...
            'sdelta',1./(delta.*nspts),'nspts',nspts,'sb',0,...
            'kuser0','MASTER','user0',m,'kuser1','SLAVE','user1',s,...
            'depmin',depmin,'depmen',depmen,'depmax',depmax,...
            'a utc',b1(m),'f utc',e1(m),'t0 utc',b2(s),'t1 utc',e2(s),...
            't2',~reltime,'t3',b02(s),'t4',e02(s),...
            'ev',st1(m,:),'st',st2(s,:),'kt0',kname1(m,1),...
            'kt1',kname1(m,2),'kt2',kname1(m,3),'kt3',kname1(m,4),...
            'kname',kname2(s,:),'user2',cmp1(m,1),'user3',cmp1(m,2),...
            'cmp',cmp2(s,:),'nxsize',npts1(m),'nysize',npts2(s));
    else
        data=changeheader(data,'z',z1(m),...
            'delta',delta,'npts',npts0,'b',b0,'e',b0+delta.*(npts0-1),...
            'sdelta',1./(delta.*nspts),'nspts',nspts,'sb',0,...
            'kuser0','MASTER','user0',m,'kuser1','SLAVE','user1',s,...
            'depmin',depmin,'depmen',depmen,'depmax',depmax,...
            'a utc',b1(m),'f utc',e1(m),'t0 utc',b1(s),'t1 utc',e1(s),...
            't2',~reltime,'t3',b01(s),'t4',e01(s),...
            'ev',st1(m,:),'st',st1(s,:),'kt0',kname1(m,1),...
            'kt1',kname1(m,2),'kt2',kname1(m,3),'kt3',kname1(m,4),...
            'kname',kname1(s,:),'user2',cmp1(m,1),'user3',cmp1(m,2),...
            'cmp',cmp1(s,:),'nxsize',npts1(m),'nysize',npts1(s));
    end
end

% update delaz info
oldcheckheaderstate=checkheader_state(true);
data=checkheader(data,'all','ignore','old_delaz','fix');
checkheader_state(oldcheckheaderstate);

% dephase fd correlations
if(fdout)
    oldseizmocheckstate=seizmocheck_state(false);
    oldcheckheaderstate=checkheader_state(false);
    data=dephase_correlations(data);
    seizmocheck_state(oldseizmocheckstate);
    checkheader_state(oldcheckheaderstate);
end

end


function [n,s,a,h,f]=check_getpeaks_options(varargin)
n=1; s=1; a=0; h=true; f=true;
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:getpeaks:badInput',...
        'Parameters must be specified with strings!');
end
valid=['npeaks  '; 'spacing '; 'adjacent'; 'highest ';'fast    '];
for i=1:2:nargin-1
    switch strmatch(lower(varargin{i}),valid)
        case 1 % npeaks
            n=varargin{i+1};
            if(~isscalar(n) || ~isnumeric(n) || ~isreal(n) || n<=0 ...
                    || n~=fix(n))
                error('seizmo:getpeaks:badInput',...
                    'NPEAKS must be an integer >=1 !');
            end
        case 2 % spacing
            s=varargin{i+1};
            if(~isscalar(s) || ~isnumeric(s) || ~isreal(s) || s<=0 ...
                    || s~=fix(s))
                error('seizmo:getpeaks:badInput',...
                    'SPACING must be an integer >=1 !');
            end
        case 3 % adjacent
            a=varargin{i+1};
            if(~isscalar(a) || ~isnumeric(a) || ~isreal(a) || a<=0 ...
                    || a~=fix(a))
                error('seizmo:getpeaks:badInput',...
                    'ADJACENT must be an integer >=0 !');
            end
        case 4 % highest
            h=varargin{i+1};
            if(~isscalar(h) || ~islogical(h))
                error('seizmo:getpeaks:badInput',...
                    'HIGHEST must be TRUE or FALSE!');
            end
        case 5 % fast
            f=varargin{i+1};
            if(~isscalar(f) || ~islogical(f))
                error('seizmo:getpeaks:badInput',...
                    'FAST must be TRUE or FALSE!');
            end
        otherwise
            error('seizmo:getpeaks:badInput',...
                'Unknown Parameter: %s !',varargin{i});
    end
end
end


function [xp,xi]=getpeaks(x,n,s,a,h,f)
if(f && n==1)
    [xp,xi]=max(x);
else
    % force column vector
    x=x(:);
    
    % all peaks toward positive infinity
    xi=find(diff([-inf; x])>0 & diff([x; -inf])<0);
    xp=x(xi);
    
    % sort by peak height
    if(h); [xp,idx]=sort(xp,'descend');
    else [xp,idx]=sort(xp,'ascend');
    end
    xi=xi(idx);
    
    % eliminate by spacing
    if(s>1)
        peak=1;
        while(peak<min(n+1,numel(xi)))
            % eliminate those within s units of current peak
            xi([false(peak,1); abs(xi(peak)-xi(peak+1:end))<=s])=[];
            peak=peak+1;
        end
    end
    
    % pad/clip if necessary
    xi(max(end,1):n,1)=nan;
    xi(n+1:end)=[];
    xp(max(end,1):n,1)=nan;
    xp(n+1:end)=[];
end

% grab adjacent points
if(a>0)
    adj=-a:a;
    xi=xi(:,ones(1,2*a+1))+adj(ones(n,1),:);
    xp=nan(size(xi));
    xp(xi>0 & xi<=numel(x))=x(xi(xi>0 & xi<=numel(x)));
end
end

