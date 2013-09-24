function [data]=correlate(master,varargin)
%CORRELATE    Compute cross correlograms of SEIZMO data records
%
%    Usage:    correlograms=correlate(master,slave)
%              correlograms=correlate(master)
%              correlograms=correlate(...,'normxc',...)
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
%     frequency-domain but the results are returned in the time domain.
%     The output CORRELOGRAMS will match the size of the non-scalar input
%     (if both are scalar then CORRELOGRAMS is scalar).
%
%     CORRELOGRAMS=CORRELATE(MASTER) returns the autocorrelations for each
%     record in MASTER.  The output CORRELOGRAMS is the same size as the
%     input.
%
%     CORRELOGRAMS=CORRELATE(...,'NORMXC',...) outputs correlograms that
%     are normalized by the zero lag value of the autocorrelations.  The
%     values for a normalized correlogram are in the range from -1 to 1,
%     with 1 being a perfect correlation between the two records and -1
%     a perfect anticorrelation.  This requires slightly more computation
%     time as the zero-lag values of the autocorrelations must be computed.
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
%     symmetric range.
%
%     CORRELOGRAMS=CORRELATE(...,'ABSXC',...) takes the absolute value of
%     the correlograms.  This is only useful for peak picking (see the next
%     option) which forces the peak picker to look at troughs in the
%     negative correlation range.  Note that troughs in the positive
%     correlation range will not become peaks by this action.
%
%     PEAKS=CORRELATE(...,'PEAKS',{'PARAM1',VAL1,...},...) passes the
%     correlograms to a peak picker.  The output struct PEAKS contains the
%     fields .cg, .lg, & .pg of size NPAIRSx1xNPEAKSxNADJACENT where NPEAKS
%     & NADJACENT are defined by the peaks options in the cell array.
%     There are two additional fields: .m & .s which give the indices of
%     the correlated records for each pair as NPAIRSx1 vectors.  See the
%     function GETPEAKS for more info on paramter/value pairs that may be
%     given.  The .pg field indicates the polarity of the peak (always 1
%     unless the 'ABSXC' parameter is passed in which case the elements are
%     either 1 or -1).
%
%    Notes:
%     - All records are required to have a common sample rate (DELTA field
%       should be equal), be evenly sampled (LEVEN field should be TRUE),
%       and single component (NCMP field should be 1).  All records are
%       also passed through CHECKHEADER so sane settings for all header
%       fields are enforced.
%     - The correlograms are given filenames using the following format:
%       CORR_-_MASTER_-_REC<idx>_-_<kname>_-_SLAVE_-_REC<idx>_-_<kname>
%       where <idx> is the index of the record in MASTER/SLAVE and <kname>
%       is the fields knetwk, kstnm, khole, kcmpnm of record <idx> joined
%       with periods ('.') in between.  The path is set to the current
%       directory ('.'), while byte-order uses that which is native to the
%       current system.  Filetype is SAC v6 binary file.  See the Header
%       changes section for details on info retained in the header.
%
%    Header Changes:
%     DEPMEN, DEPMIN, DEPMAX, NPTS
%     Z is the reference time of the master record.
%     A, F is the time limits (B, E) of the master record.
%     T0, T1 is the time limits (B, E) of the slave record.
%     B, E give the lag range.
%     T2=1 if lags are in absolute time, T2=0 if lags are in relative time.
%     USER0 is the index of master record & KUSER0 is 'MASTER'.
%     USER1 is the index of slave record & KUSER1 is 'SLAVE'.
%
%     The following info is retained:
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
%              HORZ_CORRELATIONS_SETS, IS_FULL_MATRIX_OF_CORRELATIONS
%              

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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 20, 2013 at 15:05 GMT

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
            'NONTIME_IFTYPE','ERROR',...
            'MULTIPLE_DELTA','ERROR',...
            'MULCMP_DEP','ERROR');
        slave=checkheader(slave,...
            'FALSE_LEVEN','ERROR',...
            'NONTIME_IFTYPE','ERROR',...
            'MULTIPLE_DELTA','ERROR',...
            'MULCMP_DEP','ERROR');
    else
        twodata=false;
        master=checkheader(master,...
            'FALSE_LEVEN','ERROR',...
            'NONTIME_IFTYPE','ERROR',...
            'MULTIPLE_DELTA','ERROR',...
            'MULCMP_DEP','ERROR');
    end
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

% get sample spacing
delta=getheader(master(1),'delta');
if(twodata)
    % check delta matches between the two datasets
    if(delta~=getheader(slave(1),'delta'))
        error('seizmo:correlate:badInput',...
            'DELTA field must match between MASTER & SLAVE datasets!');
    end
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
[mcxc,normxc,noauto,reltime,absxc,peaks]=deal(false);
valid=['mcxc   ';'normxc ';'noauto ';'reltime';'absxc  ';'peaks  '];
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
        otherwise
            error('seizmo:correlate:badInput',...
                'Unknown option: %s !',varargin{i});
    end
end

% getting correlation indices
if(twodata)
    if(mcxc) % all
        [s,m]=find(true(nrecs2,nrecs1));
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
    else % self
        [m,s]=deal((1:nrecs1)');
    end
end
maxm=num2str(fix(log10(max(m))+1));
maxs=num2str(fix(log10(max(s))+1));

% number of correlations
npairs=numel(m);

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

% get header info
% + allocate other info
if(twodata)
    [z1,b01,b1,e1,npts1,st1,kname1,cmp1]=getheader(master,...
        'z','b','b utc','e utc','npts','st','kname','cmp');
    [b02,b2,e2,npts2,st2,kname2,cmp2]=getheader(slave,...
        'b','b utc','e utc','npts','st','kname','cmp');
    [zlac1]=deal(nan(nrecs1,1));
    [zlac2]=deal(nan(nrecs2,1));
    cmplx1=false(nrecs1,1);
    cmplx2=false(nrecs2,1);
    nspts=2^nextpow2(2*max([npts1; npts2])-1);
    oclass1=cell(nrecs1,1);
else
    [z1,b01,b1,e1,npts1,st1,kname1,cmp1]=getheader(master,...
        'z','b','b utc','e utc','npts','st','kname','cmp');
    [zlac1]=deal(nan(nrecs1,1));
    cmplx1=false(nrecs1,1);
    nspts=2^nextpow2(2*max(npts1)-1);
    oclass1=cell(nrecs1,1);
end

% get 1st lag time & npts of correlogram
if(twodata)
    b0=-delta*(npts2(s)-1);
    if(reltime); b0=b02(s)-b01(m)+b0;
    else b0=timediff(cell2mat(b1(m)),cell2mat(b2(s)),'utc')+b0;
    end
    npts0=npts1(m)+npts2(s)-1;
else
    b0=-delta*(npts1(s)-1);
    if(reltime); b0=b01(s)-b01(m)+b0;
    else b0=timediff(cell2mat(b1(m)),cell2mat(b1(s)),'utc')+b0;
    end
    npts0=npts1(m)+npts1(s)-1;
end

% .name
if(~peaks)
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

% loop over pairs
[depmin,depmen,depmax]=deal(nan(npairs,1));
if(verbose)
    disp('Correlating Records');
    print_time_left(0,npairs);
end
for i=1:npairs
    % skip dataless
    if((twodata && (npts1(m(i))==0 || npts2(s(i))==0)) ...
            || (~twodata && (npts1(m(i))==0 || npts1(s(i))==0)))
        if(peaks)
            [p(i).cg,p(i).lg]=getpeaks(zeros(0,1),peaksopt{:});
            p(i).pg=nan(size(lg));
        else
            data(i).dep=oclass(zeros(0,1));
        end
        
        % shortcut
        if(verbose); print_time_left(i,npairs); end
        continue;
    end
    
    % convert to fd cmplx, get zlac if normalizing
    if(~cmplx1(m(i)))
        oclass1{m(i)}=class(master(m(i)).dep);
        if(normxc)
            zlac1(m(i))=sqrt(sum(double(master(m(i)).dep).^2));
            if(zlac1(m(i))==0); zlac1(m(i))=eps; end
        end
        if(twodata)
            master(m(i)).dep=conj(fft(double(master(m(i)).dep),nspts,1));
        else
            master(m(i)).dep=fft(double(master(m(i)).dep),nspts,1);
        end
        cmplx1(m(i))=true;
    end
    if(twodata)
        if(~cmplx2(s(i)))
            if(normxc)
                zlac2(s(i))=sqrt(sum(double(slave(s(i)).dep).^2));
                if(zlac2(s(i))==0); zlac2(s(i))=eps; end
            end
            slave(s(i)).dep=fft(double(slave(s(i)).dep),nspts,1);
            cmplx2(s(i))=true;
        end
    else
        if(~cmplx1(s(i)))
            oclass1{s(i)}=class(master(s(i)).dep);
            if(normxc)
                zlac1(s(i))=sqrt(sum(double(master(s(i)).dep).^2));
                if(zlac1(s(i))==0); zlac1(s(i))=eps; end
            end
            master(s(i)).dep=fft(double(master(s(i)).dep),nspts,1);
            cmplx1(s(i))=true;
        end
    end
    
    % returning peaks or correlations
    if(peaks)
        % correlate
        if(twodata)
            tmp=ifft(master(m(i)).dep.*slave(s(i)).dep,[],1,'symmetric');
            if(normxc); tmp=tmp./(zlac1(m(i))*zlac2(s(i))); end
            tmp=tmp([nspts-npts2(s(i))+2:nspts 1:npts1(m(i))],1);
        else
            tmp=ifft(conj(master(m(i)).dep).*master(s(i)).dep,[],1,...
                'symmetric');
            if(normxc); tmp=tmp./(zlac1(m(i))*zlac1(s(i)));
            end
            tmp=tmp([nspts-npts1(s(i))+2:nspts 1:npts1(m(i))],1);
        end
        
        % convert back to original class
        oclass=str2func(oclass1{m(i)});
        tmp=oclass(tmp);
        
        % lagrng
        if(~isempty(lagrng))
            lidx=ceil(...
                (lagrng(1)-b0(i))/delta):floor((lagrng(2)-b0(i))/delta);
            tmp=[zeros(sum(lidx<=0),1); ...
                tmp(lidx(lidx>0 & lidx<=npts0(i))); ...
                zeros(sum(lidx>npts0(i)),1)];
            npts0(i)=numel(tmp);
            b0(i)=b0(i)+(min(lidx)-1)*delta;
        end
        
        % peaks
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
        p(i).lg=b0(i)+delta*(p(i).lg-1);
    else
        % correlate
        if(twodata)
            data(i).dep=ifft(master(m(i)).dep.*slave(s(i)).dep,...
                [],1,'symmetric');
            if(normxc)
                data(i).dep=data(i).dep./(zlac1(m(i))*zlac2(s(i)));
            end
            data(i).dep=data(i).dep(...
                [nspts-npts2(s(i))+2:nspts 1:npts1(m(i))],1);
        else
            data(i).dep=ifft(conj(master(m(i)).dep).*master(s(i)).dep,...
                [],1,'symmetric');
            if(normxc)
                data(i).dep=data(i).dep./(zlac1(m(i))*zlac1(s(i)));
            end
            data(i).dep=data(i).dep(...
                [nspts-npts1(s(i))+2:nspts 1:npts1(m(i))],1);
        end
        
        % convert back to original class
        oclass=str2func(oclass1{m(i)});
        data(i).dep=oclass(data(i).dep);
        
        % lagrng
        if(~isempty(lagrng))
            lidx=ceil(...
                (lagrng(1)-b0(i))/delta):floor((lagrng(2)-b0(i))/delta);
            data(i).dep=[zeros(sum(lidx<=0),1); ...
                data(i).dep(lidx(lidx>0 & lidx<=npts0(i))); ...
                zeros(sum(lidx>npts0(i)),1)];
            npts0(i)=numel(data(i).dep);
            b0(i)=b0(i)+(min(lidx)-1)*delta;
        end
        
        % absolute?
        if(absxc); data(i).dep=abs(data(i).dep); end
        
        % dep*
        depmen(i)=mean(data(i).dep(~isnan(data(i).dep)));
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
if(twodata)
    data=changeheader(data,'z',z1(m),...
        'delta',delta,'npts',npts0,'b',b0,'e',b0+delta*(npts0-1),...
        'kuser0','MASTER','user0',m,'kuser1','SLAVE','user1',s,...
        'depmin',depmin,'depmen',depmen,'depmax',depmax,...
        'a utc',b1(m),'f utc',e1(m),'t0 utc',b2(s),'t1 utc',e2(s),...
        'ev',st1(m,:),'st',st2(s,:),'kt0',kname1(m,1),'kt1',kname1(m,2),...
        'kt2',kname1(m,3),'kt3',kname1(m,4),'kname',kname2(s,:),...
        'user2',cmp1(m,1),'user3',cmp1(m,2),'cmp',cmp2(s,:));
else
    data=changeheader(data,'z',z1(m),...
        'delta',delta,'npts',npts0,'b',b0,'e',b0+delta*(npts0-1),...
        'kuser0','MASTER','user0',m,'kuser1','SLAVE','user1',s,...
        'depmin',depmin,'depmen',depmen,'depmax',depmax,...
        'a utc',b1(m),'f utc',e1(m),'t0 utc',b1(s),'t1 utc',e1(s),...
        'ev',st1(m,:),'st',st1(s,:),'kt0',kname1(m,1),'kt1',kname1(m,2),...
        'kt2',kname1(m,3),'kt3',kname1(m,4),'kname',kname1(s,:),...
        'user2',cmp1(m,1),'user3',cmp1(m,2),'cmp',cmp1(s,:));
end

% update delaz info
checkheader_state(true);
data=checkheader(data,'all','ignore','old_delaz','fix');
checkheader_state(false);

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

