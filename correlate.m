function [data]=correlate(data1,varargin)
%CORRELATE    Compute cross correlograms of SEIZMO data records
%
%    Usage:    correlograms=correlate(data)
%              correlograms=correlate(data1,data2)
%              correlograms=correlate(...,'normxc',true|false,...)
%              correlograms=correlate(...,'lags',range,...)
%              peaks=correlate(...,'npeaks',n,...)
%              peaks=correlate(...,'spacing',seconds,...)
%              peaks=correlate(...,'lags',range,...)
%              peaks=correlate(...,'adjacent',halfwidth,...)
%              peaks=correlate(...,'absxc',true|false,...)
%              peaks=correlate(...,'normxc',true|false,...)
%
%    Description: CORRELOGRAMS=CORRELATE(DATA) cross correlates every
%     possible pairing of records in DATA.  If the number of records in
%     DATA is N then CORRELOGRAMS is a SEIZMO dataset of N*(N-1)/2 correlo-
%     grams.  Correlograms are normalized using zero-lag autocorrelation
%     values.  The correlograms are named using the following format:
%      CORR_-_MASTER_-_REC<idx>_-_<stninfo>_-_SLAVE_-_REC<idx>_-_<stninfo>
%     where <idx> is the index of the record in DATA (zero padded) and
%     <stninfo> is the fields knetwk, kstnm, khole, kcmpnm of record <idx>
%     joined with periods ('.') in between.  The path is set to the current
%     directory ('.'), while byte-order uses that which is native to the
%     current system.  Default filetype is SAC v6 binary file.  See the
%     Header changes section for details on info inserted into the header.
%     See Notes for details on requirements for the input dataset.
%     
%     CORRELOGRAMS=CORRELATE(DATA1,DATA2) cross correlates every record in
%     DATA1 against every record in DATA2 (meaning that records in DATA1
%     are the 'master' records).  Thus if there are N1 records in DATA1 and
%     N2 records in DATA2, then CORRELOGRAMS contains N1*N2 correlograms.
%     See the usage description above and in Header changes section for
%     further output details.  See Notes for details on requirements for
%     the input dataset(s).
%
%     *******************************************************************
%     THE FOLLOWING 2 OPTION DESCRIPTIONS ASSUME YOU HAVE NPEAKS UNSET!
%     *******************************************************************
%
%     CORRELOGRAMS=CORRELATE(...,'NORMXC',TRUE|FALSE,...) allows selecting
%     if the output correlograms are normalized or not.  The default is
%     TRUE, which returns correlograms that are normalized by their
%     autocorrelations.  Setting the NORMXC option to FALSE outputs
%     unnormalized correlograms.
%
%     CORRELOGRAMS=CORRELATE(...,'LAGS',RANGE,...) allows specifying the
%     lag range in seconds.  By default RANGE is set to the maximum
%     possible range (which is from -1*length_of_longest_master_record to
%     length_of_longest_slave_record).  This can be specified by setting
%     RANGE to [].  Giving a single-valued RANGE specifies a symmetric
%     lag window while a two-element RANGE gives an asymmetric lag window.
%
%     *******************************************************************
%     THE FOLLOWING OPTION DESCRIPTIONS ASSUME YOU HAVE SET NPEAKS>0!
%     *******************************************************************
%
%     PEAKS=CORRELATE(...,'NPEAKS',NPEAKS,...) allows turning on the
%     correlogram peak picker in the underlying function MCXC, returning
%     peak info rather than correlograms.  If NPEAKS==0 (the default), no
%     picking is done and correlograms are returned.  Setting NPEAKS>0 will
%     return info for the specified number of peaks from each correlogram.
%     PEAKS is a structure containing 3 fields: 'cg' contains the
%     correlation value of the peaks, 'lg' contains the lag time (in
%     seconds) of peaks, and 'pg' contains the polarity of peaks (useful
%     when ABSXC is set to TRUE -- the default).  Output dimensions depend
%     on the number of input datasets.  If a single dataset is given to
%     correlate, then the matrices are column vectors with strongest peak
%     info for all correlograms in the first sheet of the third dimension.
%     The second strongest peak info, if NPEAKS>=2, is returned in the
%     second sheet of the third dimension and so on for all lesser peaks so
%     that PEAKS.cg(:,:,N), PEAKS.lg(:,:,N), and PEAKS.pg(:,:,N) give info
%     about the Nth highest peak of every correlogram.  Get info on a
%     specific correlogram peak utilizing the following indexing example:
%      PEAKS.cg(sum(N+1-MASTER_IDX:N-1)+SLAVE_IDX,1,NTHPEAK)
%     where MASTER_IDX and SLAVE_IDX give the indices of the records used
%     in generating the correlogram, N is the number of records in DATA,
%     and NTHPEAK is the index of the peak of interest.  Correlating two
%     datasets against each other returns a similar structure with 
%     PEAKS.cg(:,:,N), PEAKS.lg(:,:,N), and PEAKS.pg(:,:,N) giving info
%     about the Nth highest peak of every correlogram.  The difference is
%     that the peak info for a specific correlogram is easy to extract in
%     this case: PEAKS.cg(SLAVE_IDX,MASTER_IDX,NTHPEAK).
%
%     PEAKS=CORRELATE(...,'SPACING',SECONDS,...) controls the minimum
%     spacing between returned peaks.  By default SPACING is set to 1
%     sample interval which requires that peaks be at least 1 sample
%     interval apart.  Setting SPACING==0 will cause the peak picker to
%     always return the same peak for all the peaks.  Setting SPACING too
%     high can cause the peak picker to return a peak with value=0, lag=0,
%     polarity=0 when no points are left satisfying the SPACING request.
%
%     PEAKS=CORRELATE(...,'LAGS',RANGE,...) controls the lag range that the
%     peak picker is allowed to search.  By default LAGS is set to the
%     maximum possible range.  LAGS can be a one (symmetric window) or two
%     (asymmetric window) element vector indicating the range of lags (in 
%     seconds!) to search.  Specify the default with [].
%
%     PEAKS=CORRELATE(...,'ADJACENT',HALFWIDTH,...) controls the number of
%     adjacent points to be returned in addition to each pick picked.  By
%     default HALFWIDTH is set to 0 which requires that all points within 0
%     seconds of each peak be returned too.  Setting ADJACENT to 2 seconds
%     will return the peak and all points at and within 2 seconds.  The
%     peaks and their adjacent points info are stored down the 4th
%     dimension of the fields in PEAKS and are kept in their original order
%     such that the peak is always in the middle.  For example, 
%     CG(:,:,1,ceil(end/2)) would give the highest peak cross correlation
%     value for each correlogram (note that ceil(end/2) will always give
%     the 4th dimension index of the peak values).  Adjacent points may be
%     set to value=0, lag=?, polarity=0 when the point comes within SPACING
%     of another peak or if the point falls outside the allowed lag range.
%
%     PEAKS=CORRELATE(...,'ABSXC',TRUE|FALSE,...) controls whether the peak
%     picker looks at absolute peaks or just positive peaks.  By default
%     the ABSXC option is set to TRUE which causes the peak picker to first
%     take the absolute value of the correlograms before picking the
%     highest peaks.  Polarities are returned in the PG matrix.  Setting
%     ABSXC to FALSE will cause the peak picker to search for the highest
%     peak in the unaltered correlograms and so PG will always just contain
%     ones.
%
%     PEAKS=CORRELATE(...,'NORMXC',TRUE|FALSE,...) controls whether the
%     cross correlations are to be normalized (using the autocorrelations)
%     or not before being picked.  By default the NORMXC option is set to
%     TRUE.  Correlation values for a normalized correlogram are in the
%     range from -1 to 1, with 1 being a perfect correlation and -1 being a
%     perfect anticorrelation.
%
%    Notes:
%     - All records are required to have a common sample rate (DELTA field
%       should be equal), be evenly sampled (LEVEN field should be TRUE),
%       be time series or xy files (IFTYPE field should be itime or ixy),
%       and single component (GETNCMP should return all ones).  All records
%       are also passed through CHECKHEADER so sane settings for all header
%       fields are enforced.
%     - CORRELATE ignores the reference timing of records when computing
%       lag times (or rather it assumes records are synced to the same
%       reference time).  Synchronize the dataset to a single reference
%       time to get lags based on absolute timing.
%
%    Header Changes: ALL (see BSEIZMO for defaults not listed below)
%     DEPMEN, DEPMIN, DEPMAX are updated for the correlogram.
%     B, E give the lag range.
%     DELTA is the sample spacing.
%     Index of master record: USER0
%     Index of slave record:  USER1
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
%       NZ                   NZ
%      MASTER RECORD FIELD  CORRELOGRAM FIELD
%       STLA                 EVLA
%       STLO                 EVLO
%       STEL                 EVEL
%       STDP                 EVDP
%       KNETWK               KT0
%       KSTNM                KT1
%       KHOLE                KT2
%       KCMPNM               KT3
%
%    Examples:
%     Basically equivalent to 'correlate' in SAC:
%      correlograms=correlate(data(1),data,'normxc',false)
%
%     Align all correlograms on their strongest peak:
%      p0(timeshift(correlate(data),...
%       -getfield(correlate(data,'npeaks',1),'lg')))
%
%    See also: CONVOLVE, DFT, IDFT, MCXC

%     Version History:
%        June 27, 2009 - first fully functional version
%        Oct.  7, 2009 - slave records NZ* info passed on, fixed record
%                        names in 2 dataset correlogram case, fixed ST/EV
%                        info in 1 dataset correlogram case
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct.  7, 2009 at 18:15 GMTs

% todo:

% check nargin
if(nargin<1)
    error('seizmo:selectrecords:notEnoughInputs',...
        'Not enough input arguments.');
end

% check data structure
msg=seizmocheck(data1,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data1=checkheader(data1);

% 1 or 2 datasets
nrecs1=numel(data1); onedata=true;
if(nargin>1 && isseizmo(varargin{1},'dep'))
    data2=varargin{1};
    varargin(1)=[];
    nrecs2=numel(data2); 
    onedata=false;
    data2=checkheader(data2);
    delta=cross_check_data(data1,data2);
else
    delta=cross_check_data(data1);
end

% turn off header checking
oldcheckheaderstate=get_checkheader_state;
set_checkheader_state(false);

% check varargin
nvarargin=numel(varargin);
if(mod(nvarargin,2))
    error('seizmo:correlate:OptionMustBePaired',...
        'Options must be paired with a value!');
end

% how many peaks
npeaks=0; % default is no peak picking
for i=1:2:nvarargin
    if(ischar(varargin{i}))
        switch lower(varargin{i})
            case 'npeaks'
                npeaks=varargin{i+1};
            case 'lags'
                varargin{i+1}=varargin{i+1}/delta;
            case 'spacing'
                varargin{i+1}=varargin{i+1}/delta;
            case 'adjacent'
                varargin{i+1}=varargin{i+1}/delta;
        end
    end
end

% check npeaks
if(~isscalar(npeaks) && ~isnumeric(npeaks) && fix(npeaks)~=npeaks)
    error('seizmo:correlate:badInput',...
        'Option NPEAKS must be a scalar integer')
end

% split based on npeaks
% npeaks>0 ==> assign grids to new fields
if(npeaks)
    % split based on ndatasets
    if(onedata)
        % get relative start times (disregarding absolute timing)
        b=getheader(data1,'b');
        bdiff=b(:,ones(nrecs1,1))-b(:,ones(nrecs1,1)).';
        bdiff=bdiff(tril(true(nrecs1),-1));
        % extract records
        data1=records2mat(data1);
        % get correlation peaks
        [data.cg,data.lg,data.pg]=mcxc(data1,varargin{:});
        % adjust lags
        s=size(data.lg); s(numel(s)+1:4)=1;
        data.lg=data.lg*delta+bdiff(:,:,ones(s(3),1),ones(s(4),1));
    else % two datasets
        % get relative start times (disregarding absolute timing)
        b1=getheader(data1,'b');
        b2=getheader(data2,'b');
        bdiff=b2(:,ones(nrecs1,1))-b1(:,ones(nrecs2,1)).';
        % extract records
        data1=records2mat(data1);
        data2=records2mat(data2);
        % get correlation peaks
        [data.cg,data.lg,data.pg]=mcxc(data1,data2,varargin{:});
        % make life easier (essentially an output bug)
        s=size(data.lg); ndims=numel(s);
        data.cg=permute(data.cg,[2 1 3:ndims]);
        data.lg=permute(data.lg,[2 1 3:ndims]);
        data.pg=permute(data.pg,[2 1 3:ndims]);
        % adjust lags
        s(ndims+1:4)=1;
        data.lg=data.lg*delta+bdiff(:,:,ones(s(3),1),ones(s(4),1));
    end
% npeaks==0 ==> replace dataset with correlograms
else
    % split based on ndatasets
    if(onedata)
        % extract header info needed
        [knetwk,kstnm,khole,kcmpnm,b,stla,stlo,stdp,stel,...
            nzyear,nzjday,nzhour,nzmin,nzsec,nzmsec]=...
            getheader(data1,'knetwk','kstnm','khole','kcmpnm','b',...
            'stla','stlo','stdp','stel','nzyear','nzjday','nzhour',...
            'nzmin','nzsec','nzmsec');
        
        % extract records
        data1=records2mat(data1);
        
        % get correlation peaks
        [cg,lg]=mcxc(data1,varargin{:});
        
        % separate records
        [ncors,one,nlags]=size(cg);
        if(ncors==1)
            cg=squeeze(cg);
        else
            cg=squeeze(cg).';
        end
        cg=mat2cell(cg,nlags,ones(ncors,1));
        
        % make lags for each record
        % 1 to nrecs-1 outer loop is master
        % 2 to nrecs inner loop is slave
        bdiff=b(:,ones(nrecs1,1)).'-b(:,ones(nrecs1,1));
        bdiff=bdiff(tril(true(nrecs1),-1));
        lg=lg(:,ones(ncors,1))*delta-bdiff(:,ones(nlags,1)).';
        lg=mat2cell(lg,nlags,ones(ncors,1));
        
        % make a new dataset
        data=[lg; cg];
        data=bseizmo(data{:});
        
        % get names
        i=1:nrecs1;
        idx=cellstr(num2str(i.',...
            ['%0' num2str(ceil(log10(nrecs1+1))) 'd']));
        idx=idx(:,ones(nrecs1,1)).';
        knetwk=knetwk(:,ones(nrecs1,1));
        kstnm=kstnm(:,ones(nrecs1,1));
        khole=khole(:,ones(nrecs1,1));
        kcmpnm=kcmpnm(:,ones(nrecs1,1));
        names=strcat({'CORR_-_MASTER_-_REC'},idx,{'_-_'},...
            knetwk',{'.'},kstnm',{'.'},khole',{'.'},kcmpnm',...
            {'_-_SLAVE_-_REC'},idx',{'_-_'},...
            knetwk,{'.'},kstnm,{'.'},khole,{'.'},kcmpnm);
        names=names(tril(true(nrecs1),-1));
        [data.name]=deal(names{:});
        
        % add record numbers to header
        midx=i(ones(nrecs1,1),:);
        sidx=midx.';
        midx=midx(tril(true(nrecs1),-1));
        sidx=sidx(tril(true(nrecs1),-1));
        
        % add reference times to header
        nzyear=nzyear(:,ones(nrecs1,1));
        nzjday=nzjday(:,ones(nrecs1,1));
        nzhour=nzhour(:,ones(nrecs1,1));
        nzmin=nzmin(:,ones(nrecs1,1));
        nzsec=nzsec(:,ones(nrecs1,1));
        nzmsec=nzmsec(:,ones(nrecs1,1));
        nzyear=nzyear(tril(true(nrecs1),-1));
        nzjday=nzjday(tril(true(nrecs1),-1));
        nzhour=nzhour(tril(true(nrecs1),-1));
        nzmin=nzmin(tril(true(nrecs1),-1));
        nzsec=nzsec(tril(true(nrecs1),-1));
        nzmsec=nzmsec(tril(true(nrecs1),-1));
        
        % add station id to header
        mknetwk=knetwk'; mknetwk=mknetwk(tril(true(nrecs1),-1));
        mkstnm=kstnm'; mkstnm=mkstnm(tril(true(nrecs1),-1));
        mkhole=khole'; mkhole=mkhole(tril(true(nrecs1),-1));
        mkcmpnm=kcmpnm'; mkcmpnm=mkcmpnm(tril(true(nrecs1),-1));
        knetwk=knetwk(tril(true(nrecs1),-1));
        kstnm=kstnm(tril(true(nrecs1),-1));
        khole=khole(tril(true(nrecs1),-1));
        kcmpnm=kcmpnm(tril(true(nrecs1),-1));
        
        % add station locations to header
        stla=stla(:,ones(nrecs1,1)); evla=stla.';
        stlo=stlo(:,ones(nrecs1,1)); evlo=stlo.';
        stel=stel(:,ones(nrecs1,1)); evel=stel.';
        stdp=stdp(:,ones(nrecs1,1)); evdp=stdp.';
        stla=stla(tril(true(nrecs1),-1)); evla=evla(tril(true(nrecs1),-1));
        stlo=stlo(tril(true(nrecs1),-1)); evlo=evlo(tril(true(nrecs1),-1));
        stel=stel(tril(true(nrecs1),-1)); evel=evel(tril(true(nrecs1),-1));
        stdp=stdp(tril(true(nrecs1),-1)); evdp=evdp(tril(true(nrecs1),-1));
        
        % dep*
        depmen=nan(ncors,1); depmin=depmen; depmax=depmen;
        for i=1:ncors
            depmen(i)=mean(data(i).dep(:));
            depmin(i)=min(data(i).dep(:));
            depmax(i)=max(data(i).dep(:));
        end
        
        % update header
        data=changeheader(data,...
            'nzyear',nzyear,'nzjday',nzjday,'nzhour',nzhour,...
            'nzmin',nzmin,'nzsec',nzsec,'nzmsec',nzmsec,...
            'user0',midx,'kuser0','MASTER',...
            'user1',sidx,'kuser1','SLAVE',...
            'depmen',depmen,'depmin',depmin,'depmax',depmax,...
            'knetwk',knetwk,'kstnm',kstnm,'khole',khole,'kcmpnm',kcmpnm,...
            'kt0',mknetwk,'kt1',mkstnm,'kt2',mkhole,'kt3',mkcmpnm,...
            'stla',stla,'stlo',stlo,'stel',stel,'stdp',stdp,...
            'evla',evla,'evlo',evlo,'evel',evel,'evdp',evdp);
    else
        % extract header info needed
        [mknetwk,mkstnm,mkhole,mkcmpnm,mb,mstla,mstlo,mstdp,mstel]=...
            getheader(data1,'knetwk','kstnm','khole','kcmpnm','b',...
            'stla','stlo','stdp','stel');
        [sknetwk,skstnm,skhole,skcmpnm,sb,sstla,sstlo,sstdp,sstel,...
            nzyear,nzjday,nzhour,nzmin,nzsec,nzmsec]=...
            getheader(data2,'knetwk','kstnm','khole','kcmpnm','b',...
            'stla','stlo','stdp','stel','nzyear','nzjday','nzhour',...
            'nzmin','nzsec','nzmsec');
        
        % extract records
        data1=records2mat(data1);
        data2=records2mat(data2);
        
        % get correlation peaks
        [cg,lg]=mcxc(data1,data2,varargin{:});
        
        % separate records
        nlags=size(cg,3); 
        ncors=nrecs1*nrecs2;
        cg=permute(cg,[3 2 1]);
        cg=mat2cell(cg(:,:),nlags,ones(ncors,1));
        
        % make lags for each record
        bdiff=mb(:,ones(nrecs2,1)).'-sb(:,ones(nrecs1,1));
        bdiff=bdiff(:);
        lg=lg(:,ones(ncors,1))*delta-bdiff(:,ones(nlags,1)).';
        lg=mat2cell(lg,nlags,ones(ncors,1));
        
        % make a new dataset
        data=[lg; cg];
        data=bseizmo(data{:});
        
        % get names
        i=1:nrecs1;
        midx=cellstr(num2str(i.',...
            ['%0' num2str(ceil(log10(nrecs1+1))) 'd']));
        midx=midx(:,ones(nrecs2,1)).';
        j=1:nrecs2;
        sidx=cellstr(num2str(j.',...
            ['%0' num2str(ceil(log10(nrecs2+1))) 'd']));
        sidx=sidx(:,ones(nrecs1,1));
        mknetwk=mknetwk(:,ones(nrecs2,1)).';
        mkstnm=mkstnm(:,ones(nrecs2,1)).';
        mkhole=mkhole(:,ones(nrecs2,1)).';
        mkcmpnm=mkcmpnm(:,ones(nrecs2,1)).';
        sknetwk=sknetwk(:,ones(nrecs1,1));
        skstnm=skstnm(:,ones(nrecs1,1));
        skhole=skhole(:,ones(nrecs1,1));
        skcmpnm=skcmpnm(:,ones(nrecs1,1));
        names=strcat({'CORR_-_MASTER_-_REC'},midx,{'_-_'},...
            mknetwk,{'.'},mkstnm,{'.'},mkhole,{'.'},mkcmpnm,...
            {'_-_SLAVE_-_REC'},sidx,{'_-_'},...
            sknetwk,{'.'},skstnm,{'.'},skhole,{'.'},skcmpnm);
        [data.name]=deal(names{:});
        
        % add record numbers to header
        midx=i(ones(nrecs2,1),:); midx=midx(:);
        sidx=j(ones(nrecs1,1),:).'; sidx=sidx(:);
        
        % add station id to header
        sknetwk=sknetwk(:); skstnm=skstnm(:);
        skhole=skhole(:); skcmpnm=skcmpnm(:);
        mknetwk=mknetwk(:); mkstnm=mkstnm(:);
        mkhole=mkhole(:); mkcmpnm=mkcmpnm(:);
        
        % add station locations to header
        sstla=sstla(:,ones(nrecs1,1));
        sstlo=sstlo(:,ones(nrecs1,1));
        sstel=sstel(:,ones(nrecs1,1));
        sstdp=sstdp(:,ones(nrecs1,1));
        mstla=mstla(:,ones(nrecs2,1)).';
        mstlo=mstlo(:,ones(nrecs2,1)).';
        mstel=mstel(:,ones(nrecs2,1)).';
        mstdp=mstdp(:,ones(nrecs2,1)).';
        sstla=sstla(:); sstlo=sstlo(:);
        sstel=sstel(:); sstdp=sstdp(:);
        mstla=mstla(:); mstlo=mstlo(:);
        mstel=mstel(:); mstdp=mstdp(:);
        
        % add reference times to header
        nzyear=nzyear(:,ones(nrecs1,1)); nzyear=nzyear(:);
        nzjday=nzjday(:,ones(nrecs1,1)); nzjday=nzjday(:);
        nzhour=nzhour(:,ones(nrecs1,1)); nzhour=nzhour(:);
        nzmin=nzmin(:,ones(nrecs1,1)); nzmin=nzmin(:);
        nzsec=nzsec(:,ones(nrecs1,1)); nzsec=nzsec(:);
        nzmsec=nzmsec(:,ones(nrecs1,1)); nzmsec=nzmsec(:);
        
        % dep*
        depmen=nan(ncors,1); depmin=depmen; depmax=depmen;
        for i=1:ncors
            depmen(i)=mean(data(i).dep(:));
            depmin(i)=min(data(i).dep(:));
            depmax(i)=max(data(i).dep(:));
        end
        
        % update header
        data=changeheader(data,...
            'nzyear',nzyear,'nzjday',nzjday,'nzhour',nzhour,...
            'nzmin',nzmin,'nzsec',nzsec,'nzmsec',nzmsec,...
            'user0',midx,'kuser0','MASTER',...
            'user1',sidx,'kuser1','SLAVE',...
            'knetwk',sknetwk,'kstnm',skstnm,...
            'khole',skhole,'kcmpnm',skcmpnm,...
            'depmen',depmen,'depmin',depmin,'depmax',depmax,...
            'kt0',mknetwk,'kt1',mkstnm,'kt2',mkhole,'kt3',mkcmpnm,...
            'stla',sstla,'stlo',sstlo,'stel',sstel,'stdp',sstdp,...
            'evla',mstla,'evlo',mstlo,'evel',mstel,'evdp',mstdp);
    end
end

% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);
set_checkheader_state(oldcheckheaderstate);

end

function [delta]=cross_check_data(data1,data2)
% assure datasets are capable of being sensibly cross correlated

% if 2 datasets, datasets must be consistent with one another
if(nargin==2)
    delta=getheader([data1(:); data2(:)],'delta');
    iftype=getenumdesc([data1(:); data2(:)],'iftype');
    leven=getlgc([data1(:); data2(:)],'leven');
    ncmp=getncmp([data1(:); data2(:)]);
% if 1 dataset, needs to be consistent with itself
else
    [delta]=getheader(data1,'delta');
    iftype=getenumdesc(data1,'iftype');
    ncmp=getncmp(data1);
    leven=getlgc(data1,'leven');
end

% check delta,iftype,leven,ncmp
if(~isscalar(unique(delta)))
    error('seizmo:correlate:mixedSampleRates',...
        'Mixed sample rates not allowed!');
elseif(any(~strcmpi(iftype,'Time Series File')...
        & ~strcmpi(iftype,'General X vs Y file')))
    error('seizmo:correlate:illegalOperation',...
        'Illegal operation on spectral/xyz record(s)!');
elseif(any(~strcmpi(leven,'true')))
    error('seizmo:correlate:evenlySpacedOnly',...
        'Illegal operation on unevenly spaced data!');
elseif(any(ncmp>1))
    error('seizmo:correlate:tooManyComponents',...
        'Multi-component data is not supported!')
end

% return common sample rate
delta=delta(1);

end
