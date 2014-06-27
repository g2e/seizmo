function [badpz,sacpzdb]=bad_sacpz(sacpzdb,issue)
%BAD_SACPZ    Finds SAC Polezero files with bad responses
%
%    Usage:    badpz=bad_sacpz(sacpzdb,issue)
%              [badpz,sacpzdb]=bad_sacpz(sacpzdb,issue)
%
%    Description:
%     BADPZ=BAD_SACPZ(SACPZDB,ISSUE) will search the SAC PoleZero database
%     SACPZDB (made via IRIS_SACPZDB_BUILD) for responses that have an
%     issue as defined by ISSUE.  Allowed values for ISSUE:
%      'CONSTANT' - find responses with an Inf, NaN, or Zero constant
%      'CPLXPAIR' - find responses with 1+ unpaired complex poles/zeros
%      'FEWZEROS' - find seismometer responses with <3 zeros at the origin
%      'FLAT'     - find perfect seismometer responses (no pz, constant=1)
%      'ALL'      - find responses with any of the above problems
%      'BAD'      - CONSTANT+CPLXPAIR+FLAT (default)
%     BADPZ is a scalar struct with the fields:
%      .net   - Nx1 cellstr array of network associated to the bad polezero
%      .idx   - Nx1 double array of the indices to the bad polezeros
%      .sacpz - scalar struct containing the bad polezeros
%
%     [BADPZ,SACPZDB]=BAD_SACPZ_CPLXPAIR(SACPZDB,ISSUE) returns the SACPZDB
%     with an additional field .bad that is true for PoleZero responses
%     deemed bad.
%
%    Notes:
%     - Currently a staggering 24% of the SAC PoleZero files from IRIS have
%       issues ('BAD' option)!
%     - All poles & zeros with an imaginary component are expected to have
%       a corresponding complex conjugate to pair with.  If that complex
%       conjugate is missing then the response is deemed bad.  Usually
%       these are due to typos in the response info (e.g., forgetting to
%       flip the sign of the imaginary portion, switching 2 digits, or the
%       exponent is off by 1).  Note that we cannot find typos for
%       poles/zeros that have no imaginary component.
%
%    Examples:
%     % How to check & clean your own set of annotated polezero files:
%     [badpz,pz]=bad_sacpz(readsacpz_rdseed('my/pz/dir'));
%     goodpz=ssidx(pz,~pz.bad);
%
%    See also: IRIS_SACPZDB_FIXES, IRIS_SACPZDB_BUILD, FIX_BAD_SACPZ

%     Version History:
%        May  28, 2010 - initial version
%        Feb.  3, 2012 - doc update, return sacpzdb w/o bad
%        Mar.  6, 2014 - update for new sacpz struct format
%        Mar.  7, 2014 - change name to bad_sacpz, 2nd input for issue
%                        type, returned db has .bad rather than reduced
%        May  28, 2014 - added FLAT type, added FLAT to BAD, check that
%                        input is network database, update example, allow
%                        unnetworked polezero db
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  28, 2014 at 01:25 GMT

% todo:

% check nargin
error(nargchk(1,2,nargin));

% require sacpzdb is scalar struct
if(~isstruct(sacpzdb) || ~isscalar(sacpzdb))
    error('seizmo:bad_sacpz:badInput',...
        'SACPZDB must be a scalar struct!');
end

% convert polezero db to network polezero db
if(issacpz_rdseed(sacpzdb))
    netdb=false;
    tmp=sacpzdb;
    clear sacpzdb;
    sacpzdb.CUSTOM=tmp;
    clear tmp;
else
    netdb=true;
end

% default/check issue
valid.ISSUE={'CPLXPAIR' 'CONSTANT' 'FEWZEROS' 'FLAT' 'ALL' 'BAD'};
if(nargin==1 || isempty(issue)); issue='bad'; end
if(~ismember(upper(issue),valid.ISSUE))
    error('seizmo:bad_sacpz:badInput',...
        'Unknown ISSUE: %s',issue);
end

% make channel codes for translational seismometers
valid.SEISCODE=cell(6*2*13,1);
cnt=0;
for i='BHMSUV'
    for j='HL'
        for k='0123456789ENZ'
            cnt=cnt+1;
            valid.SEISCODE{cnt}=[i j k];
        end
    end
end

% find bad sac polezeros
nets=fieldnames(sacpzdb);
totnpz=0;
cnt=0;
for m=1:numel(nets)
    if(netdb); fprintf(nets{m}); end
    if(~issacpz_rdseed(sacpzdb.(nets{m})))
        error('seizmo:bad_sacpz:badInput',...
            'Input is not a valid SACPZ network database!');
    end
    npz=numel(sacpzdb.(nets{m}).k);
    totnpz=totnpz+npz;
    
    switch upper(issue)
        case 'CPLXPAIR'
            bad=false(npz,1);
            for n=1:npz
                try
                    cplxpair(sacpzdb.(nets{m}).p{n});
                    cplxpair(sacpzdb.(nets{m}).z{n});
                catch
                    bad(n)=true;
                end
            end
        case 'CONSTANT'
            bad=isnan(sacpzdb.(nets{m}).k) | isinf(sacpzdb.(nets{m}).k) ...
                | sacpzdb.(nets{m}).k==0;
        case 'FEWZEROS'
            bad=false(npz,1);
            for n=1:npz
                if(ismember(sacpzdb.(nets{m}).kcmpnm(n),valid.SEISCODE) ...
                        && sum(sacpzdb.(nets{m}).z{n}==0)<3)
                    bad(n)=true;
                end
            end
        case 'FLAT'
            bad=ismember(sacpzdb.(nets{m}).kcmpnm,valid.SEISCODE) & ...
                sacpzdb.(nets{m}).k==1 ...
                & cellfun('isempty',sacpzdb.(nets{m}).z) ...
                & cellfun('isempty',sacpzdb.(nets{m}).p);
        case 'ALL'
            bad=isnan(sacpzdb.(nets{m}).k) | isinf(sacpzdb.(nets{m}).k) ...
                | sacpzdb.(nets{m}).k==0 ...
                | (ismember(sacpzdb.(nets{m}).kcmpnm,valid.SEISCODE) ...
                & sacpzdb.(nets{m}).k==1 ...
                & cellfun('isempty',sacpzdb.(nets{m}).z) ...
                & cellfun('isempty',sacpzdb.(nets{m}).p));
            for n=1:npz
                if(ismember(sacpzdb.(nets{m}).kcmpnm(n),valid.SEISCODE) ...
                        && sum(sacpzdb.(nets{m}).z{n}==0)<3)
                    bad(n)=true;
                else
                    try
                        cplxpair(sacpzdb.(nets{m}).p{n});
                        cplxpair(sacpzdb.(nets{m}).z{n});
                    catch
                        bad(n)=true;
                    end
                end
            end
        case 'BAD'
            bad=isnan(sacpzdb.(nets{m}).k) | isinf(sacpzdb.(nets{m}).k) ...
                | sacpzdb.(nets{m}).k==0 ...
                | (ismember(sacpzdb.(nets{m}).kcmpnm,valid.SEISCODE) ...
                & sacpzdb.(nets{m}).k==1 ...
                & cellfun('isempty',sacpzdb.(nets{m}).z) ...
                & cellfun('isempty',sacpzdb.(nets{m}).p));
            for n=1:npz
                try
                    cplxpair(sacpzdb.(nets{m}).p{n});
                    cplxpair(sacpzdb.(nets{m}).z{n});
                catch
                    bad(n)=true;
                end
            end
    end
    netcnt=sum(bad);
    if(netcnt>0)
        badpz.net(cnt+1:cnt+netcnt,1)=nets(m);
        badpz.idx(cnt+1:cnt+netcnt,1)=find(bad);
        if(cnt==0)
            badpz.sacpz=ssidx(sacpzdb.(nets{m}),bad);
        else
            badpz.sacpz=sscat(badpz.sacpz,...
                ssidx(sacpzdb.(nets{m}),bad));
        end
    end
    cnt=cnt+netcnt;
    sacpzdb.(nets{m}).bad=bad;
    if(netdb)
        fprintf([' --> Found ' num2str(netcnt) ' out of ' ...
            num2str(npz) ' are Bad SAC PoleZero(s)\n']);
    end
end
if(netdb); fprintf('TOTAL'); end
fprintf([' --> Found ' num2str(cnt) ' out of ' ...
    num2str(totnpz) ' are Bad SAC PoleZero(s)\n']);

% flatten unnetworked polezero db
if(~netdb); sacpzdb=sacpzdb.CUSTOM; end

% no bad sacpz?
if(~exist('badpz','var'))
    badpz=[];
end

end
