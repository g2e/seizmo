function [badpz,sacpzdb]=bad_sacpz_cplxpair(sacpzdb)
%BAD_SACPZ_CPLXPAIR    Finds SAC Polezero files with bad complex pairs
%
%    Usage:    badpz=bad_sacpz_cplxpair(sacpzdb)
%              [badpz,sacpzdb]=bad_sacpz_cplxpair(sacpzdb)
%
%    Description:
%     BADPZ=BAD_SACPZ_CPLXPAIR(SACPZDB) will search the SAC PoleZero
%     database SACPZDB (made via IRIS_SACPZDB_BUILD) for responses that
%     have an unpaired complex pole/zero.  Basically, all poles & zeros
%     with an imaginary component are expected to have a corresponding
%     complex conjugate to pair with.  If that complex conjugate is missing
%     then the response is deemed bad.  Usually these are due to typos in
%     the response info (forgetting to flip the sign, switching 2 digits,
%     or the exponent is off by 1).  This will not find typos for
%     poles/zeros that have no imaginary component.
%
%     [BADPZ,SACPZDB]=BAD_SACPZ_CPLXPAIR(SACPZDB) returns SACPZDB with any
%     PoleZero responses deemed bad removed.
%
%    Notes:
%
%    Examples:
%     % Check the db included in SEIZMO:
%     sacpzdb=load('sacpzdb');
%     net=fieldnames(sacpzdb);
%     cnt=0; for i=1:numel(net); cnt=cnt+numel(sacpzdb.(net{i})); end
%     bad=numel(bad_sacpz_cplxpair(sacpzdb));
%     disp([num2str(bad) ' of ' num2str(cnt) ' POLEZEROS ARE BAD (' ...
%         num2str(bad/cnt*100,'%5.2g') '%)']);
%     % See IRIS_SACPZDB_FIXES for a list of some found in the past.
%
%     % Checking your own db is simple:
%     mysacpzdb=load('mysacpzdb');
%     mydb.CUSTOM=mysacpzdb; % if you don't have a network layer...
%     badpz=bad_sacpz_cplxpair(mydb);
%
%    See also: IRIS_SACPZDB_FIXES, IRIS_SACPZDB_BUILD, GETSACPZ,
%              MAKESACPZDB

%     Version History:
%        May  28, 2010 - initial version
%        Feb.  3, 2012 - doc update, return sacpzdb w/o bad
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  3, 2012 at 01:25 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% find bad sac polezeros
nets=fieldnames(sacpzdb);
cnt=0;
for m=1:numel(nets)
    netcnt=0;
    fprintf(nets{m})
    npz=numel(sacpzdb.(nets{m}));
    keep=true(npz,1);
    for n=1:npz
        try
            cplxpair(sacpzdb.(nets{m})(n).p);
            cplxpair(sacpzdb.(nets{m})(n).z);
        catch
            cnt=cnt+1;
            netcnt=netcnt+1;
            badpz(cnt).net=nets{m};
            badpz(cnt).idx=n;
            badpz(cnt).sacpz=sacpzdb.(nets{m})(n);
            keep(n)=false;
        end
    end
    sacpzdb.(nets{m})=sacpzdb.(nets{m})(keep);
    fprintf([' --> Found ' num2str(netcnt) ' Bad SAC PoleZero(s)\n']);
end

end
