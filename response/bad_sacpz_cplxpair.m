function [badpz]=bad_sacpz_cplxpair(sacpzdb)
%BAD_SACPZ_CPLXPAIR    Finds SAC Polezero files with bad complex pairs
%
%    Usage:    badpz=bad_sacpz_cplxpair(sacpzdb)
%
%    Description: BADPZ=BAD_SACPZ_CPLXPAIR(SACPZDB) will search the SAC
%     PoleZero database SACPZDB (made via IRIS_SACPZDB_HOWTO) for responses
%     that have an unpaired complex pole or zero.  Basically, all poles &
%     zeros with an imaginary component are expected to have a
%     corresponding complex conjugate to pair with.  If that complex
%     conjugate is missing then the response is deemed bad.  Usually these
%     are due to typos in the response info (forgetting to flip the sign,
%     switching 2 digits, or the exponent is off by 1).  This will not find
%     typos for poles/zeros that have no imaginary component.
%
%    Notes:
%
%    Examples:
%     See IRIS_SACPZDB_FIXES for a list of some that I have found in the
%     past.  Checking your own db is simple:
%      mysacpzdb=load('mysacpzdb');
%      mydb.CUSTOM=mysacpzdb; % if you don't have a network layer...
%      badpz=bad_sacpz_cplxpair(mydb);
%
%    See also: IRIS_SACPZDB_HOWTO, IRIS_SACPZDB_FIXES, GETSACPZ

%     Version History:
%        May  28, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  28, 2010 at 01:25 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% find bad sac polezeros
nets=fieldnames(sacpzdb);
cnt=0;
for m=1:numel(nets)
    netcnt=0;
    fprintf(nets{m})
    for n=1:numel(sacpzdb.(nets{m}))
        try
            cplxpair(sacpzdb.(nets{m})(n).p);
            cplxpair(sacpzdb.(nets{m})(n).z);
        catch
            cnt=cnt+1;
            netcnt=netcnt+1;
            badpz(cnt).net=nets{m};
            badpz(cnt).idx=n;
            badpz(cnt).sacpz=sacpzdb.(nets{m})(n);
        end
    end
    fprintf([' --> Found ' num2str(netcnt) ' Bad SAC PoleZero(s)\n']);
end

end
