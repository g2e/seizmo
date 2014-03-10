function [sacpzdb]=fix_bad_sacpz(sacpzdb)
%FIX_BAD_SACPZ    Fixes SOME cases of bad SAC PoleZeros
%
%    Usage:    sacpzdb=fix_bad_sacpz(sacpzdb)
%
%    Description:
%     SACPZDB=FIX_BAD_SACPZ(SACPZDB) implements simple fixes to the IRIS
%     SAC PoleZero database as mentioned in IRIS_SACPZDB_FIXES.  These may
%     or may not be correct though - see the Notes section below.
%
%    Notes:
%     - These should be reported to IRIS but I'm not really familiar with
%       how to do that, nor willing to put in the time required.  Once I
%       finish my thesis expect the thousands of problems I've found (many
%       of which are not corrected here) to be reported.
%     - IRIS fixed the bad polezeros of the GE network by providing
%       completely different responses whereas mine was just a sign
%       change to a complex value.
%
%    Examples:
%     % See if we reduce the bad count:
%     numel(bad_sacpz(sacpzdb))
%     numel(bad_sacpz(fix_bad_sacpz(sacpzdb)))
%
%    See also: BAD_SACPZ, IRIS_SACPZDB_FIXES, IRIS_SACPZDB_BUILD

%     Version History:
%        May  28, 2010 - initial version
%        Feb.  3, 2012 - doc update
%        Mar.  6, 2014 - update for new sacpz struct format
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  6, 2014 at 02:25 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% the fixes
% - looks like the GE net is fixed! The responses are completely different.
badnet={'GE' 'PN' 'CZ' 'AF' 'AZ'};
badfield={'z' 'p' 'p' 'p' 'p'};
badvalue={
[-15.15; -176.6; -463.1 + 430.5i; -463.1 - 430.5i; 0; 0; 0]
[-0.1486 + 0.1486i; -0.1486 - 0.1486i; -414.69; -999.027 - 999.027i]
[-0.0055; -0.0173 + 0.0179i; -0.0173 - 0.0179i; 
 -7.212 + 17.415i; -7.212 - 17.415i; -17.415 + 7.212i; -17.415 - 7.212i;
 -23.4459; -37.0379; -45.4483; -8.7307 + 42.4144i; -8.7307 - 41.4144i]
[-0.9211 + 0.94i; -0.9211 + 0.94i]
[-981 + 1009i; -981 - 1009i; -3290 - 1263i; -3290 - 1263i]};
goodvalue={
[-15.15; -176.6; -463.1 + 430.5i; -463.1 - 430.5i; 0; 0; 0]
[-0.1486 + 0.1486i; -0.1486 - 0.1486i; -414.69; -999.027]
[-0.0055; -0.0173 + 0.0179i; -0.0173 - 0.0179i; 
 -7.212 + 17.415i; -7.212 - 17.415i; -17.415 + 7.212i; -17.415 - 7.212i;
 -23.4459; -37.0379; -45.4483; -8.7307 + 42.4144i; -8.7307 - 42.4144i]
[-0.9211 + 0.94i; -0.9211 - 0.94i]
[-981 + 1009i; -981 - 1009i; -3290 + 1263i; -3290 - 1263i]};

% fix bad sac polezeros
for a=1:numel(badnet)
    if(~isfield(sacpzdb,badnet{a})); continue; end
    fprintf(badnet{a})
    cnt=0;
    for b=1:numel(sacpzdb.(badnet{a}).k)
        if(isequal(sacpzdb.(badnet{a}).(badfield{a})(b),badvalue(a)))
            % found one!
            cnt=cnt+1;
            sacpzdb.(badnet{a}).(badfield{a})(b)=goodvalue(a);
        end
    end
    fprintf([' --> Fixed ' num2str(cnt) ' SAC PoleZero(s)!\n'])
end

end
