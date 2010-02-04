function [data]=fixdelta(data,tol)
%FIXDELTA    Fix sample spacing for SEIZMO records
%
%    Usage:    data=fixdelta(data)
%              data=fixdelta(data,tol)
%
%    Description: DATA=FIXDELTA(DATA) modifies the sample spacing (DELTA
%     header field) of records in DATA to be the decimal equivalent of a
%     fraction of 2 small integers.  This is particularly useful for
%     upgrading the sample spacing from single to double precision when the
%     original sample spacing can be expressed as the fraction of two small
%     integers, as it extends the precision significantly.
%
%     DATA=FIXDELTA(DATA,TOL) allows specifying the maximum tolerance TOL
%     that the fraction of 2 small integers must match DELTA within.  For
%     example, a tolerance of 1e-4 requires the new sample spacing to be
%     within 0.01% of the original sample spacing.  The default tolerance
%     is 1e-5 (0.001% of the original sample spacing).
%
%    Notes:
%
%    Header changes: DELTA, E
%
%    Examples:
%     Force double precision and update the delta field:
%      data=fixdelta(changeclass(data,'double'));
%
%    See also: RAT, CHANGECLASS

%     Version History:
%        Feb. 21, 2008 - initial version
%        Feb. 23, 2008 - minor doc update
%        Feb. 28, 2008 - minor doc update
%        Mar.  4, 2008 - minor doc update
%        Oct.  8, 2008 - doc update, add history
%        Nov. 16, 2008 - update for new name schema, doc and history update
%        Nov. 22, 2008 - doc update
%        Dec.  5, 2008 - tolerance option added
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        Oct.  5, 2009 - doc update
%        Nov. 26, 2009 - document RAT bug, change default to not vary with
%                        input (RAT's default does so we choose a default)
%        Dec.  4, 2009 - update E, only fix evenly spaced
%        Jan. 25, 2010 - fixed LEVEN bug (only fixed first record)
%        Jan. 29, 2010 - minor code cleaning
%        Jan. 30, 2010 - one less call to SEIZMOCHECK, use RRAT
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 30, 2010 at 19:45 GMT

% todo:

% check nargin
msg=nargchk(1,2,nargin);
if(~isempty(msg)); error(msg); end

% check data structure & header
data=checkheader(data);

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt delta fix
try
    % grab header info
    leven=getlgc(data,'leven');
    [delta,b,npts]=getheader(data,'delta','b','npts');
    
    % only work on evenly spaced
    es=strcmpi(leven,'false');
    if(seizmoverbose && any(es))
        warning('seizmo:fixdelta:Uneven',...
            ['Record(s):\n' sprintf('%d ',find(es)) ...
            '\nNot fixing DELTA field for unevenly sampled record(s)!']);
    end
    es=~es;

    % default tolerance
    if(nargin==1 || ~isscalar(tol) || ~isreal(tol))
        [n,d]=rrat(delta,1e-5);
    else
        [n,d]=rrat(delta,tol);
    end

    % fix delta
    data(es)=changeheader(data(es),'delta',n(es)./d(es),...
        'e',b(es)+(npts(es)-1).*n(es)./d(es));

    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

end
