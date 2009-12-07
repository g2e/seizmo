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
%     - Matlab r2007b function RAT has a bug in it - change line 116ish to:
%        if(x==0) || (abs((C(1,1)/C(2,1)-X(j))/X(j))<=max(tol,eps(X(j))))
%       This will force RAT to function as described.
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
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Dec.  4, 2009 at 06:35 GMT

% todo:

% check nargin
msg=nargchk(1,2,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data);
if(~isempty(msg)); error(msg.identifier,msg.message); end

% turn off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% attempt header check
try
    % check header
    data=checkheader(data);
    
    % turn off header checking
    oldcheckheaderstate=get_checkheader_state;
    set_checkheader_state(false);
catch
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror)
end

% attempt rest
try
    % grab header info
    [leven,delta,b,npts]=getheader(data,'leven','delta','b','npts');
    
    % only work on evenly spaced
    es=strcmpi(leven,'false');
    v=seizmoverbose;
    if(any(es) && v)
        warning('seizmo:fixdelta:Uneven',...
            ['Record(s):\n' sprintf('%d ',find(es)) ...
            '\nNot fixing DELTA field for unevenly sampled record(s)!']);
    end
    es=~es;

    % default tolerance
    if(nargin==1 || ~isscalar(tol) || ~isreal(tol))
        [n,d]=rat(delta,1e-5);
    else
        [n,d]=rat(delta,tol);
    end

    % fix delta
    data(es)=changeheader(data(es),'delta',n(es)./d(es),...
        'e',b(es)+(npts(es)-1).*n(es)./d(es));

    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    set_checkheader_state(oldcheckheaderstate);
catch
    % toggle checking back
    set_seizmocheck_state(oldseizmocheckstate);
    set_checkheader_state(oldcheckheaderstate);
    
    % rethrow error
    error(lasterror)
end

end
