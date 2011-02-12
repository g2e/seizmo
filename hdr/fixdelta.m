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
%     integers, as it extends the precision significantly.  NOTE THAT
%     EXTENDING THE PRECISION IS USELESS IF YOU DESTROY THE ACCURACY IN THE
%     PROCESS (IE DON'T CHANGE THE SAMPLE RATE TO A VALUE NOT FAITHFUL TO
%     THE DATA).  See the Notes section below for more ranting.
%
%     DATA=FIXDELTA(DATA,TOL) allows specifying the maximum tolerance TOL
%     that the fraction of 2 small integers must match DELTA within.  For
%     example, a tolerance of 1e-4 requires the new sample spacing to be
%     within 0.01% of the original sample spacing.  The default tolerance
%     is 1e-6 (0.0001% of the original sample spacing).  See the Notes
%     section below to understand the consequences of changing TOL.
%
%    Notes:
%     - FIXDELTA adjusts the sample spacing slightly to remove inaccuracy
%       caused by storing it as a single precision value.  This also can be
%       used to 'clean up' sample rates that can not be expressed nicely
%       (like 40sps, etc) but THIS DOES LEAD TO LARGE TIMING ERRORS.
%       For instance:
%        TOL=1e-4 on a day file   = up to +/-8.64s of error
%        TOL=1e-5 on a day file   = up to +/-.864s of error
%        TOL=1e-6 on a day file   = up to +/-.0864s of error (default)
%        TOL=1e-4 on an hour file = up to +/-.36s of error
%        TOL=1e-5 on an hour file = up to +/-.036s of error
%        TOL=1e-6 on an hour file = up to +/-.0036s of error (default)
%       Use FIXDELTA wisely by setting TOL to a level that keeps the error
%       introduced to an acceptable level.  SYNCRATES & INTERPOLATE may
%       be used to modify the sample spacing while getting new data values
%       for each of the new time points.
%
%    Header changes: DELTA, E
%
%    Examples:
%     Force double precision and update the delta field:
%      data=fixdelta(changeclass(data,'double'));
%
%    See also: RRAT, CHANGECLASS

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
%        Mar. 24, 2010 - added rants to docs, default tol set to 1e-6
%        Feb. 11, 2011 - mass nargchk fix, point to rrat not rat, todo list
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 11, 2011 at 15:05 GMT

% todo:
% - option to set tolerance as max time skew

% check nargin
error(nargchk(1,2,nargin));

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
        [n,d]=rrat(delta,1e-6);
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
    error(lasterror);
end

end
