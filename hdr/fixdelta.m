function [data]=fixdelta(data,tol,units)
%FIXDELTA    Fix sample spacing for SEIZMO records
%
%    Usage:    data=fixdelta(data)
%              data=fixdelta(data,tol)
%              data=fixdelta(data,tol,units)
%
%    Description:
%     DATA=FIXDELTA(DATA) modifies the sample spacing (DELTA header field)
%     of records in DATA to be the decimal equivalent of a fraction of 2
%     small integers.  This is particularly useful for upgrading the sample
%     spacing from single to double precision when the original sample
%     spacing can be expressed as the fraction of two small integers, as it
%     extends the precision significantly.  NOTE THAT EXTENDING THE
%     PRECISION IS USELESS IF YOU DESTROY THE ACCURACY IN THE PROCESS (IE
%     DON'T CHANGE THE SAMPLE RATE TO A VALUE NOT FAITHFUL TO THE DATA).
%     See the Notes section below for more ranting.
%
%     DATA=FIXDELTA(DATA,TOL) specifies the maximum skew in record length
%     allowed by changing the sample spacing.  The default is 0.01s.  See
%     the Notes section below to understand the consequences of changing
%     TOL.
%
%     DATA=FIXDELTA(DATA,TOL,UNITS) changes the units of TOL.  UNITS can be
%     either 'ratio', 'abs' or 'sec'.  The default is 'sec'.
%      RATIO - New DELTA must be within +/- TOL*DELTA of old DELTA
%        ABS - New DELTA must be within +/- TOL of old DELTA
%        SEC - New DELTA cannot skew the record length more than +/- TOL
%
%    Notes:
%     - FIXDELTA adjusts the sample spacing slightly to remove inaccuracy
%       caused by storing it as a single precision value.  This also can be
%       used to 'clean up' sample rates that can not be expressed nicely
%       (like 40sps, etc) but THIS DOES LEAD TO LARGE TIMING ERRORS.  Use
%       FIXDELTA wisely by setting TOL to a level that keeps the error
%       introduced to an acceptable level.  SYNCRATES & INTERPOLATE may
%       be used to modify the sample spacing while getting new data values
%       for each of the new time points.
%
%    Header changes: DELTA, E
%
%    Examples:
%     % Force double precision and update the delta field:
%     data=fixdelta(changeclass(data,'double'));
%
%     % Require all samples to be within 1msec of their original timing:
%     data=fixdelta(data,.001,'s');
%
%    See also: RAT, RRAT, CHANGECLASS

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
%        Jan. 30, 2012 - units option, some debug code, change defaults
%        Feb.  4, 2012 - no warning for multiple delta
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  4, 2012 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% check data structure & header
data=checkheader(data,'multiple_delta','ignore');

% turn off struct checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt delta fix
try
    % grab header info
    [delta,b,npts,leven]=getheader(data,'delta','b','npts','leven lgc');
    
    % only work on evenly spaced
    es=strcmpi(leven,'false');
    if(seizmoverbose && any(es))
        warning('seizmo:fixdelta:Uneven',...
            ['Not fixing DELTA field for unevenly sampled record(s):' ...
            sprintf('%d ',find(es))]);
    end
    es=~es;

    % defaults
    if(nargin<2 || isempty(tol)); tol=.01; end
    if(nargin<3 || isempty(units)); units='sec'; end
    
    % check tol
    if(~isscalar(tol) || ~isreal(tol))
        error('seizmo:fixdelta:badInput',...
            'TOL must be a real-valued scalar!');
    elseif(~isscalar(cellstr(units)))
        error('seizmo:fixdelta:badInput',...
            'UNITS must be a string!');
    end
    
    % get fractions
    switch lower(units)
        case {'ratio' 'r'}
            [n,d]=rrat(delta,tol);
            if(seizmodebug)
                disp('OLD_DELTA  NEW_DELTA  SKEW  SKEW_ALLOWED');
                disp([delta n./d (delta-n./d).*(npts-1) delta.*tol.*npts]);
            end
        case {'abs' 'a' 'absolute'}
            [n,d]=rat(delta,tol);
            if(seizmodebug)
                disp('OLD_DELTA  NEW_DELTA  SKEW  SKEW_ALLOWED');
                disp([delta n./d (delta-n./d).*(npts-1) tol.*npts])
            end
        case {'sec' 's' 'seconds'}
            if(~isscalar(unique(npts)))
                nrecs=numel(data); n=nan(nrecs,1); d=n;
                for i=1:nrecs
                    [n(i),d(i)]=rat(delta(i),tol/npts(i));
                end
            else
                [n,d]=rat(delta,tol/unique(npts));
            end
            if(seizmodebug)
                disp('OLD_DELTA  NEW_DELTA  SKEW  SKEW_ALLOWED');
                disp([delta n./d (delta-n./d).*(npts-1) tol.*n.^0])
            end
        otherwise
            error('seizmo:fixdelta:badInput',...
                'Unknown UNITS: %s',units);
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
