function [dates,leaps]=getleaps(option)
%GETLEAPS    Returns leapsecond info and stores in memory for later calls
%
%    Description:  [DATES,LEAPS]=GETLEAPS() will return leap second info
%     and store that info in a global variable to speed up subsequent
%     calls.  DATES contains the serial dates on which the 'leaps' occur
%     and LEAPS contains the offset (either 1 or -1).
%
%     [DATES,LEAPS]=GETLEAPS(OPTION) controls the behavior of GETLEAPS.
%     OPTION must be a string equal to either 'REFRESH' or 'REUSE'.
%     REFRESH will recalculate the leapsecond info by calling LEAPSECONDS.
%     REUSE will use the leap info in memory if possible, otherwise
%     calculating the leap info using LEAPSECONDS.  The default is 'REUSE'.
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Usage:    [dates,leaps]=getleaps()
%              [dates,leaps]=getleaps('refresh'|'reuse')
%
%    Examples:
%     To see the benefit of 'reuse':
%      tic; for i=1:100; [dates,leaps]=getleaps('reuse'); end; toc
%      tic; for i=1:100; [dates,leaps]=getleaps('refresh'); end; toc
%
%    See also: leapseconds, totalleaps, leapsinday

%     Version History:
%        Nov. 10, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 10, 2008 at 16:00 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% check option
if(nargin==0 || isempty(option))
    option='reuse';
elseif(~ischar(option) ...
        || (~strcmpi(option,'refresh') && ~strcmpi(option,'reuse')))
    error('SAClab:getleaps:optionBad',...
        'Bad OPTION -- must be ''refresh'' or ''reuse''');
end

% global
global SACLAB

% retrieve from memory or recalculate
switch lower(option)
    case 'refresh'
        leapstr=leapseconds;
        dates=datenum(leapstr);
        leaps=2*strcmp('+',cellstr(leapstr(:,end)))-1;
        [dates,idx]=sort(dates);
        leaps=leaps(idx);
        SACLAB.GETLEAPS.DATES=dates;
        SACLAB.GETLEAPS.LEAPS=leaps;
    case 'reuse'
        if(isfield(SACLAB,'GETLEAPS')...
                && isfield(SACLAB.GETLEAPS,'DATES')...
                && isfield(SACLAB.GETLEAPS,'LEAPS')...
                && isnumeric(SACLAB.GETLEAPS.DATES)...
                && isnumeric(SACLAB.GETLEAPS.LEAPS)...
                && isequal(size(SACLAB.GETLEAPS.LEAPS),...
                size(SACLAB.GETLEAPS.DATES)))
            dates=SACLAB.GETLEAPS.DATES;
            leaps=SACLAB.GETLEAPS.LEAPS;
        else
            [dates,leaps]=getleaps('refresh');
        end
end
