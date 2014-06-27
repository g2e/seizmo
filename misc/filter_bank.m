function [bank]=filter_bank(range,option,width,offset)
%FILTER_BANK    Makes a set of narrow-band bandpass filters
%
%    Usage:    bank=filter_bank(range,option,width,offset)
%
%    Description:
%     BANK-FILTER_BANK(RANGE,OPTION,WIDTH,OFFSET) returns a set (also known
%     as a "bank") of bandpass filters.  Each row in the bank corresponds
%     to a separate filter, with the first column giving the center
%     frequency and the 2nd and 3rd columns giving the low- and high-
%     frequency passband corners for that filter:
%        BANK=[Fc Flo Fhi; Fc Flo Fhi; ...]
%     RANGE is a two element vector giving the frequency range for the
%     filter bank in Hz.  OPTION is either 'constant' or 'variable' and 
%     indicates whether the width and offset for filters in the bank are 
%     constant over the frequency range (WIDTH and OFFSET are taken in Hz)
%     or if the width and offset varies (WIDTH and OFFSET are assumed to be
%     given as a fraction of the center frequency).
%
%    Notes:
%     - The 'variable' case gives a constant width and offset in
%       logarithmic space.
%
%    Examples:
%     % Build a filter bank over the range 0.01-0.1 Hz with filter widths
%     % equal to 20% the center frequency and adjacent filters separated by 
%     % 10% of the larger center frequency:
%     bank=filter_bank([0.01 0.1],'variable',0.2,0.1)
%
%     % Build a filter bank over the range 0.01-0.1 Hz with filter widths
%     % of constant 10 mHz and offset by 5 mHz:
%     bank=filter_bank([0.01 0.1],'constant',0.010,0.005)
%
%     % You can also build banks based on desired periods!  To build a
%     % filter bank over the range 18s-182s with filter widths 1s wide and
%     % stepping at 1s increments:
%     bank=1./filter_bank([18 182],'constant',1,1)
%
%    See also: IIRFILTER, MULTIBANDALIGN

%     Version History:
%        Mar.  6, 2009 - initial version
%        Apr. 23, 2009 - octave compat fix, minor doc fix
%        Sep.  7, 2009 - doc update
%        Feb.  5, 2010 - fixed constant width filter bug, added note
%        Feb. 11, 2011 - mass nargchk fix
%        Mar.  2, 2011 - doc update
%        June  4, 2014 - doc update, added example with period ranges after
%                        discussion with Ghassan
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June  4, 2014 at 15:05 GMT

% todo:

% check number of arguments
error(nargchk(4,4,nargin));

% check arguments
if(numel(range)>2 || ~isnumeric(range) || any(range<0))
    error('seizmo:filter_bank:badInput',...
        'RANGE must be an array of 2 positive frequencies (Hz)')
end
if(~isscalar(width) || ~isnumeric(width) || width<0)
    error('seizmo:filter_bank:badInput','Filter width must be >0')
end
if(~isscalar(offset) || ~isnumeric(offset) || offset<0)
    error('seizmo:filter_bank:badInput','Filter offset must be >0')
end

% fix range
range=sort(range(:));

% decide how to make bank based on option
if(isequal(option,'constant'))
    bank(:,1)=range(1):offset:range(2);
    bank(:,2)=bank(:,1)-width/2;
    bank(:,3)=bank(:,1)+width/2;
elseif(isequal(option,'variable'))
    bank(1,1)=range(1);
    width1=1-width/2;
    width2=1+width/2;
    offset=1+offset;
    count=1;
    while(1)
        bank(count,2)=bank(count,1)*width1;
        bank(count,3)=bank(count,1)*width2;
        count=count+1;
        bank(count,1)=bank(count-1,1)*offset;
        if(bank(count,1)>range(2))
            bank(count,:)=[];
            break;
        end
    end
else
    error('seizmo:filter_bank:badInput','Unknown option: %s',option)
end

end
