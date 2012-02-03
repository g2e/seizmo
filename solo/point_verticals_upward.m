function [data,down]=point_verticals_upward(data)
%POINT_VERTICALS_UPWARD    Adjusts vertically downward SEIZMO records
%
%    Usage:    data=point_verticals_upward(data)
%              [data,down]=point_verticals_upward(data)
%
%    Description:
%     DATA=POINT_VERTICALS_UPWARD(DATA) will reorient verticals pointed
%     downward to upward.  This is done by flipping the polarity of the
%     data in the records and setting the CMPINC field from 180/-180 to 0.
%     All records with CMPINC not set to 180/-180 are unaltered.
%
%     [DATA,DOWN]=POINT_VERTICALS_UPWARD(DATA) also returns the logical
%     matrix DOWN which is equal sized to DATA with TRUE elements
%     indicating the records adjusted.
%
%    Notes:
%     - Only checks that the CMPINC field is 180/-180 to assess orientation
%
%    Header changes: DEP*, CMPINC
%
%    Examples:
%     % See if you have any downward pointing data:
%     plot2(data,point_verticals_upward(data));
%
%    See also: VERTCMP, ROTATE

%     Version History:
%        Jan. 12, 2011 - initial version
%        Jan. 28, 2012 - seizmocheck state bugfix, doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 28, 2012 at 13:35 GMT

% check nargin
error(nargchk(1,1,nargin));

% check data structure
error(seizmocheck(data,'dep'));

% turn off checking
oldseizmocheckstate=seizmocheck_state(false);

% attempt reorientation
try
    % get orientation
    cmpinc=getheader(data,'cmpinc');
    
    % find vertically down (either -180 or 180)
    down=abs(cmpinc)==180;
    
    % only proceed if downward exists
    if(sum(down))
        % flip data polarity
        data(down)=multiply(data(down),-1);
        
        % set orientation to up
        data(down)=changeheader(data(down),'cmpinc',0);
    end
    
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
catch
    % toggle checking back
    seizmocheck_state(oldseizmocheckstate);
    
    % rethrow error
    error(lasterror);
end

end
