function [diff]=timediff(times1,times2,option)
%TIMEDIFF    Return number of seconds between times
%
%    Usage:    diff=timediff(times1,times2)
%              diff=timediff(times1,times2,'utc'|'tai')
%
%    Description:
%     TIMEDIFF(TIMES1,TIMES2) returns the difference in time between times
%     TIMES1 and TIMES2 in units of seconds.  The formula is essentially
%     TIMES2-TIMES1.  So if TIMES2 is later than TIMES1 the number will be
%     positive.  TIMES1 and TIMES2 must be a Nx2, Nx3, Nx5, or Nx6 array of
%     [yr dayofyr], [yr mon dayofmon], [yr dayofyr hr min sec] or
%     [yr mon dayofmon hr min sec].  Only the seconds portion of TIMES1 and
%     TIMES2 is allowed to be non-integer (ie. you cannot have 1.5 minutes
%     etc).  TIMES1 and TIMES2 must have an equal number of times or be a
%     single time.
%
%     TIMEDIFF(TIMES1,TIMES2,'UTC'|'TAI') allows finding the difference
%     between UTC times (which may have leap seconds occasionally inserted
%     on certain dates -- see LEAPSECONDS).  The default is 'TAI' and does
%     not account for UTC leap seconds.  Using '' or [] will also give the
%     default behavior.
%
%    Notes:
%
%    Examples:
%     % Find the number of seconds in 2005:
%     timediff([2005 1],[2006 1],'utc')
%
%    See also: FIXTIMES, ISLEAPYEAR, LEAPSECONDS, UTC_OFFSET, UTC_LOD,
%              LEAPSECONDS_UPDATE, CAL2DOY, DOY2CAL, FIXDATES, UTC2TAI,
%              TAI2UTC, GREGORIAN2MODSERIAL, GREGORIAN2SERIAL,
%              SERIAL2GREGORIAN, MODSERIAL2GREGORIAN

%     Version History:
%        Nov. 12, 2008 - initial version
%        Apr. 23, 2009 - fix nargchk for octave, move usage up
%        June 10, 2009 - minor doc fix
%        June 24, 2009 - added scalar datetime expansion, minor doc update
%        Sep.  5, 2009 - added SUBMAT as a subfunction so the time package
%                        is free from 3rd party dependencies (I think),
%                        minor doc update
%        Sep. 18, 2009 - relaxed checks 
%                         + returns empty matrix on empty input
%                         + allow more cases with differing formats
%        Sep. 25, 2009 - description update
%        Feb. 11, 2011 - mass nargchk fix
%        Nov.  1, 2011 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov.  1, 2011 at 15:05 GMT

% todo:

% check nargin
error(nargchk(2,3,nargin));

% check times
sz1=size(times1);
sz2=size(times2);
if(isempty(times1) || isempty(times2)); diff=[]; return; end
if(~isnumeric(times1) || ~isnumeric(times2)...
        || ~any(sz1(2)==[2 3 5 6]) || ~any(sz2(2)==[2 3 5 6])...
        || (~isequal(sz1([1 3:end]),sz2([1 3:end])) && ...
        prod(sz1([1 3:end]))~=1 && prod(sz2([1 3:end]))~=1))
    error('seizmo:timediff:badInput',...
        'TIMES1 and TIMES2 must be numeric date arrays!');
end

% expand scalar
if(prod(sz1([1 3:end]))==1)
    times1=repmat(times1,[sz2(1) 1 sz2(3:end)]);
elseif(prod(sz2([1 3:end]))==1)
    times2=repmat(times2,[sz1(1) 1 sz1(3:end)]);
end

% check option
if(nargin==2 || isempty(option))
    option='tai';
elseif(~ischar(option) || ~any(strcmpi(option,{'utc' 'tai'})))
    error('seizmo:timediff:optionBad',...
        'OPTION must be ''utc'' or ''tai''!');
end

% proceed by option
switch lower(option)
    case 'tai'
        % TAI => MODIFIED SERIAL
        modserial=gregorian2modserial(times2)-gregorian2modserial(times1);
    case 'utc'
        % add zeros (date => time)
        if(any(sz1(2)==[2 3])); times1(:,end+3,:)=0; end
        if(any(sz2(2)==[2 3])); times2(:,end+3,:)=0; end
        
        % UTC => TAI => MODIFIED SERIAL
        modserial=gregorian2modserial(utc2tai(times2))...
            -gregorian2modserial(utc2tai(times1));
end

% take difference
diff=submat(modserial,2,1)*86400+submat(modserial,2,2);

end


function [X]=submat(X,varargin)
%SUBMAT    Returns a submatrix reduced along indicated dimensions
%
%    Usage:    Y=submat(X,DIM1,LIST1,DIM2,LIST2,...)
%
%    Description:
%     Y=SUBMAT(X,DIM,LIST) creates a matrix Y that is the matrix X reduced
%     along dimension DIM to the indices in LIST.  If DIM is a list of
%     dimensions, LIST is used to reduce each dimension.  DIM may be ':' to
%     indicate all dimensions of X (up to the maximum non-singleton).
%
%     Y=SUBMAT(X,DIM1,LIST1,DIM2,LIST2,...) allows for access to
%     multiple dimensions independently.  Later inputs take indexing
%     preference if a dimension is indexed more than once.
%
%    Notes:
%
%    Examples:
%     % Return x reduced to only the elements in index 1 of dimension 5:
%     x=submat(x,5,1)
%
%     % Reduce dimensions 1 thru 3 to return a matrix that only contains
%     % elements that were in row/column/page 4 or 5 for those dimensions:
%     x=submat(x,1:3,4:5)
%
%     % Reduce to elements in the 3rd row, 4th column:
%     x=submat(x,1,3,2,4)
%
%     % Return the element with all subscripts of 2:
%     x=submat(x,':',2)
%
%    See also: SUBMAT_EVAL, COLON OPERATOR (:), REPMAT

%     Version History:
%        Nov. 12, 2008 - initial version
%        Apr. 23, 2009 - move usage up
%        Sep.  8, 2009 - fixed error message
%        Sep. 21, 2009 - updated examples, removed unnecessary brackets
%        Aug.  9, 2010 - ':' dimension now functions, 15% slower though
%        Nov.  1, 2011 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov.  1, 2011 at 20:45 GMT

% todo:

% CHECK VARARGIN
if(~mod(nargin,2))
    error('misc:submat:badNumArgs','Unpaired DIMENSION,INDICES!');
end

% FIX ':' DIMENSION
varargin(2.*find(strcmp(':',varargin(1:2:end)))-1)={1:ndims(X)};

% DEFAULT TO ENTIRE MATRIX AND EXPAND TO MAX INPUT DIMENSION
list(1:max([ndims(X) [varargin{1:2:end}]]))={':'};

% REDUCED/REPLICATED DIMENSIONS
for i=1:2:nargin-2
    [list{varargin{i}}]=deal(varargin{i+1});
end

% SUBSET
X=X(list{:});

end
