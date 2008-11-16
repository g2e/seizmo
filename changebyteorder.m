function [data]=changebyteorder(data,endianness)
%CHANGEBYTEORDER    Change the byteorder of SEIZMO data records
%
%    Description: CHANGEBYTEORDER(DATA,ENDIANNESS) changes the byte-order
%     that the records in the SEIZMO struct DATA will be written as to
%     ENDIANNESS.  ENDIANNESS must be the string 'ieee-le' or 'ieee-be' or
%     it may be a char/cellstr array of those strings to define each
%     record's endianness individually.
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Usage:    data=changebyteorder(data,endianness)
%
%    Examples:
%     Change records in current directory to the platform's byte-ordering:
%      writeseizmo(changebyteorder(readseizmo('*'),nativebyteorder))
%
%    See also: nativebyteorder, writeseizmo, readseizmo, bseizmo,
%              readdata, readdatawindow, readheader, writeheader

%     Version History:
%        Sep. 25, 2008 - initial version
%        Nov. 16, 2008 - rename from CENDIAN to CHANGEBYTEORDER
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov. 16, 2008 at 06:30 GMT

% todo:

% check number of inputs
error(nargchk(2,2,nargin))

% check data structure
error(seizmocheck(data))

% check and fix type
if(~iscellstr(endianness))
    if(~ischar(endianness))
        error('seizmo:changebyteorder:badInput',...
            'ENDIANNESS must be a char or cellstr array!');
    else
        endianness=cellstr(endianness);
    end
end

% expand scalar
if(isscalar(endianness))
    endianness=endianness(ones(numel(data),1));
elseif(numel(data)~=numel(endianness))
    error('seizmo:cendian:badInput',...
        'ENDIANNESS must be scalar or match the size of DATA!');
end

% check endianness
endianness=lower(endianness);
if(any(~strcmp(endianness,'ieee-le') & ~strcmp(endianness,'ieee-be')))
    error('seizmo:changebyteorder:badEndian',...
        'ENDIANNESS must be ''ieee-le'' or ''ieee-be''!');
end

% change endianness
[data.endian]=deal(endianness{:});

end