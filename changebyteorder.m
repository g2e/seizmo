function [data]=changebyteorder(data,endianness)
%CHANGEBYTEORDER    Change the byteorder of SEIZMO data records
%
%    Usage:    data=changebyteorder(data,endianness)
%
%    Description: CHANGEBYTEORDER(DATA,ENDIANNESS) changes the byte-order
%     that the records in the SEIZMO struct DATA will be written as to
%     ENDIANNESS.  ENDIANNESS must be the string 'ieee-le' or 'ieee-be' or
%     it may be a char/cellstr array of those strings to define each
%     record's endianness individually.
%
%    Notes:
%     - empty string ('') will preserve the byteorder of records
%
%    Examples:
%     Change records in current directory to the platform's byte-ordering:
%      writeseizmo(changebyteorder(readseizmo('*'),nativebyteorder))
%
%    See also: NATIVEBYTEORDER, WRITESEIZMO, READSEIZMO, BSEIZMO,
%              READDATA, READDATAWINDOW, READHEADER, WRITEHEADER

%     Version History:
%        Sep. 25, 2008 - initial version
%        Nov. 16, 2008 - rename from CENDIAN to CHANGEBYTEORDER
%        Apr. 23, 2009 - fix nargchk and seizmocheck for octave,
%                        move usage up
%        May  28, 2009 - fixed to update byteorder field, minor doc update
%        May  29, 2009 - allow empty endianness (no change)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 17, 2009 at 20:45 GMT

% todo:

% check number of inputs
msg=nargchk(2,2,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data);
if(~isempty(msg)); error(msg.identifier,msg.message); end

% fast exit
if(isempty(endianness)); return; end

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
    error('seizmo:changebyteorder:badInput',...
        'ENDIANNESS must be scalar or match the size of DATA!');
end

% find those to preserve
preserve=strcmp(endianness,'');
[endianness{preserve}]=deal(data(preserve).byteorder);

% check endianness
endianness=lower(endianness);
if(any(~strcmp(endianness,'ieee-le') & ~strcmp(endianness,'ieee-be')))
    error('seizmo:changebyteorder:badEndian',...
        'ENDIANNESS must be ''ieee-le'' or ''ieee-be''!');
end

% change endianness
[data.byteorder]=deal(endianness{:});

end
