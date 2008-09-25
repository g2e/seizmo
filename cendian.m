function [data]=cendian(data,endianness)
%CENDIAN    Change the endianness of SAClab data records
%
%    Description: CENDIAN(DATA,ENDIANNESS) changes the byte-order that
%     records in the SAClab struct DATA will be written in to ENDIANNESS.
%     ENDIANNESS must be the string 'ieee-le' or 'ieee-be' or it may be a
%     char/cellstr array with rows/elements defining the endianness for 
%     each record.
%
%    Notes:
%
%    System requirements: Matlab 7
%
%    Input/Output requirements: DATA must be a valid SAClab struct,
%     ENDIANNESS must be a char or cellstr array
%
%    Header changes: N/A
%
%    Usage:    data=cendian(data,endianness)
%
%    Examples:
%     Change SAClab files in current directory to current platform's
%     byte-ordering:
%      wseis(cendian(rseis('*'),nativeendian))
%
%    See also: nativeendian, wseis, rseis, wh, rh, rdata, rpdw

%     Version History:
%        Sep. 25, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 25, 2008 at 06:20 GMT

% todo:

% check number of inputs
error(nargchk(2,2,nargin))

% check data structure
error(seischk(data))

% check and fix type
if(~iscellstr(endianness))
    if(~ischar(endianness))
        error('SAClab:cendian:badInput',...
            'ENDIANNESS must be a char or cellstr array!');
    else
        endianness=cellstr(endianness);
    end
end

% expand scalar
if(isscalar(endianness))
    endianness=endianness(ones(numel(data),1));
elseif(numel(data)~=numel(endianness))
    error('SAClab:cendian:badInput',...
        'ENDIANNESS must be scalar or match the size of DATA!');
end

% check endianness
endianness=lower(endianness);
if(any(~strcmp(endianness,'ieee-le') & ~strcmp(endianness,'ieee-be')))
    error('SAClab:cendian:badEndian',...
        'ENDIANNESS must be ieee-le or ieee-be!');
end

% change endianness
[data.endian]=deal(endianness{:});

end