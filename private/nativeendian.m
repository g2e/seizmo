function [endianness]=nativeendian()
%NATIVEENDIAN    Returns native endianness of present platform
%
%    Description: NATIVEENDIAN returns the endianness of the current
%     computer using the builtin Matlab function COMPUTER.  Output is a
%     string of either 'ieee-le' or 'ieee-be' (works with FREAD/FWRITE
%     statements).
%
%    Notes:
%
%    Tested on: Matlab r2007b
%
%    Header changes: N/A
%
%    Usage:    endianness=nativeendian
%
%    Examples:
%     Change SAClab files in current directory to current platform's
%     byte-ordering:
%      wseis(cendian(rseis('*'),nativeendian))
%
%    See also: cendian, wseis, rseis, wh, rh, rdata, rpdw

%     Version History:
%        Sep. 25, 2008 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 25, 2008 at 05:45 GMT

% todo:

% get this platform's native byte-order
[platform,maxint,endianness]=computer;
clear platform maxint

% put in a SAClab consumable format (for fread/fwrite)
if(strcmpi(endianness,'L')); endianness='ieee-le';
else endianness='ieee-be'; end

end