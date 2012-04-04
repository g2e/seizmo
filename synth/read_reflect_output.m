function [raw]=read_reflect_output(file)
%READ_REFLECT_OUTPUT    Reads ascii reflectivity output file (DEPRECATED!)
%
%    Usage:    raw=read_reflect_output(file)
%
%    Description:
%     RAW=READ_REFLECT_OUTPUT(FILE) reads an ascii reflectivity output file
%     produced by running Brian Kennett's reflect code.  This output file
%     usually ends with the extension '.sei' and contains the synthetic
%     seismogram data concatenated together.  This will return the values
%     as a single vector of numbers.  Use the output of READ_REFLECT_INPUT
%     to parse this or just use REFLECT2SEIZMO to have it done for you.
%
%    Notes:
%     - The included reflect codes no longer produce ascii output so this
%       function is deprecated.  Why?  Reading ascii into matlab is slow.
%
%    Examples:
%
%    See also: READ_REFLECT_INPUT, MAKE_REFLECT_INPUT, REFLECT2SEIZMO

%     Version History:
%        Aug. 10, 2010 - initial version
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 23:00 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% get raw unparsed data
if(nargin<1); file=[]; end
raw=str2double(getwords(readtxt(file,...
    {'*.sei;*.SEI' 'Reflect Output Files (*.sei,*.SEI)';
    '*.*' 'All Files (*.*)'})));

end
