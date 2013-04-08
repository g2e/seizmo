function [int]=str2int(str)
%STR2INT    Converts simple integer strings to integers (fast)
%
%    Usage:    int=str2int(str)
%
%    Description:
%     INT=STR2INT(STR) converts the string STR which is all digits from 0-9
%     to the corresponding integer.  INT is a double array containing the
%     integers.  No input checking!.  For more details see the Notes & 
%     Examples sections.
%
%    Notes:
%     - STR2INT is a speedy, checkless string conversion function that
%       comes in handy for some parsing jobs.
%     - NO SPACES ALLOWED.  Convert any leading spaces to zeros.
%
%    Examples:
%     % Convert a simple number string to the number:
%     str2int('1234567890')
%
%     % Convert several integer strings (as a char array):
%     str2int(['7840';'0081';'0639'])
%
%    See also: STR2DOUBLE, STR2NUM

%     Version History:
%        Apr.  5, 2013 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  5, 2013 at 11:15 GMT

% todo:

str=double(str)-48;
pt=10.^(size(str,2)-1:-1:0)';
int=str*pt;

end
