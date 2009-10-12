function [words]=getwords(str,delimiter)
%GETWORDS    Returns a cell array of words from a string
%
%    Usage:    words=getwords('str')
%              words=getwords('str',delimiter)
%
%    Description: WORDS=GETWORDS('STR') extracts words in STR and returns
%     them separated into a cellstr array WORDS without any whitespace.
%
%     WORDS=GETWORDS('STR',DELIMITER) separates words in STR using the
%     single character DELIMITER.
%
%    Notes:
%
%    Examples:
%     Break up a sentence:
%      getwords('This example is pretty dumb!')
%       ans = 
%       'This'    'example'    'is'    'pretty'    'dumb!'
%
%    See also: JOINWORDS, STRTOK, ISSPACE

%     Version History:
%        June 11, 2009 - initial version
%        Sep. 13, 2009 - minor doc update, added input check
%        Sep. 16, 2009 - add delimiter option
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 16, 2009 at 03:20 GMT

% todo:

% check nargin
msg=nargchk(1,2,nargin);
if(~isempty(msg)); error(msg); end;

% check str
if(~ischar(str) || ~isvector(str))
    error('seizmo:getwords:badInput','STR must be a char array!');
end

% force str to row vector
str=str(:).';

% highlight
if(nargin==2)
    % check delimiter
    if(~ischar(delimiter) || ~isscalar(delimiter))
        error('seizmo:getwords:badInput','DELIMITER must be a char!');
    end
    
    % highlight word boundaries
    idx=diff([false str~=delimiter false]);
else
    % highlight word boundaries
    idx=diff([false ~isspace(str) false]);
end

% get word boundaries
s=find(idx==1);
e=find(idx==-1)-1;

% number of words
nw=numel(s);

% get words
words=cell(1,nw);
for i=1:nw; words{i}=str(s(i):e(i)); end

end
