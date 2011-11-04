function [words]=getwords(str,delimiter,collapse)
%GETWORDS    Returns a cell array of words from a string
%
%    Usage:    words=getwords('str')
%              words=getwords('str',delimiter)
%              words=getwords('str',delimiter,collapse)
%
%    Description:
%     WORDS=GETWORDS('STR') extracts words in STR and returns them
%     separated into a cellstr array WORDS without any whitespace.
%
%     WORDS=GETWORDS('STR',DELIMITER) separates words in STR using the
%     single character DELIMITER. The default is '' or [] which indicates
%     any whitespace character.
%
%     WORDS=GETWORDS('STR',DELIMITER,COLLAPSE) toggles treating multiple
%     delimiters as a single delimiter.  Setting COLLAPSE to TRUE treats
%     multiple delimiters between words as a single delimiter.  Setting
%     COLLAPSE to FALSE will always return the word between a delimiter
%     pair, even if the word is '' (ie no characters).  The default is
%     TRUE.
%
%    Notes:
%
%    Examples:
%     % Break up a sentence:
%     getwords('This example is pretty dumb!')
%       ans = 
%       'This'    'example'    'is'    'pretty'    'dumb!'
%
%     % Turn off multi-delimiter collapsing to allow handling empty words:
%     getwords('the answer is  !',[],false)
%       ans = 
%       'the'    'answer'    'is'    [1x0 char]    '!'
%
%    See also: JOINWORDS, STRTOK, ISSPACE

%     Version History:
%        June 11, 2009 - initial version
%        Sep. 13, 2009 - minor doc update, added input check
%        Sep. 16, 2009 - add delimiter option
%        Nov. 20, 2009 - make multi-delimiter collapse optional
%        July 30, 2010 - nargchk fix
%        Jan.  3, 2011 - use isstring
%        Nov.  1, 2011 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov.  1, 2011 at 13:00 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% check str
if(~isstring(str))
    error('seizmo:getwords:badInput','STR must be a char array!');
end

% check collapse
if(nargin<3 || isempty(collapse)); collapse=true; end
if(~islogical(collapse) || ~isscalar(collapse))
    error('seizmo:getwords:badInput','COLLAPSE must be a logical!');
end

% force str to row vector
str=str(:).';

% highlight word boundaries
if(nargin>1 && ~isempty(delimiter))
    % check delimiter
    if(~ischar(delimiter) || ~isscalar(delimiter))
        error('seizmo:getwords:badInput','DELIMITER must be a char!');
    end
    idx=[true str==delimiter true];
else
    idx=[true isspace(str) true];
end

% get word boundaries
if(collapse)
    idx=diff(idx);
    s=find(idx==-1);
    e=find(idx==1)-1;
else
    s=find(idx(1:end-1));
    e=find(idx(2:end))-1;
end

% number of words
nw=numel(s);

% get words
words=cell(1,nw);
for i=1:nw; words{i}=str(s(i):e(i)); end

end
