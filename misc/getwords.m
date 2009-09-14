function [words]=getwords(str)
%GETWORDS    Returns a cell array of words from a string
%
%    Usage:    words=getwords('str')
%
%    Description: WORDS=GETWORDS('STR') extracts words in STR and returns
%     them separated into a cellstr array WORDS without any whitespace.
%
%    Notes:
%     - punctuation is not handled (although it could be)
%
%    Examples:
%     Break up a sentence:
%      getwords('This example is pretty dumb!')
%       ans = 
%       'This'    'example'    'is'    'pretty'    'dumb!'
%
%    See also: joinwords, strtok, isspace

%     Version History:
%        June 11, 2009 - initial version
%        Sep. 13, 2009 - minor doc update, added input check
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 13, 2009 at 15:25 GMT

% todo:

% check str
if(~ischar(str))
    error('seizmo:getwords:badInput','STR must be a char array!');
end

% highlight word boundaries
idx=diff([false ~isspace(str) false]);

% get word boundaries
s=find(idx==1);
e=find(idx==-1)-1;

% number of words
nw=numel(s);

% get words
words=cell(1,nw);
for i=1:nw; words{i}=str(s(i):e(i)); end

end
