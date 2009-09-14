function [str]=joinwords(words)
%JOINWORDS    Combines a cellstr into a space-separated string
%
%    Usage:    str=joinwords(words)
%
%    Description: STR=JOINWORDS(WORDS) combines the cellstr array WORDS and
%     returns a single rowed char array STR with the elements of WORDS
%     delimited by spaces.
%
%    Notes:
%     - punctuation should be part of the preceeding word
%
%    Examples:
%     Make a sentence:
%      joinwords({'This' 'example' 'is' 'pretty' 'dumb!'})
%       ans = 
%       This example is pretty dumb!
%
%    See also: getwords, strtok, isspace

%     Version History:
%        Sep. 13, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 13, 2009 at 15:15 GMT

% todo:

% check words
if(~iscellstr(words))
    error('seizmo:joinwords:badInput','WORDS must be a cellstr array!');
end

str=sprintf('%s ',words{:});
str=str(1:end-1);

end
