function [str]=joinwords(words,delimiter)
%JOINWORDS    Combines a cellstr into a space-separated string
%
%    Usage:    str=joinwords(words)
%              str=joinwords(words,delimiter)
%
%    Description: STR=JOINWORDS(WORDS) combines the cellstr array WORDS and
%     returns a single rowed char array STR with the elements of WORDS
%     delimited by spaces.
%
%     STR=JOINWORDS(WORDS,DELIMITER) uses the single character DELIMITER to
%     separate words in the output string.
%
%    Notes:
%     - punctuation should be part of the preceeding word!
%
%    Examples:
%     Make a sentence but with underscores:
%      joinwords({'This' 'example' 'is' 'pretty' 'dumb!'},'_')
%       ans = 
%       This_example_is_pretty_dumb!
%
%    See also: getwords, strtok, isspace

%     Version History:
%        Sep. 13, 2009 - initial version
%        Sep. 16, 2009 - add delimiter option
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 16, 2009 at 04:50 GMT

% todo:

% check words
if(~iscellstr(words))
    error('seizmo:joinwords:badInput','WORDS must be a cellstr array!');
end

% delimit based on nargin
if(nargin==2)
    % check delimiter
    if(~ischar(delimiter) || ~isscalar(delimiter))
        error('seizmo:joinwords:badInput','DELIMITER must be a char!');
    end
    str=sprintf(['%s' delimiter],words{:});
else
    str=sprintf('%s ',words{:});
end
str=str(1:end-1);

end
