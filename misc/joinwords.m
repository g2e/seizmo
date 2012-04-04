function [str]=joinwords(words,delimiter)
%JOINWORDS    Combines a cellstr into a space-separated string
%
%    Usage:    str=joinwords(words)
%              str=joinwords(words,delimiter)
%
%    Description:
%     STR=JOINWORDS(WORDS) combines the cellstr array WORDS and returns a
%     single rowed char array STR with the elements of WORDS delimited by
%     spaces.
%
%     STR=JOINWORDS(WORDS,DELIMITER) uses the string DELIMITER to separate
%     words in the output string.
%
%    Notes:
%     - Punctuation should be part of the preceeding word!
%
%    Examples:
%     % Make a sentence but with underscores:
%     joinwords({'This' 'example' 'is' 'pretty' 'dumb!'},'_')
%      ans = 
%      This_example_is_pretty_dumb!
%
%    See also: GETWORDS, STRTOK, ISSPACE

%     Version History:
%        Sep. 13, 2009 - initial version
%        Sep. 16, 2009 - add delimiter option
%        Aug. 14, 2010 - allow multi-character delimiter
%        Apr.  4, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  4, 2012 at 04:50 GMT

% todo:

% check words
if(~iscellstr(words))
    error('seizmo:joinwords:badInput','WORDS must be a cellstr array!');
end

% delimit based on nargin
if(nargin==2)
    % check delimiter
    if(~ischar(delimiter) || ndims(delimiter)~=2 || size(delimiter,1)~=1)
        error('seizmo:joinwords:badInput','DELIMITER must be a string!');
    end
    str=sprintf(['%s' delimiter],words{:});
    n=numel(delimiter);
else
    str=sprintf('%s ',words{:});
    n=1;
end
str=str(1:end-n);

end
