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
%    See also: strtok

%     Version History:
%        June 11, 2009 - initial version
%
%     Testing Table:
%                                  Linux    Windows     Mac
%        Matlab 7       r14        
%               7.0.1   r14sp1
%               7.0.4   r14sp2
%               7.1     r14sp3
%               7.2     r2006a
%               7.3     r2006b
%               7.4     r2007a
%               7.5     r2007b
%               7.6     r2008a
%               7.7     r2008b
%               7.8     r2009a
%        Octave 3.2.0
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June 11, 2009 at 09:55 GMT

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
