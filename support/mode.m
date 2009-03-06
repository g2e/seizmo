function m = mode(x)
%
% MODE statistical mode of an array
%   MODE(X)  returns the most frequent element(s) of array X. X can be a
%   charracter array or a cell array of strings.
%
%   Example      
%     mode([3 9 4 3 2 3 2])  returns 3,
%     mode([3 9 4 2 3 2]) returns [2;3].
%
%   See also FREQTABLE

% Mukhtar Ullah
% December 28, 2004
% mukhtar.ullah@informatik.uni-rostock.de

[y,h] = freqtable(x);
m = y(h==max(h));