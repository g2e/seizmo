function [B,N] = freqtable(A)
%
% FREQTABLE Frequency table.
%   [Y,N] = FREQTABLE(X) takes a vector X and returns the unique values of
%   X in the output Y, and the number of instances of each value in the 
%   output N. X can be a charachter array or cell array of strings. 
%

% Mukhtar Ullah
% December 28, 2004
% mukhtar.ullah@informatik.uni-rostock.de

if isnumeric(A)     % use of built-in functions to avoid UNIQUE
    S = sort(A(~isnan(A)));
    B = S([find(diff(S));end]);
    N = histc(S,B);
else
    [B,m,n] = unique(A);
    N = histc(n,1:numel(B));
end