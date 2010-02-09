function [i1,i2,i3]=horzpairs(data)
%HORZPAIRS    Returns indice arrays for pairing horizontal SEIZMO records

% - first check for non-horizontals
%   - cmpinc == 90
%   - how can we skip?
%
% - next group by stream
%   - getstreamidx
%
% - then check that there are no more than 2 horz channels per stream
%   - error if we find more
%
% - then check that there is only 1 orientation per channel
%   - error if we find more
%
% - if only 1 channel in stream
%   - how can we skip?
%
% - make sure 2 channels are orthogonal
%   - error if not
%
% - ok now return this info in a useful manner
%   i1 = index in data
%   i2 = pair number (stream)
%   i3 = cmp number

end
