function [data]=convolve(data,varargin)
%CONVOLVE    Convolve SEIZMO records with a function
%
%

% todo:
% - need options to create convolve function
%   - type (all windows?)
%     - triangle
%     - square
%     - gaussian
%     - custom
%   - halfduration
%   - maxamplitude
%   - offset
%
% - notes: 
%   - CONV is causal - need to offset by halfduration (just timing)
%   - study specfem3d - how do they work with kinematic sources
%   - 
% check nargin
error(nargchk(1,2,nargin))

% check data structure
error(seizmocheck(data,'dep'))

% toggle off struct checking
oldseizmocheckstate=get_seizmocheck_state;
set_seizmocheck_state(false);

% check headers
data=checkheader(data);

% check options


% toggle checking back
set_seizmocheck_state(oldseizmocheckstate);

end