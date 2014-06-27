function []=noise_coda(indir,outdir,varargin)
%NOISE_CODA    Process surface wave coda for interstation Green's functions
%
%    Usage:    noise_coda(indir,outdir)
%              noise_coda(indir,outdir,steps)
%              noise_coda(indir,outdir,steps,'opt1',val,...,'optN',val)
%
%    Description:
%
%    Notes:
%     - 
%
%    Header changes: Varies with steps chosen...
%
%    Examples:
%     % 
%
%    See also: NOISE_C3, NOISE_PROCESS, NOISE_STACK, STACK2STACK,
%              NOISE_STACK_ARBITRARY, NOISE_STACK_DELAZ, NOISE_SETUP,
%              NOISE_OVERVIEW

%     Version History:
%        May  30, 2014 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  30, 2014 at 11:15 GMT

% todo:
% - event directory setup?
% - input is a directory of event directories with displacement seismograms
% - output is a directory of event directories with cross correlations
% - options:
%   - tdnorm
%   - fdnorm

end

