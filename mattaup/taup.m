function []=taup()
%TAUP    Invokes the TauP java applet
%
%    Usage:    taup
%
%    Description:
%     TAUP calls the TauP toolkit's java applet, which allows the user to
%     change models, event depth, station distance, and phases
%     interactively to calculate travel times, pierce points and ray paths.
%     Unfortunately none of the calculated data is returned.  See the
%     functions listed below to get usable output from TauP.
%
%    Notes:
%     - TauP toolkit developed by:
%        H. Philip Crotwell, Thomas J. Owens, Jeroen Ritsema
%        Department of Geological Sciences
%        University of South Carolina
%        http://www.seis.sc.edu
%        crotwell@seis.sc.edu
%
%    Examples:
%     % Only possible usage example:
%     taup
%
%    See also: TAUPCURVE, TAUPPATH, TAUPPIERCE, TAUPTIME, TAUPCREATE

%     Version History:
%        Sep.  5, 2009 - major doc update
%        Feb. 24, 2012 - doc update, make work with Octave
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 24, 2012 at 21:40 GMT

% todo:

% try adding *.jar if no TauP class exists
if(~exist('edu.sc.seis.TauP.TauP','class'))
    fs=filesep;
    mypath=fileparts(mfilename('fullpath'));
    javaaddpath([mypath fs 'lib' fs 'MatTauP-1.2beta4.jar']);
    javaaddpath([mypath fs 'lib' fs 'TauP-1.2beta4.jar']);
    javaaddpath([mypath fs 'lib' fs 'seisFile-1.0.8.jar']);
end

t=javaObject('edu.sc.seis.TauP.TauP');
t.show;

end
