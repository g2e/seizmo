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
%        Jan. 26, 2014 - no longer need to update jar filenames
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 26, 2014 at 21:40 GMT

% todo:

% try adding *.jar if no TauP class exists
if(~exist('edu.sc.seis.TauP.TauP','class'))
    fs=filesep;
    mypath=fileparts(mfilename('fullpath'));
    jars=dir([mypath fs 'lib' fs '*.jar']);
    for i=1:numel(jars)
        if(~ismember([mypath fs 'lib' fs jars(i).name],javaclasspath))
            javaaddpath([mypath fs 'lib' fs jars(i).name]);
        end
    end
end

t=javaObject('edu.sc.seis.TauP.TauP');
t.show;

end
