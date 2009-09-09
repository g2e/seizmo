function []=taup()
%TAUP    Invokes the TauP java applet
%
%    Usage:    taup
%
%    Description: TAUP calls the TauP toolkit's java applet, which allows
%     the user to change models, event depth, station distance, and phases
%     interactively to calculate travel times, pierce points and ray paths.
%     Unfortunately none of the calculated data is returned.  See the
%     functions listed below to get usable output from TauP.
%
%    Notes:
%     - These scripts require the included file:
%          mattaup/lib/matTaup.jar
%       to be added in Matlab's javaclasspath.  You may use the functions
%       javaaddpath and javarmpath to alter the dynamic portion of the path
%       or you will need to add the jar file to the Matlab system file
%       'classpath.txt'.  Use the command 'edit classpath.txt' in Matlab to
%       add the the jar file (be careful and use the full path!).  This may
%       require administrator privileges.
%
%     - MatTauP is only a wrapping program for TauP toolkit, which is
%       developed by:
%        H. Philip Crotwell, Thomas J. Owens, Jeroen Ritsema
%        Department of Geological Sciences
%        University of South Carolina
%        http://www.seis.sc.edu
%        crotwell@seis.sc.edu
%
%     - MatTauP was written by:
%        Qin Li 
%        Unverisity of Washington
%        qinli@u.washington.edu
%        Nov, 2002
%
%    Examples:
%     Only possible usage example:
%      taup
%
%    See also: taupcurve, tauppath, tauppierce, tauptime

%     Version History:
%        Sep.  5, 2009 - major doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep.  5, 2009 at 21:40 GMT

% todo:

import edu.sc.seis.TauP.*
t=TauP;
t.show;

end
