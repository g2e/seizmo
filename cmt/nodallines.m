function [neu1,neu2]=nodallines(strike,dip,rake,step)
%NODALLINES    Returns nodal lines for a focal mechanism
%
%    Usage:    [neu1,neu2]=nodallines(sdr)
%              [neu1,neu2]=nodallines(strike,dip,rake)
%              [...]=nodallines(...,step)
%
%    Description:
%     [NEU1,NEU2]=NODALLINES(SDR) returns points along the nodal lines for
%     a given strike, dip & rake (and the calculated auxiliary plane).  The
%     nodal lines are the directions from which the P-wave energy from the
%     focal mechanism is theoretically expected to be zero.  SDR must be a
%     Nx3 array formatted as [strike dip rake] where N allows for multiple
%     focal mechanisms to be processed simultaneously.  Strike is clockwise
%     from North, dip is positive downward from the horizontal and rake is
%     counter-clockwise in the fault plane from the strike direction.  Note
%     that the strike must be such that when you look along the direction
%     of the strike the fault dips to your right.  The outputs NEU1 & NEU2
%     are 73xNx3 arrays where there are 73 points along each line, N lines,
%     and the 3 pages give the North, East & Up coordinates.  See the
%     Examples section below for an example of how to plot this info.
%
%     [NEU1,NEU2]=NODALLINES(STRIKE,DIP,RAKE) allows strike, dip & rake to
%     be given separately.
%
%     [...]=NODALLINES(...,STEP) alters the angular step size between
%     points on a nodal line (from the center of the focal mechanism).  The
%     default is 5 degrees making 73 points for each line.
%
%    Notes:
%
%    Examples:
%     % Plot the nodal lines for the first 10
%     % moment tensors in the GlobalCMT catalog:
%     cmts=findcmt('n',10);
%     sdr=[cmts.strike1 cmts.dip1 cmts.rake1];
%     [neu1,neu2]=nodallines(sdr);
%     figure;
%     plot3(neu1(:,:,1),neu1(:,:,2),neu1(:,:,3));
%     plot3(neu2(:,:,1),neu2(:,:,2),neu2(:,:,3));
%     xlabel('n'); ylabel('e'); zlabel('u');
%
%    See also: STRIKEDIP2NORM, NORM2STRIKEDIP, SDR2SLIP, SDR2NULL, SDR2TPB,
%              TPB2SDR, NORMSLIP2SDR, NEU2VPA, VPA2NEU, MT2TPB, TPB2MT

%     Version History:
%        Mar. 18, 2013 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 18, 2013 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,4,nargin));

% one or both inputs
switch nargin
    case {1 2}
        if(size(strike,2)~=3 || ndims(strike)>2)
            error('seizmo:nodallines:badInput',...
                'SDR must be a Nx3 array as [STRIKE DIP RAKE] !');
        elseif(~isnumeric(strike) || ~isreal(strike))
            error('seizmo:sdr2tpb:badInput',...
                'SDR must be a real-valued Nx3 array!');
        end
        if(nargin==2)
            if(~isscalar(dip) || ~isnumeric(dip) || ~isreal(dip) ...
                    || dip<=0 || dip>360)
                error('seizmo:sdr2tpb:badInput',...
                    'STEP must be a positive scalar (default is 5)!');
            end
            step=dip;
        else
            step=5;
        end
        [strike,dip,rake]=deal(strike(:,1),strike(:,2),strike(:,3));
    case {3 4}
        if(~isnumeric(strike) || ~isreal(strike) ...
                || ~isnumeric(dip) || ~isreal(dip) ...
                || ~isnumeric(rake) || ~isreal(rake))
            error('seizmo:nodallines:badInput',...
                'STRIKE/DIP/RAKE must be real-valued arrays!');
        end
        if(nargin==4)
            if(~isscalar(step) || ~isnumeric(step) || ~isreal(step) ...
                    || step<=0 || step>360)
                error('seizmo:sdr2tpb:badInput',...
                    'STEP must be a positive scalar (default is 5)!');
            end
        else
            step=5;
        end
        [strike,dip,rake]=expandscalars(strike,dip,rake);
    otherwise
        error('seizmo:nodallines:badNumInputs',...
            'Incorrect number of inputs (only 1 or 3)!');
end

% force row vectors for strike, dip, rake
strike=strike(:)';
dip=dip(:)';
rake=rake(:)';

% sweep through rakes for nodal lines
rakes=(0:step:360)';

% how am i going to get all this info into 2 outputs?
% NSTEPxNx3  x2
n=numel(dip);
ns=numel(rakes);
[neu1,neu2]=deal(nan(ns,n,3));
[neu1(:,:,1),neu1(:,:,2),neu1(:,:,3)]=sdr2slip(strike(ones(ns,1),:),...
    dip(ones(ns,1),:),rakes(:,ones(1,n)));
[strike2,dip2]=auxplane(strike,dip,rake);
[neu2(:,:,1),neu2(:,:,2),neu2(:,:,3)]=sdr2slip(strike2(ones(ns,1),:),...
    dip2(ones(ns,1),:),rakes(:,ones(1,n)));

end
