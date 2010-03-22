function [strike,dip,rake]=auxplane(varargin)
%AUXPLANE    Returns strike-dip-slip of 2nd (auxiliary) focal plane
%
%    Usage:  [strike,dip,rake]=auxplane(strike,dip,rake)
%            [strike,dip,rake]=auxplane(sdr)
%            sdr=auxplane(...)
%
%    Description: [STRIKE,DIP,RAKE]=AUXPLANE(STRIKE,DIP,RAKE) finds the
%     corresponding auxilary plane(s) and rake(s) that are able to match
%     the strain(s) produced by the input plane(s) and rake(s).  STRIKE,
%     DIP, & RAKE should be column vectors of equal length or scalars
%     (scalar expansion is performed).  The auxilary plane is perpendicular
%     to the input fault plane.  The plane that is perpendicular to both
%     the fault plane and the auxilary plane has a normal along the null
%     axis of the focal mechanism.
%
%     [STRIKE,DIP,RAKE]=AUXPLANE(SDR) does the above with a combined array
%     SDR that is Nx3 where the array is equal to [STRIKE DIP RAKE].
%
%     SDR=AUXPLANE(...) outputs a combined array of [STRIKE DIP RAKE] for
%     the auxilary plane.
%
%    Notes:
%     - The auxilary plane is perpendicular to the input fault plane.
%     - The plane perpendicular to both the fault plane and the auxilary
%       plane has its normal along the null axis of the focal mechanism.
%
%    Examples:
%     Corresponding fault plane and rake for dip-slip on a normal fault
%     that strikes North-South:
%      [strike,dip,rake]=auxplane(0,45,-90)
%
%    See also: STRIKEDIP, NS2SDR, SDR2NS, MT2SDR, SDR2MT

%     Version History:
%        Mar.  8, 2010 - initial version
%        Mar. 22, 2010 - added docs
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 22, 2010 at 15:05 GMT

% todo:

% for radians/degrees
R2D=180/pi;

if(nargin==1)
    % check input
    sz=size(varargin{1});
    if(~isreal(varargin{1}) || ~isequal(3,sz(2)))
        error('seizmo:auxplane:badInput',...
            'SDR must be a real-valued Nx3 array!');
    end
    
    % convert to radians
    varargin{1}=varargin{1}/R2D;
    
    % make strike relative to west (why?)
    varargin{1}(:,1)=varargin{1}(:,1)+pi/2;
    
    % cos & sin
    cos1=cos(varargin{1});
    sin1=sin(varargin{1});
    
    % vector pointing in direction of slip
    % is normal to the auxiliary plane
    sl1=-cos1(:,3).*cos1(:,1)-sin1(:,3).*sin1(:,1).*cos1(:,2);
    sl2=cos1(:,3).*sin1(:,1)-sin1(:,3).*cos1(:,1).*cos1(:,2);
    sl3=sin1(:,3).*sin1(:,2);
    
    % strike & dip of auxiliary plane
    [strike,dip]=strikedip(sl2,sl1,sl3);
    
    % rake direction is normal to the primary focal plane
    % - getting normal to primary focal plane
    n1=sin1(:,1).*sin1(:,2);
    n2=cos1(:,1).*sin1(:,2);
    %n3=cos1(:,2);
    
    % a vector along intersection of aux. plane with horizontal
    h1=-sl2;
    h2=sl1;
    %h3=0*sl1;
    
    % angle between normal vector and horizontal vector gives rake (+/-)
    rake=R2D*acos((h1.*n1+h2.*n2)./sqrt(h1.*h1+h2.*h2));
    
    % fix sense of orientation
    j=find(sl3<=0);
    rake(j)=-rake(j);
elseif(nargin==3)
    % check inputs
    if(any(~cellfun('isreal',varargin) | cellfun('size',varargin,2)~=1))
        error('seizmo:auxplane:badInput',...
            'All inputs must be real-valued Nx1 vectors!');
    end
    % expand scalars
    n=cellfun('prodofsize',varargin);
    sz=size(varargin{find(n==max(n),1,'first')});
    for i=1:3
        if(n(i)==1)
            varargin{i}=varargin{i}(ones(sz));
        else
            if(~isequal(size(varargin{i}),sz))
                error('seizmo:auxplane:badInput',...
                    'Non-scalar inputs must be equal sized!');
            end
        end
    end
    % combine
    varargin{1}=cat(2,varargin{:});
    
    % convert to radians
    varargin{1}=varargin{1}/R2D;
    
    % make strike relative to west (why?)
    varargin{1}(:,1)=varargin{1}(:,1)+pi/2;
    
    % cos & sin
    cos1=cos(varargin{1});
    sin1=sin(varargin{1});
    
    % vector pointing in direction of slip
    % is normal to the auxiliary plane
    sl1=-cos1(:,3).*cos1(:,1)-sin1(:,3).*sin1(:,1).*cos1(:,2);
    sl2=cos1(:,3).*sin1(:,1)-sin1(:,3).*cos1(:,1).*cos1(:,2);
    sl3=sin1(:,3).*sin1(:,2);
    
    % strike & dip of auxiliary plane
    [strike,dip]=strikedip(sl2,sl1,sl3);
    
    % rake direction is normal to the primary focal plane
    % - getting normal to primary focal plane
    n1=sin1(:,1).*sin1(:,2);
    n2=cos1(:,1).*sin1(:,2);
    %n3=cos1(:,2);
    
    % a vector along intersection of aux. plane with horizontal
    h1=-sl2;
    h2=sl1;
    %h3=0*sl1;
    
    % angle between normal vector and horizontal vector gives rake (+/-)
    rake=R2D*acos((h1.*n1+h2.*n2)./sqrt(h1.*h1+h2.*h2));
    
    % fix sense of orientation
    j=find(sl3<=0);
    rake(j)=-rake(j);
else
    error('seizmo:auxplane:badNumInputs',...
        'Incorrect number of inputs!');
end

% split up if wanted
if(nargout<=1)
    % Nx3 form
    strike=cat(2,strike,dip,rake);
elseif(nargout<=3)
    % 3 Nx1 form already done
else
    error('seizmo:auxplane:badNumOutputs',...
        'Incorrect number of outputs!');
end

end
