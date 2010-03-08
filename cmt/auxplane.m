function [strike,dip,rake]=auxplane(varargin)
%AUXPLANE    Returns strike-dip-slip of 2nd (auxiliary) focal plane
%
%    Usage:
%
%    Description:
%
%    Notes:
%
%    Examples:
%
%    See also:

%     Version History:
%        Mar.  8, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  8, 2010 at 13:50 GMT

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
