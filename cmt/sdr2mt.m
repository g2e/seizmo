function [varargout]=sdr2mt(varargin)
%SDR2MT    Convert strike-dip-rake to moment tensor in Aki & Richards form
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

% convert NEU strike-dip-rake to A&R mt
% d = dip
% f = strike
% l = rake
% Mxx = - Mo(sind cosl sin2f + sin2d sinl sin2f )
% Myy = Mo(sind cosl sin2f - sin2d sinl cos2f )
% Mzz = Mo sin2d sinl
% Mxy = Mo(sind cosl cos2f + 0.5 sin2d sinl sin2f )
% Mxz = - Mo(cosd cosl cosf + cos2d sinl sinf )
% Myz = - Mo(cosd cosl sinf - cos2d sinl cosf )

% conversion
R2D=180/pi;

if(nargin==1)
    % check input
    sz=size(varargin{1});
    if(~isreal(varargin{1}) || ~isequal(3,sz(2)))
        error('seizmo:sdr2mt:badInput',...
            'SDR must be a real-valued Nx3 array!');
    end
    sz(2)=6;
    varargin{1}=varargin{1}/R2D;
    varargout{1}=nan(sz);
    varargout{1}(:,1,:)=sin(varargin{1}(:,2,:)).*cos(varargin{1}(:,3,:))...
        .*sin(2*varargin{1}(:,1,:))+sin(2*varargin{1}(:,2,:))...
        .*sin(varargin{1}(:,3,:)).*sin(varargin{1}(:,1,:)).^2;
    varargout{1}(:,2,:)=sin(varargin{1}(:,2,:)).*cos(varargin{1}(:,3,:))...
        .*sin(2*varargin{1}(:,1,:))-sin(2*varargin{1}(:,2,:))...
        .*sin(varargin{1}(:,3,:)).*cos(varargin{1}(:,1,:)).^2;
    varargout{1}(:,3,:)=sin(2*varargin{1}(:,2,:)).*sin(varargin{1}(:,3,:));
    varargout{1}(:,4,:)=sin(varargin{1}(:,2,:)).*cos(varargin{1}(:,3,:))...
        .*cos(2*varargin{1}(:,1,:))+0.5*sin(2*varargin{1}(:,2,:))...
        .*sin(varargin{1}(:,3,:)).*sin(2*varargin{1}(:,1,:));
    varargout{1}(:,5,:)=cos(varargin{1}(:,2,:)).*cos(varargin{1}(:,3,:))...
        .*cos(varargin{1}(:,1,:))+cos(2*varargin{1}(:,2,:))...
        .*sin(varargin{1}(:,3,:)).*sin(varargin{1}(:,1,:));
    varargout{1}(:,6,:)=cos(varargin{1}(:,2,:)).*cos(varargin{1}(:,3,:))...
        .*sin(varargin{1}(:,1,:))-cos(2*varargin{1}(:,2,:))...
        .*sin(varargin{1}(:,3,:)).*cos(varargin{1}(:,1,:));
    varargout{1}(:,[1 5 6],:)=-varargout{1}(:,[1 5 6],:);
elseif(nargin==3)
    % check inputs
    if(any(~cellfun('isreal',varargin) | cellfun('size',varargin,2)~=1))
        error('seizmo:sdr2mt:badInput',...
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
                error('seizmo:sdr2mt:badInput',...
                    'Non-scalar inputs must be equal sized!');
            end
        end
    end
    % combine
    varargin{1}=cat(2,varargin{:});
    % convert
    sz(2)=6;
    varargin{1}=varargin{1}/R2D;
    varargout{1}=nan(sz);
    varargout{1}(:,1,:)=sin(varargin{1}(:,2,:)).*cos(varargin{1}(:,3,:))...
        .*sin(2*varargin{1}(:,1,:))+sin(2*varargin{1}(:,2,:))...
        .*sin(varargin{1}(:,3,:)).*sin(varargin{1}(:,1,:)).^2;
    varargout{1}(:,2,:)=sin(varargin{1}(:,2,:)).*cos(varargin{1}(:,3,:))...
        .*sin(2*varargin{1}(:,1,:))-sin(2*varargin{1}(:,2,:))...
        .*sin(varargin{1}(:,3,:)).*cos(varargin{1}(:,1,:)).^2;
    varargout{1}(:,3,:)=sin(2*varargin{1}(:,2,:)).*sin(varargin{1}(:,3,:));
    varargout{1}(:,4,:)=sin(varargin{1}(:,2,:)).*cos(varargin{1}(:,3,:))...
        .*cos(2*varargin{1}(:,1,:))+0.5*sin(2*varargin{1}(:,2,:))...
        .*sin(varargin{1}(:,3,:)).*sin(2*varargin{1}(:,1,:));
    varargout{1}(:,5,:)=cos(varargin{1}(:,2,:)).*cos(varargin{1}(:,3,:))...
        .*cos(varargin{1}(:,1,:))+cos(2*varargin{1}(:,2,:))...
        .*sin(varargin{1}(:,3,:)).*sin(varargin{1}(:,1,:));
    varargout{1}(:,6,:)=cos(varargin{1}(:,2,:)).*cos(varargin{1}(:,3,:))...
        .*sin(varargin{1}(:,1,:))-cos(2*varargin{1}(:,2,:))...
        .*sin(varargin{1}(:,3,:)).*cos(varargin{1}(:,1,:));
    varargout{1}(:,[1 5 6],:)=-varargout{1}(:,[1 5 6],:);
else
    error('seizmo:sdr2mt:badNumInputs',...
        'Incorrect number of inputs!');
end

% split up if wanted
if(nargout<=1)
    % Nx6 form (already done)
elseif(nargout<=6)
    % 6 Nx1 form
    varargout=num2cell(varargout{1},1);
else
    error('seizmo:sdr2mt:badNumOutputs',...
        'Incorrect number of outputs!');
end

end
