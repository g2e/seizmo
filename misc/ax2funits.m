function [varargout]=ax2funits(varargin)
%AX2FUNITS    Converts from data units to normalized figure units
%
%    Usage:    [xf,yf]=ax2funits(x,y)
%              posf=ax2funits(pos)
%              [...]=ax2funits(ax,...)
%
%    Description:
%     [Xf,Yf]=AX2FUNITS(X,Y) converts X,Y coordinates to normalized figure
%     units using the current axes.  This is useful for ANNOTATION input.  
%
%     POSf=AX2FUNITS(POS) converts a 4-element position vector POS from
%     data units to normalized figure units using the current axes.  The
%     position vector has the form [Xo Yo Width Height].
%
%     [...]=AX2FUNITS(AX,...) converts the units using the specified axes
%     handle AX.
%
%    Notes:
%
%    Examples:
%     % Make a simple plot with a double arrow:
%     figure;
%     [x,y]=ax2funits([.1 .9],[.1 .9]);
%     annotation('doublearrow',x,y);
%
%    See also: ANNOTATION

%     Version History:
%        Oct. 15, 2012 - edited from FEX #10656, less memory usage,
%                        simplify docs & error msg, allow Nx4 POS vector,
%                        use TRUE_AXIS_POSITION to handle axes tight etc
%
%     Written by Scott Hirsch (shirsch at mathworks dot com)
%                Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 15, 2012 at 13:35 GMT

% todo:

% check number of inputs
error(nargchk(1,3,nargin))

% get axes handle
if(numel(varargin{1})==1 ...
        && ishandle(varargin{1}) ...
        && strcmp(get(varargin{1},'type'),'axes'))	
    ax=varargin{1};
    varargin=varargin(2:end);
else
    ax=gca;
end

% check input
if(numel(varargin)==1)
    % POS vector
    if(size(varargin{1},2)~=4 && ndims(varargin{1})==2)
        error('seizmo:ax2funits:badInput',...
            'POS must be Nx4 vector as [X Y Width Height]!');
    end
else
    % X/Y matrices
    if(isscalar(varargin{1}))
        varargin{1}=varargin{1}(ones(size(varargin{2})));
    end
    if(isscalar(varargin{2}))
        varargin{2}=varargin{2}(ones(size(varargin{1})));
    end
    if(~isequal(size(varargin{1}),size(varargin{2})))
        error('seizmo:ax2funits:badInput',...
            'X/Y must be scalar or equal sized!');
    end
end

	
% get axes limits
axun=get(ax,'units');
set(ax,'units','normalized');
axpos=true_axis_position(ax); % works with axes tight
%axpos=get(ax,'position');
axlim=axis(ax);
axwidth=diff(axlim(1:2));
axheight=diff(axlim(3:4));

% transform data
if(numel(varargin)==1)
    varargin{1}(1)=(varargin{1}(1)-axlim(1))/axwidth*axpos(3)+axpos(1);
    varargin{1}(2)=(varargin{1}(2)-axlim(3))/axheight*axpos(4)+axpos(2);
    varargin{1}(3)=varargin{1}(3)*axpos(3)/axwidth;
    varargin{1}(4)=varargin{1}(4)*axpos(4)/axheight;
    varargout{1}=varargin{1};
else
    varargout{1}=(varargin{1}-axlim(1))*axpos(3)/axwidth+axpos(1);
    varargout{2}=(varargin{2}-axlim(3))*axpos(4)/axheight+axpos(2);
end

% restore axes units
set(ax,'units',axun)

end
