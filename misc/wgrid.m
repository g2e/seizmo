function []=wgrid(varargin)
%WGRID    Wedge grid lines
%
%    Usage:    wgrid on
%              wgrid off
%              wgrid
%              wgrid(ax,...)
%
%    Description:
%     WGRID ON adds grid lines to the current axes if drawn by WEDGE.
%     
%     WGRID OFF removes the grid lines from the current axes if it
%     was drawn with WEDGE.
%
%     WGRID toggles the grid lines of the current axes if it was
%     drawn with WEDGE.
%
%     WGRID(AX,...) works with axes AX rather than the current axes.
%
%    Notes:
%     - wedge(ax,'rgrid','on','azgrid','on');
%
%    Examples:
%     % Plot and clear the grid:
%     wedge(30:30:300,1:10);
%     wgrid off;
%
%    See also: WEDGE, AZLABEL, RLABEL

%     Version History:
%        May   1, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May   1, 2012 at 18:35 GMT

% todo:

% check number of inputs
error(nargchk(0,2,nargin));

% extract handle
if(numel(varargin) && isnumeric(varargin{1}))
    if(all(ishghandle(varargin{1},'axes')))
        ax=varargin{1};
        varargin(1)=[];
    else
        error('seizmo:wgrid:badHandle',...
            'Given handle is not of type Axes!');
    end
else
    ax=gca;
end

% check handle is wedge
if(~all(strcmp(get(ax,'createfcn'),'wedge')))
    error('seizmo:wgrid:badHandle',...
        'Axes handle was not created by function WEDGE!');
end

% loop over each axis
for a=1:numel(ax)
    % get wedge axes current parameters
    w=getappdata(ax(a));

    % get all appropriate handles
    rh=findall(ax(a),'tag','wedge_grid_r');
    ah=findall(ax(a),'tag','wedge_grid_az');
    
    % state given?
    if(numel(varargin))
        % check state
        if(~ischar(varargin{1}) || ~any(strcmpi(varargin{1},{'on' 'off'})))
            error('seizmo:wgrid:badState',...
                'Grid state must be ''ON'' or ''OFF''!');
        end
        
        % wgrid off turns off the minorgrid
        if(strcmpi(varargin{1},'off'))
            rh=[rh; findall(ax(a),'tag','wedge_minorgrid_r')];
            ah=[ah; findall(ax(a),'tag','wedge_minorgrid_az')];
            setappdata(ax(a),'rminorgrid','off');
            setappdata(ax(a),'azminorgrid','off');
        end
        
        % force state of all handles
        set([rh;ah],'visible',onoffstate(w.axis,varargin{1}));
        
        % set parameter field
        setappdata(ax(a),'rgrid',varargin{1})
        setappdata(ax(a),'azgrid',varargin{1});
    else % toggle state
        state={'on' 'off'}';
        if(numel([rh;ah])==1)
            set([rh;ah],'visible',onoffstate(w.axis,...
                state{strcmp(get([rh;ah],'visible'),'on')+1}));
        else
            set([rh;ah],{'visible'},onoffstate(w.axis,...
                state(strcmp(get([rh;ah],'visible'),'on')+1)));
        end
        setappdata(ax(a),'rgrid',...
            state{strcmp(getappdata(ax(a),'rgrid'),'on')+1});
        setappdata(ax(a),'azgrid',...
            state{strcmp(getappdata(ax(a),'azgrid'),'on')+1});
    end
end

end


function [s]=onoffstate(varargin)
% childonoff=onoffstate(parentonoff,childonoff)
for i=1:nargin
    if(strcmpi(varargin{i},'off')); s=varargin{i}; return; end
    s=varargin{end};
end
end

