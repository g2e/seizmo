function [ax,vararg]=axparse(varargin)
%AXPARSE    Strips out axes arguments (leading or p/v pair)
%
%    Usage:    [ax,vararg]=axparse(varargin{:})
%
%    Description:
%     [AX,VARARGIN]=AXPARSE(VARARGIN{:}) processes leading axes inputs and
%     trailing p/v pair axes inputs, returning the axes and any remaining
%     inputs.  Valid property strings are 'parent', 'axes' & 'ax'.  This is
%     a helper function for other functions with complex argument lists.
%
%    Notes:
%
%    Examples:
%     % A silly example:
%     ax=axes;
%     [ax,varargin]=axparse(ax,2:20,3:21,'parent',ax);
%     plot(ax,varargin{:})
%
%    See also: AXESCHECK, __PLT_GET_AXIS_ARG__

%     Version History:
%        Feb. 24, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 24, 2012 at 17:15 GMT

% todo:

% defaults
ax=[]; vararg=varargin;

% strip leading single axis
if(numel(vararg) && isscalar(vararg{1}) && ishghandle(vararg{1},'axes'))
    ax=vararg{1};
    vararg=vararg(2:end);
end

% strip pertinent prop/val pairs
if(numel(vararg))
    axprop=find(strcmpi('parent',vararg) | strcmpi(vararg,'ax') ...
        | strcmpi(vararg,'axes'));
    if(numel(axprop))
        axprop=unique([axprop axprop+1]);
        if(axprop(end)<=numel(vararg))
            if(all(ishghandle(vararg{axprop(end)},'axes')))
                ax=vararg{axprop(end)};
                vararg(axprop)=[];
            end
        else
            % missing the axes value?
        end
    end
end

end
