function [varargout]=clonefigure(in,out)
%CLONEFIGURE    Makes a clone of a figure
%
%    Usage:    clonefigure
%              clonefigure(in)
%              clonefigure(in,out)
%              out=clonefigure(...)
%
%    Description:
%     CLONEFIGURE clones the current figure to a new figure.
%
%     CLONEFIGURE(IN) clones the figures specified by the handles in IN
%     to new figures.
%
%     CLONEFIGURE(IN,OUT) clones the figures given by handles in IN
%     to the figures given by handles in OUT.  If a figure corresponding to
%     one of the handles in OUT does not exist it is created.
%
%     OUT=CLONEFIGURE(...) returns the cloned figure handles.
%
%    Notes:
%     - Will not clone a figure to itself successfully.
%     - Will not switch figures successfully ie. clonefigure([1 2],[2 1])
%
%    Examples:
%     % clone figure 1 to figure 2 & 3
%     clonefigure(1,2:3)
%
%     % clone figures 1,2,3
%     clonefigure(1:3)
%
%    See also: COPYOBJ, FIGURE

%     Version History:
%        Aug.  5, 2010 - initial version
%        Oct. 11, 2010 - skip 'FixedColors', better colorbar cloning, fix
%                        crash when adjusting text visibility (may cause
%                        regression!)
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 11, 2010 at 12:25 GMT

% todo:

% check nargin
error(nargchk(0,2,nargin));

% defaults
if(nargin<1 || isempty(in)); in=gcf; end
if(nargin<2 || isempty(out)); out=nan(size(in)); end

% check figure handles
if(~isreal(in) || any(~ishandle(in(:))) ...
        || any(~strcmp('figure',get(in(:),'type'))))
    error('seizmo:clonefigure:badInput',...
        'IN must be a valid figure handle!');
elseif(nargin>1 && (~isreal(out) || any(out<1) || any(out~=fix(out))))
    error('seizmo:clonefigure:badInput',...
        'OUT must be a valid figure handle!');
elseif(~isequalsizeorscalar(in,out))
    error('seizmo:clonefigure:badInput',...
        'Inputs must be equal sized or scalar!');
elseif(isscalar(out) && ~isscalar(in))
    error('seizmo:clonefigure:badInput',...
        'Can not clone multiple figures to a single figure!');
end

% expand scalars
[in,out]=expandscalars(in,out);

% loop over each output figure
hin=cell(numel(out),1); hout=hin;
for i=1:numel(out)
    % bring to front
    if(isnan(out(i)))
        out(i)=figure;
    else
        figure(out(i));
    end
    
    % clear output figure
    clf(out(i));
    
    % recursive object copy
    [hin{i},hout{i}]=rcopyobj(in(i),out(i));
    
    % copy all figure props minus a few fields
    set(out(i),rmfield(get(in(i)),{'BeingDeleted' 'Children' ...
        'CurrentAxes' 'CurrentCharacter' 'CurrentObject' 'CurrentPoint' ...
        'Type' 'FixedColors'}));
    
    % MATLAB BUG WORKAROUNDS
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % need to force drawing to get actual values
    drawnow;
    % colorbar position & axes fixes
    cbin=findobj(in(i),'tag','Colorbar');
    cbout=findobj(out(i),'tag','Colorbar');
    for j=1:numel(cbin)
        set(cbout(j),'position',get(cbin(j),'position'),...
            'axes',hout{i}(hin{i}==get(cbin(j),'axes')));
    end
    % text visibility fixes
    tin=findall(in(i),'type','text');
    tout=findall(out(i),'type','text');
    if(numel(tin)==numel(tout))
        for j=1:numel(tin)
            set(tout(j),'visible',get(tin(j),'visible'));
        end
    else
        % not sure yet
    end
end

% output figure handle if wanted
if(nargout); varargout{1}=out; end

end

function [kids,newkids]=rcopyobj(oldparent,newparent)
%RCOPYOBJ    Recursive object copy
kids=get(oldparent,'children');
if(~isempty(kids))
    newkids=copyobj(kids,newparent);
    nkids=numel(newkids);
    k=cell(1,nkids); nk=k;
    for i=1:numel(newkids)
        [k{i},nk{i}]=rcopyobj(kids(i),newkids(i));
        k{i}=k{i}(:)';
        nk{i}=nk{i}(:)';
    end
    kids=[kids(:); cell2mat(k)'];
    newkids=[newkids(:); cell2mat(nk)'];
else
    newkids=kids;
end
end
