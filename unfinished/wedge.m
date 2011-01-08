function [varargout]=wedge(varargin)
%

% todo
% - everything
% - allow no data (ie wedge plot adjust)

% parse parameters
% - check for axis
% - indices of p/v start
[ax,varargin,nargs]=axescheck(varargin{:});
pvstart=3; % anything else is an error
for i=1:nargs
    if(ischar(varargin{i}))
        pvstart=i;
        break;
    end
end

% complain about strange p/v start or too few inputs
if(nargs<2 || pvstart~=3)
    error('seizmo:wedge:badInput',...
        'Must have exactly 2 data arguments!');
end

% allow standalone linestyle specifier at p/v start
if(mod(nargs-pvstart+1,2))
    varargin=[varargin 'linestyle' varargin];
    nargs=nargs+1;
end

% act according to nextplot state
ax=newplot(ax);

% get nextplot & hold state
next=lower(get(cax,'NextPlot'));
held=ishold(cax);

% parse p/v info
%
% - what are wedge properties (w)
% - what are axis properties (a)
% - what are lineseries properties (l)
%
% a activepositionproperty
% w box
% w clipping (turn this OFF for all polar objects)
% w backgroundcolor
% l color
% a colororder
% a dataaspectratio
% w fontangle
% w fontname
% w fontsize
% w fontunits
% w fontweight
% w gridlinestyle
% w axislinewidth
% l linewidth
% l marker property series
% w minorgridlinestyle
% a outerposition
% w parent (already stripped)
% a position
% w tickdir
% w ticklength
% a units
% w visible
% w raxislocation/taxislocation
% w rcolor/tcolor
% w rlim
% w rminorgrid
% w rminortick
% w rtick
% w rticklabel
% w same for t series
% w tdir
% w placement
% w tunits (internally always degrees, 
[w,a,l]=parse_wedge_pv(ax,varargin{pvstart:end});

% draw grid if not held
if(~held)
    % first determine auto rlim & tlim
    if(isempty(w.rlim))
        w.rlim=[0 max(abs(varargin{2}(~isinf(varargin{2}))))];
    end
    if(isempty(w.tlim))
        w.tlim=[0 360];
    end
    
    % number of theta points based on tlim
    tpts=-540:1:720;
    tpts=tpts(tpts>w.tlim(1) & tpts<w.tlim(2));
    tpts=[w.tlim(1) tpts w.tlim(2)];
    
    % get the x/y axis limits
    switch lower(w.placement)
        case 'center'
            % extreme x/y limits
            
        case 'origin'
            
    end
    
    % auto radial & theta tick setup
    
    
    % background
    
    % "axis"
    
    % grid
    
    % ticks
    
    % labels
end

% plot data on top

% cleanup

end


function [w,a,l]=parse_wedge_pv(ax,varargin)
%PARSE_WEDGE_PV    Separates out wedge/axis/lineseries properties
%
% todo:
% - convert all radians to degrees
% - 

% default wedge properties
w.box='on';
w.clipping='off';
w.backgroundcolor='w';
w.fontangle=get(cax,'defaulttextfontangle');
w.fontname=get(cax,'defaulttextfontname');
w.fontsize=get(ax,'defaulttextfontsize');
w.fontunits=get(ax,'defaulttextfontunits');
w.fontweight=get(ax,'defaulttextfontweight');
w.textunits='data';
w.gridlinestyle=':';
w.axislinewidth=0.5;
w.minorgridlinestyle=':';
w.tickdir=[];
w.ticklen=[];
w.visible='on';
w.raxislocation='out'; % in/out
w.taxislocation='cw'; % cw/ccw
w.rcolor='k';
w.tcolor='k';
w.rlim=[];
w.rminorgrid='off';
w.rminortick='off';
w.rtick=[];
w.rticklabel=[];
w.tlim=[];
w.tminorgrid
w.tminortick
w.ttick=[];
w.tticklabel=[];
w.tdir='normal';
w.tunits='degrees'; % degrees/radians
w.placement='center'; % center/origin

% replace defaults with those from held wedge
if(ishold(ax) && strcmp(get(ax,'tag'),'wedge'))
    
end

% property names
wfields=fieldnames(w);
afields={'activepositionproperty' 'colororder' 'dataaspectratio' ...
    'outerposition' 'position' 'units'};

% demand properties are strings
if(~iscellstr(varargin(1:2:end)))
    error('seizmo:wedge:badInput',...
        'Properties must be strings!');
end

% loop over each pair and parse
w={}; a={}; l={};
for i=1:2:nargin-1
    if(any(strcmpi(wfields,varargin{i})))
        w=[w varargin(i:i+1)];
    elseif(any(strcmpi(afields,varargin{a})))
        a=[a varargin(i:i+1)];
    else
        l=[l varargin(i:i+1)];
    end
end

% check & assign polar values
for i=1:2:numel(w)
    switch lower(varargin{i})
        case 'tdir'
            if(~ischar(varargin{i+1}) ...
                    || ~any(strcmpi(varargin{i+1},{'normal' 'reverse'})))
                error('seizmo:wedge:badInput',...
                    'TDIR must be ''normal'' or ''reverse''!');
            end
        otherwise
    end
    
    % assign
    w.(lower(varargin{i}))=varargin{i+1};
end

end


