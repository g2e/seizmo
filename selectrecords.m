function [data,selected,h]=selectrecords(data,varargin)
%SELECTRECORDS    Select or delete SEIZMO data records graphically
%
%    Usage:    [data,selected,h]=selectrecords(data)
%              data=selectrecords(data,option)
%              data=selectrecords(data,option,type)
%              data=selectrecords(data,option,type,selected)
%              data=selectrecords(data,option,type,selected,plotoptions)
%
%    Description: [DATA,SELECTED,H]=SELECTRECORDS(DATA) returns the records
%     in SEIZMO data structure DATA that are graphically selected by the
%     user.  By default the plottype is PLOT0 and no records are
%     preselected.  Selection/unselection of records is performed by left-
%     clicking over a record.  Complete dataset selection by middle-
%     clicking over the plot or closing the figure.  Optional additional
%     outputs are the logical array SELECTED which indicates the records
%     that were selected and the plot handle H.
%
%     SELECTRECORDS(DATA,OPTION) sets whether selected records from DATA
%     are kept or deleted.  OPTION must be either 'keep' or 'delete'.  When
%     OPTION is 'keep', the background color for selected records is set to
%     a dark green.  For OPTION set to 'delete', the background color is
%     set to a dark red for selected records.  The default is 'keep'.
%
%     SELECTRECORDS(DATA,OPTION,TYPE) sets the plot type to be used in
%     record selection.  TYPE must be one of 'plot0','plot1','p0', or 'p1'.
%     The default is 'plot0'.
%
%     SELECTRECORDS(DATA,OPTION,TYPE,SELECTED) allows preselecting records
%     in DATA using the logical array SELECTED.  SELECTED must be either
%     true (all selected), false (all unselected), or a logical array with
%     the same number of elements as DATA.  The default is false.
%
%     SELECTRECORDS(DATA,OPTION,TYPE,SELECTED,PLOTOPTIONS) passes plotting
%     options PLOTOPTIONS (all arguments after SELECTED) to the plotting
%     function chosen with TYPE.
%
%    Notes:
%
%    Header changes: NONE
%
%    Examples:
%     To select which records to delete using plot1:
%       data=selectrecords(data,'delete','plot1')
%
%    See also: plot1, plot0

%     Version History:
%        Apr. 17, 2008 - initial version (happy bday bro!)
%        Nov. 16, 2008 - update for name changes
%        Nov. 30, 2008 - update for name changes
%        Dec. 13, 2008 - minor bug fixes
%        Mar. 13, 2009 - doc fixes
%        Apr.  3, 2009 - added preselect argument, combine keep/delete code
%                        redundancy
%        Apr. 23, 2009 - fix seizmocheck for octave, move usage up
%        May  30, 2009 - major doc update, major code cleaning
%        June  2, 2009 - fixed history, patches go on bottom now
%
%     Testing Table:
%                                  Linux    Windows     Mac
%        Matlab 7       r14        
%               7.0.1   r14sp1
%               7.0.4   r14sp2
%               7.1     r14sp3
%               7.2     r2006a
%               7.3     r2006b
%               7.4     r2007a
%               7.5     r2007b
%               7.6     r2008a
%               7.7     r2008b
%               7.8     r2009a
%        Octave 3.2.0
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated June  2, 2009 at 19:00 GMT

% todo:

% check nargin
if(nargin<1)
    error('seizmo:selectrecords:notEnoughInputs',...
        'Not enough input arguments.');
elseif(nargin>4 && mod(nargin,2))
    error('seizmo:selectrecords:plotOptionMustBePaired',...
        'Plot options must be paired with a value!');
end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% number of records
nrecs=numel(data);

% valid values for options
valid.OPERATION={'keep' 'delete'};
valid.TYPE={'plot0' 'plot1' 'p0' 'p1'};

% defaults
option.OPERATION='keep';
option.TYPE='plot1';
option.SELECTED=false;
option.DELETECOLOR=[0.3 0 0];
option.KEEPCOLOR=[0 0.3 0];

% get options from SEIZMO global
global SEIZMO
try
    fields=fieldnames(SEIZMO.SELECTRECORDS);
    for i=1:numel(fields)
        if(~isempty(SEIZMO.SELECTRECORDS.(fields{i})))
            option.(fields{i})=SEIZMO.SELECTRECORDS.(fields{i});
        end
    end
catch
end

% get options from command line
for i=1:min([3 numel(varargin)])
    if(~isempty(varargin{i}))
        switch i
            case 1
                option.OPERATION=varargin{i};
            case 2
                option.TYPE=varargin{i};
            case 3
                option.SELECTED=varargin{i};
        end
    end
end
varargin(1:min([3 numel(varargin)]))=[];

% check options
fields=fieldnames(option);
for i=1:numel(fields)
    % specific checks
    switch lower(fields{i})
        case {'operation' 'type'}
            if(~ischar(option.(fields{i})) || ...
                    size(option.(fields{i}),1)~=1 || ...
                    ~any(strcmpi(option.(fields{i}),valid.(fields{i}))))
                error('seizmo:selectrecords:badInput',...
                    ['%s option must be one of the following:\n'...
                    sprintf('%s ',valid.(fields{i}){:})],fields{i});
            end
        case 'selected'
            if(isnumeric(option.SELECTED))
                option.SELECTED=logical(option.SELECTED);
            end
            if(~islogical(option.SELECTED) || ...
                    ~any(numel(option.SELECTED)==[1 nrecs]))
                error('seizmo:selectrecords:badInput',...
                    ['SELECTED option must be a logical scalar or a\n' ...
                    'logical array with one element per record!']);
            end
            if(isscalar(option.SELECTED))
                option.SELECTED=option.SELECTED(ones(nrecs,1),1);
            end
        case {'keepcolor' 'deletecolor'}
            % let plot routines handle this
    end
end

% set color
if(strcmpi(option.OPERATION,'keep'))
    color=option.KEEPCOLOR;
else
    color=option.DELETECOLOR;
end

% proceed by plot type
button=0; handles=ones(nrecs,1)*-1;
selected=option.SELECTED;
switch option.TYPE
    case {'plot0' 'p0'}
        h=plot0(data,varargin{:});
        handle=gca;
        
        % patch under all
        xlims=xlim(handle);
        for i=1:nrecs
            handles(i)=patch([xlims(ones(1,2)) xlims(2*ones(1,2))],...
                [i+0.5 i-0.5 i-0.5 i+0.5],color);
            uistack(handles(i),'bottom');
        end
        set(handles(~selected),'visible','off')
        
        while(button~=2)
            % get mouse button pressed
            try
                [x,y,button]=ginput(1);
            catch
                h=[];
                button=2;
            end
            if(button==1 && isequal(handle,gca))
                % figure out which record from y position
                clicked=round(y);
                
                % check range
                if(clicked<1 || clicked>nrecs); continue; end
                
                % remove from selected if already selected, turn off patch
                if(selected(clicked))
                    selected(clicked)=false;
                    set(handles(clicked),'visible','off');
                % otherwise add to selected, turn on patch
                else
                    selected(clicked)=true;
                    set(handles(clicked),'visible','on');
                end
            end
        end
    case {'plot1' 'p1'}
        % plot type 1
        [h,sh]=plot1(data,varargin{:});
        
        % color preselected
        bgcolors=cell2mat(get(sh,'color'));
        set(sh(selected),'color',color);
        
        while(button~=2)
            % get mouse button pressed
            try
                [x,y,button]=ginput(1);
            catch
                h=[];
                button=2;
            end
            if(button==1)
                % grab axis handle
                handle=gca;
                
                % figure out which record
                clicked=find(handle==sh,1);
                
                % remove from list if in list and change color
                if(selected(clicked))
                    selected(clicked)=false;
                    set(handle,'color',bgcolors(clicked,:));
                % otherwise add to list and change color
                else
                    selected(clicked)=true;
                    set(handle,'color',color);
                end
            end
        end
end

% handle data
if(strcmpi(option.OPERATION,'keep'))
    data=data(selected);
else
    data(selected)=[];
end

end
