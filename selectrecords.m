function [data,selected,h]=selectrecords(data,option,plottype,selected,varargin)
%SELECTRECORDS    Select or delete SEIZMO data records graphically
%
%    Description:  SELECTRECORDS(DATA,OPTION,PLOTTYPE) allows a user to select
%     records in DATA.  Selection/unselection is done by left clicking the 
%     record and exiting the function is done by middle-clicking the plot.
%     OPTION must be either 'keep' or 'delete' and determines what will be 
%     done with the selected records.  For option 'keep' the selected 
%     records will be given a green shaded background.  For option 'delete'
%     the background will be colored a shade of red.  PLOTTYPE controls the
%     plot type to use for record selection and must be (currently) 'p1' or
%     'p0'.  
%
%     SELECTRECORDS(DATA,OPTION,PLOTTYPE,VARARGIN) sends addition arguments
%     VARARGIN to the plot function.
%
%    Usage: [data]=selectrecords(data,option,plottype,varargin)
%
%    Examples:
%     To select which records to delete using p1:
%       data=selectrecords(data,'delete','p1')
%
%    See also: plot1, plot0

% check data structure
error(seizmocheck(data,'dep'))

% number of records
nrecs=numel(data);

% defaults
if(nargin<4); selected=false(nrecs,1); end
if(nargin<3); plottype='p0'; end
if(nargin<2); option='keep'; end

% default coloring
if(strcmpi(option,'delete'))
    color=[0.3 0 0];
elseif(strcmpi(option,'keep'))
    color=[0 0.3 0];
else
    error('seizmo:selectrecs:badInput','Unknown option: %s',option)
end

% plottype selection
button=0;
handles=ones(nrecs,1)*-1;
if(isequal(plottype,'p1'))
    % plot type 1
    [h,sh]=plot1(data,varargin{:});
    
    % color preselected
    bgcolors=cell2mat(get(sh,'color'));
    set(sh(selected),'color',color);
    
    while(button~=2)
        % get mouse button pressed
        [x,y,button]=ginput(1);
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
elseif(isequal(plottype,'p0'))
    % plot type 0
    [h]=plot0(data,varargin{:});
    handle=gca;
    
    % color preselected
    xlims=xlim(handle);
    clicked=find(selected); nclicked=numel(clicked);
    for i=1:nclicked
        handles(clicked(i))=patch([xlims(ones(1,2)) xlims(2*ones(1,2))],...
            [clicked(i)+0.5 clicked(i)-0.5 clicked(i)-0.5 clicked(i)+0.5],...
            color);
    end
    alpha(handles(selected),0.99) % alpha doesn't work right
    
    while(button~=2)
        % get mouse button pressed
        [x,y,button]=ginput(1);
        if(button==1 && isequal(handle,gca))
            % figure out which record from y position
            clicked=round(y);
            
            % check range
            if(clicked<1 || clicked>nrecs); continue; end
            
            % remove from list if in list and remove patch
            if(selected(clicked))
                selected(clicked)=false;
                delete(handles(clicked));
            % otherwise add to list and add patch
            else
                selected(clicked)=true;
                xlims=xlim(handle);
                handles(clicked)=patch([xlims(1) xlims(1) xlims(2) xlims(2)],...
                    [clicked+0.5 clicked-0.5 clicked-0.5 clicked+0.5],color);
                alpha(handles(clicked),0.99) % alpha doesn't work right
            end
        end
    end
else
    error('seizmo:selectrecs:badInput','Unsupported plottype: %s',plottype)
end

% handle data
if(strcmpi(option,'keep'))
    data=data(selected);
else
    data(selected)=[];
end

end
