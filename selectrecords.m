function [data]=selectrecords(data,option,plottype,varargin)
%SELECTRECS    Select or delete SEIZMO data records graphically
%
%    Description:  SELECTRECS(DATA,OPTION,PLOTTYPE) allows a user to select
%     records in DATA.  Selection/unselection is done by left clicking the 
%     record and exiting the function is done by middle-clicking the plot.
%     OPTION must be either 'keep' or 'delete' and determines what will be 
%     done with the selected records.  For option 'keep' the selected 
%     records will be given a green shaded background.  For option 'delete'
%     the background will be colored a shade of red.  PLOTTYPE controls the
%     plot type to use for record selection and must be (currently) 'p1' or
%     'p3'.  
%
%     SELECTRECS(DATA,OPTION,PLOTTYPE,VARARGIN) sends addition arguments
%     VARARGIN to the plot function.
%
%    Usage: [data]=selectrecs(data,option,plottype,varargin)
%
%    Examples:
%     To select which records to delete using p1:
%       data=selectrecs(data,'delete','p1')
%
%    See also: p1, p3

% check data structure
error(seizmocheck(data,'dep'))

% defaults
if(nargin<3); plottype='plot0'; end
if(nargin<2); option='keep'; end

% plottype selection
nrecs=numel(data);
button=0;
selected=false(nrecs,1);
bgcolors=zeros(nrecs,3);
handles=ones(nrecs,1)*-1;
if(isequal(plottype,'plot1'))
    % plot type 1
    [h,sh]=plot1(data,varargin{:});
    
    % are we deleting or keeping
    if(isequal(option,'keep'))
        while(button~=2)
            % bring plot to focus
            figure(h);
            
            % get mouse button pressed
            [x,y,button]=ginput(1);
            if(button==1)
                % grab axis handle and background color of clicked subplot
                handle=gca;
                bgcolor=get(handle,'color');
                
                % figure out which record
                clicked=find(handle==sh,1);
                
                % remove from list if in list and change color
                if(selected(clicked))
                    selected(clicked)=false;
                    set(handle,'color',bgcolors(clicked,:));
                % otherwise add to list and change color
                else
                    selected(clicked)=true;
                    bgcolors(clicked,:)=bgcolor;
                    set(handle,'color',[0 0.3 0]);
                end
            end
        end
        data=data(selected);
    elseif(isequal(option,'delete'))
        while(button~=2)
            % bring plot to focus
            figure(h);
            
            % get mouse button pressed
            [x,y,button]=ginput(1);
            if(button==1)
                % grab axis handle and background color of clicked subplot
                handle=gca;
                bgcolor=get(handle,'color');
                
                % figure out which record
                clicked=find(handle==sh,1);
                
                % remove from list if in list and change color
                if(selected(clicked))
                    selected(clicked)=false;
                    set(handle,'color',bgcolors(clicked,:));
                % otherwise add to list and change color
                else
                    selected(clicked)=true;
                    bgcolors(clicked,:)=bgcolor;
                    set(handle,'color',[0.3 0 0]);
                end
            end
        end
        data(selected)=[];
    else
        error('seizmo:selectrecs:badInput','Unknown option: %s',option)
    end
elseif(isequal(plottype,'plot0'))
    % plot type 3
    [h]=plot0(data,varargin{:});
    
    % are we deleting or keeping
    if(isequal(option,'keep'))
        while(button~=2)
            % bring plot to focus
            figure(h);
            handle=gca;
            
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
                        [clicked+0.5 clicked-0.5 clicked-0.5 clicked+0.5],[0 0.3 0]);
                    alpha(handles(clicked),0.99)
                end
            end
        end
        data=data(selected);
    elseif(isequal(option,'delete'))
        while(button~=2)
            % bring plot to focus
            figure(h);
            handle=gca;
            
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
                        [clicked+0.5 clicked-0.5 clicked-0.5 clicked+0.5],[0.3 0 0]);
                    alpha(handles(clicked),0.99);
                end
            end
        end
        data(selected)=[];
    else
        error('seizmo:selectrecs:badInput','Unknown option: %s',option)
    end
else
    error('seizmo:selectrecs:badInput','Unsupported plottype: %s',plottype)
end

end
