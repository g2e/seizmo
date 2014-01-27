function []=gooduglycheck(indir,outdir,varargin)
%GOODUGLYCHECK    Interactive QC tool for event data
%
%    Usage:    gooduglycheck(indir,outdir)
%              gooduglycheck(...,'freqbands',corners,...)
%              gooduglycheck(...,'moveouts',kmps,...)
%              gooduglycheck(...,'option',value,...)
%
%    Description:
%     GOODUGLYCHECK(INDIR,OUTDIR) is for inspection of record quality in
%     the directory structure under directory INDIR.  The program expects
%     that under directory INDIR there are directories corresponding to
%     separate events and that within these event directories are the
%     records to be QCed.  The user is presented with a series of menus and
%     interactive plots to select which events to look at and to delete
%     poor records from those events.  By default, 3 plots are produced
%     for each event displaying the recordings filtered at some frequency
%     band.  The filter ranges are 20-40s, 95-200s, and 35-120s which is
%     works for regional and global surface wave analysis.  The last plot
%     is interactive to allow record deletion -- use left-click to select
%     or deselect & middle-click to complete selection).  Records that are
%     not deleted are written within the directory OUTDIR with the same
%     directory layout.  
%
%     GOODUGLYCHECK(...,'FREQBANDS',CORNERS,...) indicates filtered bands
%     to be inspected (in addition to the unfiltered records).  The filters
%     are implemented as 2-pass 4-pole butterworth bandpasses.  The default
%     value of CORNERS is [1/40 1/20; 1/200 1/95].  This plots the period
%     ranges of 20-40s and 95-200s.  Note that CORNERS must be in Hz!
%
%     GOODUGLYCHECK(...,'MOVEOUTS',KMPS,...) indicates the moveouts to be
%     shown in the filtered record plots.  Moveouts should be in units of
%     kilometers per second.  The default value of KMPS is 3:5.
%
%     GOODUGLYCHECK(...,'OPTION',VALUE,...) passes additional options on to
%     PLOT1 for plot manipulation.
%
%    Notes:
%     - The directory structure should look as follows (the names are
%       allowed to be different ie. EVENTDIR1 may be 2006.044.04.03.55.9):
%        INDIR
%          |
%          --> EVENTDIR1
%          --> EVENTDIR2
%          .
%          .
%          .
%          --> EVENTDIRN
%                   |
%                   --> RECORD1
%                   --> RECORD2
%                   .
%                   .
%                   .
%                   --> RECORDN
%     - If OUTDIR exists the user is presented with the opportunity to
%       overwrite or delete the contents of OUTDIR or to exit the program.
%       DELETING WILL DELETE ONLY THE EVENT DIRECTORIES SELECTED _NOT_ THE
%       ENTIRE OUTDIR.  This allows for inplace quality control by setting
%       OUTDIR equal to INDIR.  Using overwrite in this case would not do
%       anything (except change the creation time).
%
%    Examples:
%     % Use only 1 filter band:
%     gooduglycheck('unchecked','qced','freqbands',[1/40 1/20])
%
%     % No filters, no moveouts (for speed):
%     gooduglycheck('unchecked','qced','fb',[],'m',[])
%
%    See also: FREQWINDOW, MAKEKERNELS, PLOTKERNELS

%     Version History:
%        Apr. 22, 2010 - major code cleanup and added documentation
%        Apr. 23, 2010 - lots of debugging
%        Jan. 19, 2011 - updated for all the recent gui changes
%        Jan. 21, 2011 - further tweaking
%        Apr.  5, 2011 - warn on event location variation
%        June  8, 2011 - turn on fill in cut commands
%        Apr.  2, 2012 - minor doc update
%        Jan. 26, 2014 - abs path exist fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan. 26, 2014 at 10:35 GMT

% todo
% - something like freqwindow may be better?

% check nargin
error(nargchk(2,inf,nargin));
if(mod(nargin,2))
    error('seizmo:gooduglycheck:uppairedOption',...
        'One (or more) input OPTION/VALUE is unpaired!');
end

% directory separator
fs=filesep;

% check indir
if(~isstring(indir))
    error('seizmo:gooduglycheck:badInput',...
        'INDIR must be a directory location!');
end
if(~isabspath(indir)); indir=[pwd fs indir]; end
if(~isdir(indir))
    error('seizmo:gooduglycheck:badInput',...
        'INDIR must be a directory location!');
end

% check outdir
reply='o';
if(~isstring(outdir))
    error('seizmo:gooduglycheck:badInput',...
        'OUTDIR must be a valid directory path!');
end
if(~isabspath(outdir)); outdir=[pwd fs outdir]; end
if(exist(outdir,'file') && ~isdir(outdir))
    error('seizmo:gooduglycheck:badInput',...
        'OUTDIR location is a file!');
elseif(isdir(outdir))
    fprintf('Directory: %s\nDirectory Exists!\n',outdir);
    reply=input('Overwrite/Delete/Quit? O/D/Q [Q]: ','s');
    if(strncmpi(reply,'o',1))
        disp(['Overwriting (But Not Deleting) ' ...
            'To-Be-Selected Event Directories!']);
    elseif(strncmpi(reply,'d',1))
        % only delete selected event directories
        disp('Deleting Contents of To-Be-Selected Event Directories!');
        
        % the code below deletes the entire superdirectory (too dangerous)
        %if(~rmdir(outdir,'s'))
        %    error('seizmo:gooduglycheck:couldNotDelete',...
        %        'Could Not Delete Directory: %s',outdir);
        %end
    else % quiting
        disp('Quiting!');
        return;
    end
end

% default options (pass through checks)
varargin=[{'m' 3:5 'fb' 1./[40 20; 200 95]} varargin];

% check optional inputs
nvargin=numel(varargin);
if(nvargin && ~iscellstr(varargin(1:2:end)))
    error('seizmo:gooduglycheck:badInput',...
        'OPTION must be a string!');
end
keep=true(nvargin,1);
for i=1:2:nvargin
    switch lower(varargin{i})
        case {'f' 'fb' 'freq' 'freqband' 'freqbands'}
            if(~isempty(varargin{i+1}) && (~isreal(varargin{i+1}) ...
                    || size(varargin{i+1},2)~=2 ...
                    || any(varargin{i+1}(:)<=0)))
                error('seizmo:gooduglycheck:badInput',...
                    ['FREQBANDS must be an Nx2 array of positive reals' ...
                    ' specifying ranges in Hz [LOW HIGH]!']);
            end
            freqbands=varargin{i+1};
            if(~isempty(freqbands))
                % force into [LOW HIGH] (so plot names look right)
                freqbands=[min(freqbands(:,1),freqbands(:,2)) ...
                    max(freqbands(:,1),freqbands(:,2))];
            end
            nfb=size(freqbands,1);
            keep(i:i+1)=false;
        case {'m' 'mo' 'move' 'moveout' 'moveouts'}
            if(~isreal(varargin{i+1}))
                error('seizmo:gooduglycheck:badInput',...
                    'MOVEOUTS must be positive reals!');
            end
            moveouts=varargin{i+1}(:);
            nmo=numel(moveouts);
            keep(i:i+1)=false;
    end
end
varargin=varargin(keep);

% get event directories
events=dir(indir);
events(strcmp({events.name},'.') | strcmp({events.name},'..'))=[];
events(~[events.isdir])=[];
eventlist=char({events.name}.');

% get user selected events
s=listdlg('PromptString','Select Events:',...
          'InitialValue',1:numel(events),...
          'ListSize',[170 300],...
          'ListString',eventlist);

% loop over user selected events
for i=s(:)'
    % display event at terminal so we know where we are
    disp(events(i).name)
    
    % read in data
    data0=readseizmo([indir fs events(i).name fs '*']);
    nrecs=numel(data0);
    
    % check event info matches
    [ev,outc,dist]=getheader(data0,'ev','o utc','dist');
    outc=cell2mat(outc);
    if(size(unique(ev,'rows'),1)>1)
        warning('seizmo:gooduglycheck:muddledHeader',...
            'Looks like EVENT location info varies among records!');
    elseif(any(abs(timediff(outc(1,:),outc))>0.002))
        error('seizmo:gooduglycheck:oUTCFieldVaries',...
            'ORIGIN time varies among records!');
    end
    
    % implement filters
    fdata=cell(nfb,1);
    for j=1:nfb
        fdata{j}=iirfilter(data0,'bp','b','c',freqbands(j,:),'o',4,'p',2);
    end
    
    % filter for deletion window
    ddata=iirfilter(data0,'bp','b','c',[1/120 1/35],'o',4,'p',2);
    
    % loop until user satisfied
    win=[min(dist/10) max(dist/2)];
    unsatisfied=true; ffh=-1;
    while(unsatisfied)
        % plot filtered records
        for j=1:nfb
            % p1 style plots
            ax=plot1(cut(fdata{j},win(1),win(2),'fill',true),...
                'xlabel',' ','ylabel',' ','align',true,varargin{:});
            if(ishandle(ax(1))); ffh(j)=get(ax(1),'parent'); end
            drawnow;
            
            % moveout lines -- not flags
            for k=1:nrecs
                hold(ax(k),'on');
                yrng=ylim(ax(k));
                xrng=dist(k)./moveouts;
                mh=plot(ax(k),xrng(:,[1 1])',yrng','linewidth',2);
                movekids(mh,'back');
                for m=1:nmo
                    text(xrng(m),yrng(2),num2str(moveouts(m)),...
                        'color',get(mh(m),'color'),'parent',ax(k));
                end
                hold(ax(k),'off');
            end
            
            % figure title
            set(ffh(j),'name',[events(i).name '  ' ...
                num2str(1/freqbands(j,2)) '-' ...
                num2str(1/freqbands(j,1)) 's']);
        end
        
        % use specific filter for record deletion
        ufh=figure('color','k','name',[events(i).name ' 35-120s']);
        ncols=fix(sqrt(nrecs));
        nrows=ceil(nrecs/ncols);
        ax=makesubplots(nrows,ncols,1:nrecs,'align','parent',ufh);
        [deleted,deleted]=selectrecords(cut(ddata,win(1),win(2),...
            'fill',true),'delete','p1',[],...
            'xlim',win,'ax',ax,'xlabel',' ','ylabel',' ',varargin{:});
        data=data0;
        data(deleted)=[];
        
        % get user action
        choice=0;
        while(~choice)
            mymenu={'Choose An Option:'};
            choice=menu(mymenu,...
                'Write and Continue To Next Event',...
                'Skip Event (No Write)',...
                'Redo Record Deletion',...
                'Redefine Time Limits',...
                'Exit');

            % handle user choice
            switch choice
                case 1 % write and proceed
                    unsatisfied=false;
                    % delete if indicated
                    if(strncmpi(reply,'d',1))
                        if(isdir([outdir fs events(i).name]))
                            if(~rmdir([outdir fs events(i).name],'s'))
                                error('seizmo:gooduglycheck:notDeleted',...
                                    'Could Not Delete Directory: %s',...
                                    [outdir fs events(i).name]);
                            end
                        end
                    end
                    writeseizmo(data,'path',[outdir fs events(i).name fs]);
                case 2 % skip
                    unsatisfied=false;
                case 3 % redo
                    % do nothing
                case 4 % change time
                    [win,win,ax]=userwindow(data0);
                    if(ishandle(ax)); close(get(ax,'parent')); end
                    win=win.limits;
                case 5 % exit
                    return;
            end
        end
        
        % close figure handles
        close(ffh(ishandle(ffh)));
        close(ufh(ishandle(ufh)));
    end
end

end
