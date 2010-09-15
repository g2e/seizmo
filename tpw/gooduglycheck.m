function [ok]=gooduglycheck(indir,outdir,varargin)
%GOODUGLYCHECK    Interactive QC tool for event data
%
%    Usage:    gooduglycheck(indir,outdir)
%              gooduglycheck(...,'freqbands',corners,...)
%              gooduglycheck(...,'moveouts',kmps,...)
%              ok=gooduglycheck(...)
%
%    Description: GOODUGLYCHECK(INDIR,OUTDIR) is used to inspect the
%     quality of records in the directory structure under directory INDIR.
%     The program expects that under directory INDIR there are directories
%     corresponding to separate events and that within these event
%     directories are the records to be QCed.  The user is presented with a
%     series of menus and interactive plots to select which events to
%     look at and to delete poor records from those events.  By default,
%     this is oriented to the quality control of surface wave recordings as
%     for each event 5 plots are produced displaying the unfiltered and
%     filtered recordings.  The filter ranges are 20-40s, 35-60s, 55-100s,
%     and 95-200s periods which are right in line for regional and global
%     surface wave analysis.  Records that are not deleted are written
%     within the directory OUTDIR with the same directory layout.  
%
%     GOODUGLYCHECK(...,'FREQBANDS',CORNERS,...) indicates filtered bands
%     to be inspected (in addition to the unfiltered records).  The filters
%     are implemented as 2-pass 4-pole butterworth bandpasses.  The default
%     value of CORNERS is [1/40 1/20; 1/60 1/35; 1/100 1/55; 1/200 1/95].
%     This gives 4 addtional bands looking at the period ranges of 20-40s,
%     35-60s, 55-100s, and 95-200s.  Note that CORNERS must be in Hz!
%
%     GOODUGLYCHECK(...,'MOVEOUTS',KMPS,...) indicates the moveouts to be
%     shown as markers in the plots.  Moveouts should be in units of
%     kilometers per second.  The default value of KMPS is 3:5.
%
%     OK=GOODUGLYCHECK(...) returns OK as FALSE if the user exited before
%     finishing the selected events or TRUE if all selected events were
%     processed.
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
%     Use only 1 filter band:
%      gooduglycheck('unchecked','qced','freqbands',[1/100 1/50])
%
%     No filters, no moveouts:
%      gooduglycheck('unchecked','qced','fb',[],'m',[])
%
%    See also: FREQWINDOW

%     Version History:
%        Apr. 22, 2010 - major code cleanup and added documentation
%        Apr. 23, 2010 - lots of debugging
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 23, 2010 at 10:35 GMT

% todo

% check nargin
error(nargchk(2,inf,nargin));
if(mod(nargin,2))
    error('seizmo:gooduglycheck:uppairedOption',...
        'One (or more) input OPTION/VALUE is unpaired!');
end

% did we finish ok? assume no at the start
ok=false;

% directory separator
fs=filesep;

% check indir
if(~ischar(indir) || size(indir,1)~=1 || ~isdir(indir))
    error('seizmo:gooduglycheck:badInput',...
        'INDIR must be a directory location!');
end

% check outdir
reply='o';
if(~ischar(outdir) || size(outdir,1)~=1)
    error('seizmo:gooduglycheck:badInput',...
        'OUTDIR must be a valid directory path!');
elseif(exist(outdir,'file') && ~isdir(outdir))
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
        
        % this deletes the entire superdirectory (too dangerous)
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
varargin=[{'m' 3:5 'fb' 1./[40 20; 60 35; 100 55; 200 95]} varargin];

% check optional inputs
nvargin=numel(varargin);
if(nvargin && ~iscellstr(varargin(1:2:end)))
    error('seizmo:gooduglycheck:badInput',...
        'OPTION must be a string!');
end
for i=1:2:numel(varargin)
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
            % force into [LOW HIGH] (so plot names look right)
            freqbands=[min(freqbands(:,1),freqbands(:,2)) ...
                max(freqbands(:,1),freqbands(:,2))];
            nfb=size(freqbands,1);
        case {'m' 'mo' 'move' 'moveout' 'moveouts'}
            if(~isreal(varargin{i+1}) || numel(varargin{i+1})>10)
                error('seizmo:gooduglycheck:badInput',...
                    'MOVEOUTS must be real numbers (up to 10 allowed)!');
            end
            moveouts=varargin{i+1};
            nmo=numel(moveouts);
        otherwise
            error('seizmo:gooduglycheck:unknownOption',...
                'Unknown Option: %s',varargin{i});
    end
end

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
    
    % check event info matches
    [ev,outc,dist,b,e]=getheader(data0,'ev','o utc','dist','b','e');
    outc=cell2mat(outc);
    if(size(unique(ev,'rows'),1)>1)
        error('seizmo:gooduglycheck:muddledHeader',...
            'EVENT location info varies among records!');
    elseif(any(abs(timediff(outc(1,:),outc))>0.002))
        error('seizmo:gooduglycheck:oUTCFieldVaries',...
            'ORIGIN time varies among records!');
    end
    
    % insert moveouts
    for j=0:(nmo-1)
        t=['t' num2str(j)];
        kt=['kt' num2str(j)];
        kmps=[num2str(moveouts(j+1)) 'km/s'];
        data0=changeheader(data0,t,dist/moveouts(j+1),kt,kmps);
    end
    
    % implement filters
    fdata=cell(nfb,1);
    for j=1:nfb
        fdata{j}=iirfilter(data0,'bp','b','c',freqbands(j,:),'o',4,'p',2);
    end
    
    % loop until user is satisfied
    user_happy=false; win=[min(b) max(e)];
    while(~user_happy)
        % plot up filtered data
        fh=nan(nfb,1); ax=fh;
        for j=1:nfb
            fh(j)=figure('name',[num2str(1/freqbands(j,2)) '-' ...
                num2str(1/freqbands(j,1)) 's  ' events(i).name],...
                'color','k');
            ax(j)=axes('parent',fh(j));
            ax(j)=plot1(fdata{j},'ax',ax(j),...
                'xlim',win,'markers',true);
        end
        
        % get user selection
        fh(j+1)=figure('name',['UNFILTERED  ' events(i).name],'color','k');
        ax(j+1)=axes('parent',fh(j+1));
        [data,deleted,ax(j+1)]=selectrecords(data0,'delete','p1',[],...
            'xlim',win,'markers',true);
        
        % get user action
        choice=0;
        while(~choice)
            mymenu={'Choose An Option'};
            choice=menu(mymenu,...
                'Write Kept Records and Proceed To Next Event',...
                'Skip Event (Write No Records)',...
                'Redo Record Selection',...
                'Redefine Time Limits',...
                'Exit');

            % handle user choice
            switch choice
                case 1 % write and proceed
                    user_happy=true;
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
                    user_happy=true;
                case 3 % redo
                    % do nothing
                case 4 % change time
                    [win,win,ax(j+2)]=userwindow(data0);
                    fh(j+2)=get(ax(j+2),'parent');
                    win=win.limits;
                case 5 % exit
                    return;
            end
        end
        
        % close figure handles
        close(fh(ishandle(fh)));
    end
end

% success!
ok=true;

end
