function [ok]=gooduglycheck(indir,outdir)
% did we finish ok?
ok=false;

% check indir
if(~ischar(indir) || size(indir,1)~=1)
    error('INDIR must be a string giving one directory!');
elseif(~isdir(indir))
    error('INDIR must be a directory!');
end

% check outdir
if(~ischar(outdir) || size(outdir,1)~=1)
    error('OUTDIR must be a string giving one directory!');
elseif(exist(outdir,'file') && ~isdir(outdir))
    error('OUTDIR must be a directory!');
end

% get date directories
dates=dir(indir);
dates(strcmp({dates.name},'.') | strcmp({dates.name},'..'))=[];
dates(~[dates.isdir])=[];
datelist=char(strcat({dates.name}.'));

% get user selected start date
s=listdlg('PromptString','Select events:',...
          'InitialValue',1:numel(dates),...
          'ListSize',[170 300],...
          'ListString',datelist);

% loop over user selected events
cwd=pwd;
for i=s(:)'
    disp(dates(i).name)
    
    % get data
    cd([indir filesep dates(i).name]);
    data0=r('*');
    cd(cwd);
    
    % insert moveouts
    dist=gh(data0,'dist');
    data0=ch(data0,...
        't0',dist/3,'kt0','3km/s',...
        't1',dist/4,'kt1','4km/s');
    
    % filter data at 4 different bands
    % what bands? 20s-40s, 35s-60s, 55s-100s, 95s-200s
    data1=iirfilter(data0,'bandpass','butter',[1/40 1/20],4,2);
    data2=iirfilter(data0,'bandpass','butter',[1/60 1/35],4,2);
    data3=iirfilter(data0,'bandpass','butter',[1/100 1/55],4,2);
    data4=iirfilter(data0,'bandpass','butter',[1/200 1/95],4,2);
    
    % loop until user is satisfied
    still_looking=1; xwin=[];
    while(still_looking)
        % plot up filtered data
        f4=plot1(data4,...
            'name',[dates(i).name '  Period band: 95-200s'],...
            'xlimits',xwin,'position',[0 0 1 1]);
        f3=plot1(data3,...
            'name',[dates(i).name '  Period band: 55-100s'],...
            'xlimits',xwin,'position',[0 0 1 1]);
        f2=plot1(data2,...
            'name',[dates(i).name '  Period band: 35-60s'],...
            'xlimits',xwin,'position',[0 0 1 1]);
        f1=plot1(data1,...
            'name',[dates(i).name '  Period band: 20-40s'],...
            'xlimits',xwin,'position',[0 0 1 1]);
        %f3=plot1(data3,...
        %    'name',[dates(i).name '  Period band: 55-100s'],'xlimits',xwin);
        %set(f3,'units','normalized','position',[0 0 0.48 0.42]);
        %f4=plot1(data4,...
        %    'name',[dates(i).name '  Period band: 95-200s'],'xlimits',xwin);
        %set(f4,'units','normalized','position',[0.5 0 0.48 0.42]);
        %f1=plot1(data1,...
        %    'name',[dates(i).name '  Period band: 20-40s'],'xlimits',xwin);
        %set(f1,'units','normalized','position',[0 0.5 0.48 0.42]);
        %f2=plot1(data2,...
        %    'name',[dates(i).name '  Period band: 35-60s'],'xlimits',xwin);
        %set(f2,'units','normalized','position',[0.5 0.5 0.48 0.42]);
        
        % get user selection
        [data,f5,f5]=selectrecords(data0,'delete','p1',...
            'name',[dates(i).name '  Period band: UNFILTERED'],...
            'xlimits',xwin,'position',[0 0 1 1]);
        
        % get user action
        mymenu={'What to do?'};
        j=menu(mymenu,...
            'Move to next event',...
            'Redo deletion',...
            'Redefine xlimits',...
            'Skip event (event is bad)',...
            'QUIT!');
        if(j==5); return; end
        if(j==1 || j==4); still_looking=0; end
        close([f1 f2 f3 f4 f5])
        
        % redefine xlimits if desired
        if(j==3)
            [xwin,xwin,f1,f2]=userwindow(data0);
            close([f1 f2]);
        else
            xwin=[];
        end
    end
    
    % check if user wanted to skip this event or 0 files left
    if(j==4); continue; end
    
    % make output directory, write data, move back to top
    mkdir([outdir filesep dates(i).name]);
    cd([outdir filesep dates(i).name]);
    
    % write data
    w(data);
    cd(cwd);
end

ok=true;

end