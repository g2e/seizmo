function [ok]=freqwindow(indir,outdir)
%window surface wave data at several frequencies
%
% how to adjust for your data
% - learn filter_bank and adjust line 42ish
% - adjust 1D models (prem to array, custom in array) at line 50ish
% - change taper parameter below
% - change SNR cutoff/window parameters below
% - change final padded start/stop below
taperwidth=0.2; % taper 20% on each edge
snrcut=3; % trim records with signal to noise ratio < 3
padstart=-2000; % start padding at -2000s
padend=7000; % pad to 7000s
snrwin=0.5; % size of noise window relative to signal

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

% get filter bank
bank=flipud(filter_bank([0.0055 0.055],'variable',0.2,0.1));
nf=size(bank,1);
fs=num2str((1:nf)','%02d');
band=char(strcat({'   '},num2str(1./bank(:,2),'%5.1f'),{'s  to  '},...
    num2str(1./bank(:,3),'%5.1f'),'s'));

% initial moveout models
prem=[ % from a random paper
 15.0000    3.8600
 20.0000    3.9650
 25.0000    3.9900
 30.0000    4.0120
 35.0000    4.0220
 40.0000    4.0300
 45.0000    4.0350
 50.0000    4.0400
 55.0000    4.0450
 60.0000    4.0500
 65.0000    4.0580
 70.0000    4.0670
 75.0000    4.0770
 80.0000    4.0880
 85.0000    4.1010
 90.0000    4.1150
 95.0000    4.1330
100.0000    4.1500
105.0000    4.1660
110.0000    4.1820
115.0000    4.1980
120.0000    4.2140
125.0000    4.2300
130.0000    4.2460
135.0000    4.2620
140.0000    4.2780
145.0000    4.2940
150.0000    4.3100
155.0000    4.3260
160.0000    4.3420
165.0000    4.3580
170.0000    4.3740
175.0000    4.3900
180.0000    4.4060
185.0000    4.4220];
local=[ % simple fit to 1D results
 18.5000    3.6000
 35.0000    3.8800
170.0000    4.3500];

% get moveout for each filter (linear interpolation from above)
moveto=interp1(prem(:,1),prem(:,2),1./bank(:,1),[],'extrap');
movein=interp1(local(:,1),local(:,2),1./bank(:,1),[],'extrap');
adj=[1.5 1.25 1.1 1.05 1.01 .99 .95 .9 .75 .5 1];

% loop over user selected events
cwd=pwd;
for i=s(:)'
    % display event name
    disp(dates(i).name);
    
    % get data
    cd([indir filesep dates(i).name]);
    data1=sortbyfield(r('*'),'gcarc');
    cd(cwd);
    nrecs=numel(data1);
    idx=1:nrecs;
    
    % get some header info
    [dist,b,e,o,t,iztype]=gh(data1,'dist','b','e','o','t','iztype');
    mind=min(dist);
    
    % loop over each filter (short period to long)
    for j=1:nf
        % display filter
        disp(band(j,:));
        
        % filter data
        data2=iirfilter(data1,'bandpass','butter',bank(j,2:3),4,2);
        
        % loop until user is happy with this band
        ssatisfied=0;
        while(~ssatisfied)
        
        % make temp variables
        a=max(1000,1/bank(j,1)*20);
        iwin=[-a a];
        mvto=moveto(j);
        mvin=movein(j);
        good=idx;
        f2=[];
        
        % loop until user is happy with initial window/winnow
        redo=0;
        skipthis=0;
        skiprest=0;
        satisfied=0;
        while(~satisfied)
            % adjust timing of data for moveout
            z=o+mind./mvto+(dist-mind)./mvin;
            data3=timeshift(data2(good),-z(good));
            data3=ch(data3,'iztype','iunkn');
            
            % initial window about arrival
            data3=cut(data3,iwin(1),iwin(2));
            
            % make record section plot
            f1=recsec(data3,'xlimits',iwin,...
                'name',[dates(i).name '  Period band: ' band(j,:)]);
            
            % get choice from user
            choice=menu('what to do?',...
                'Window Phase',...
                ['Adjust moveout to array (Currently ' num2str(mvto) 'km/s'],...
                ['Adjust moveout within array (Currently ' num2str(mvin) 'km/s'],...
                'Adjust initial window',...
                'Remove some records',...
                'Redo this filter band',...
                'Skip this filter band',...
                'Skip all remaining filter bands',...
                'QUIT!!!');
            
            % act on choice
            switch choice
                case 1 % window phase
                    satisfied=1;
                case 2 % moveout to array
                    choice=menu(['ADJUST MOVEOUT TO ARRAY?  (CURRENTLY '...
                        num2str(mvto) 'km/s)'],...
                        ['+50% = ' num2str(mvto*1.50) 'km/s)'],...
                        ['+25% = ' num2str(mvto*1.25) 'km/s)'],...
                        ['+10% = ' num2str(mvto*1.10) 'km/s)'],...
                        [' +5% = ' num2str(mvto*1.05) 'km/s)'],...
                        [' +1% = ' num2str(mvto*1.01) 'km/s)'],...
                        [' -1% = ' num2str(mvto*0.99) 'km/s)'],...
                        [' -5% = ' num2str(mvto*0.95) 'km/s)'],...
                        ['-10% = ' num2str(mvto*0.90) 'km/s)'],...
                        ['-25% = ' num2str(mvto*0.75) 'km/s)'],...
                        ['-50% = ' num2str(mvto*0.50) 'km/s)'],...
                        'KEEP AS IS!');
                    mvto=mvto*adj(choice);
                case 3 % moveout within array
                    choice=menu(['ADJUST MOVEOUT WITHIN ARRAY?  (CURRENTLY '...
                        num2str(mvin) 'km/s)'],...
                        ['+50% = ' num2str(mvin*1.50) 'km/s)'],...
                        ['+25% = ' num2str(mvin*1.25) 'km/s)'],...
                        ['+10% = ' num2str(mvin*1.10) 'km/s)'],...
                        [' +5% = ' num2str(mvin*1.05) 'km/s)'],...
                        [' +1% = ' num2str(mvin*1.01) 'km/s)'],...
                        [' -1% = ' num2str(mvin*0.99) 'km/s)'],...
                        [' -5% = ' num2str(mvin*0.95) 'km/s)'],...
                        ['-10% = ' num2str(mvin*0.90) 'km/s)'],...
                        ['-25% = ' num2str(mvin*0.75) 'km/s)'],...
                        ['-50% = ' num2str(mvin*0.50) 'km/s)'],...
                        'KEEP AS IS!');
                    mvin=mvin*adj(choice);
                case 4 % initial window
                    choice=menu(['ADJUST INITIAL WINDOW?  (CURRENTLY '...
                        num2str(diff(iwin)) 's)'],...
                        ['+50% = ' num2str(diff(iwin)*1.50) 's)'],...
                        ['+25% = ' num2str(diff(iwin)*1.25) 's)'],...
                        ['+10% = ' num2str(diff(iwin)*1.10) 's)'],...
                        [' +5% = ' num2str(diff(iwin)*1.05) 's)'],...
                        [' +1% = ' num2str(diff(iwin)*1.01) 's)'],...
                        [' -1% = ' num2str(diff(iwin)*0.99) 's)'],...
                        [' -5% = ' num2str(diff(iwin)*0.95) 's)'],...
                        ['-10% = ' num2str(diff(iwin)*0.90) 's)'],...
                        ['-25% = ' num2str(diff(iwin)*0.75) 's)'],...
                        ['-50% = ' num2str(diff(iwin)*0.50) 's)'],...
                        'KEEP AS IS!');
                    iwin=iwin*adj(choice);
                case 5 % select
                    [data3,bad,f2]=selectrecords(data3,'delete','p0',...
                        'name',[dates(i).name '  Period band: ' band(j,:)]);
                    good(bad)=[];
                case 6 % redo
                    redo=1;
                    satisfied=1;
                case 7 % skip one
                    skipthis=1;
                    satisfied=1;
                case 8 % skip all
                    skiprest=1;
                    satisfied=1;
                case 9 % quit
                    return;
            end
            
            % close figures
            try
                close([f1 f2]);
                f2=[];
            catch
            end
        end
        
        % check redo/skip
        if(redo); continue; end
        if(skiprest || skipthis); break; end
        
        % adjust timing of data for moveout
        z=o+mind./mvto+(dist-mind)./mvin;
        data2=timeshift(data2,-z);
        data2=ch(data2,'iztype','iunkn');
        
        % initial window about arrival
        data2=cut(data2,iwin(1),iwin(2));
        
        % loop until user satisfied with final window/winnow
        f3=[];
        satisfied=0;
        while(~satisfied)
            % window
            [data3,win,f1,f2]=userwindow(data2(good));
            
            % insert window limits into header
            data3=ch(data3,'a',win(1),'ka','winbgn',...
                'f',win(2),'kf','winend');
            
            % taper
            data3=taper(data3,taperwidth);
            
            % snr cut (both sides!!)
            bad1=snrcut>quicksnr(data2(good),[win(1)-snrwin*(win(2)-win(1)) win(1)]-eps,win);
            bad2=snrcut>quicksnr(data2(good),[win(2) win(2)+snrwin*(win(2)-win(1))]+eps,win);
            data3(bad1 | bad2)=[];
            good(bad1 | bad2)=[];
            
            % plot tapered and trimmed set
            f4=recsec(data3,'xlimits',win);
            
            % undo time shift
            data3=timeshift(data3,z(good));
            data3=ch(data3,'iztype',iztype(good));
            
            % pad
            data3=cut(data3,padstart,padend,'fill',true);
            
            % plot padded data
            f5=recsec(data3,'xlimits',[padstart padend]);
            
            % get choice from user
            choice=menu('what to do?',...
                'Write files and move on!',...
                'Redo final window (w/ taper and snr cut)',...
                'Remove some records, redo final window',...
                'Redo this filter band',...
                'Skip this filter band',...
                'Skip all remaining filter bands',...
                'QUIT!!!');
            
            % act on choice
            switch choice
                case 1 % write files
                    % make output directory, move into it
                    mkdir([outdir filesep dates(i).name filesep fs(j,:)]);
                    cd([outdir filesep dates(i).name filesep fs(j,:)]);
                    
                    % write data, move back to top
                    w(data3);
                    cd(cwd);
                    satisfied=1;
                    ssatisfied=1;
                case 2 % redo final window
                    
                case 3 % remove some records
                    [data3,bad,f3]=selectrecords(data3,'delete','p0',...
                        'name',[dates(i).name '  Period band: ' band(j,:)]);
                    good(bad)=[];
                case 4 % redo initial and final window
                    satisfied=1;
                case 5 % skip one
                    satisfied=1;
                    ssatisfied=1;
                case 6 % skip all
                    skiprest=1;
                    satisfied=1;
                    ssatisfied=1;
                case 7 % quit
                    return;
            end
            
            % close figures
            try
                close([f1 f2 f3 f4 f5]);
                f3=[];
            catch
            end
        end
        
        end
        
        % check skip (skip rest of event)
        if(skiprest); break; end
    end
end

ok=true;

end
