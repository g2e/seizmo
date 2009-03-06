function [names,summary]=relarr_filterwrapper(data_init,CONF)

% FILTER DEPENDENT CONFIGURATION
CONF=filter_title(CONF);
CONF=filter_param(CONF);

% FILTER DATA
disp('FILTERING')
data_init=iirfilter(data_init,CONF.FILTTYPE,CONF.FILTSTYLE,CONF.FILTLIMITS,CONF.FILTORDER,CONF.FILTPASSES,CONF.FILTATTEN);

% POST-FILTER DATA PROCESSING
disp('PERFORMING POST-FILTER PROCESSING')
[data_init,CONF]=filter_postprep(data_init,CONF);

% QUICK SNR
if(CONF.QCKSNR)
    disp('CALCULATING SNR')
    % GET SNR FOR THE PHASE AND PUT IT IN THE HEADER
    data_init=ch(data_init,CONF.SNRFIELD,qcksnr(data_init,gh(data_init,CONF.TTFIELD),CONF.NOISWIN,CONF.SIGWIN));

    % DEFAULTS
    nrecs=length(data_init);
    snrcut=false(nrecs,1); nsnrcut=0;
    snrcutmark=ones(nrecs,1)*32; % ' '

    % SNR CUT
    if(CONF.QCKSNRCUT)
        % FIND CUT
        snrcut=(gh(data_init,CONF.SNRFIELD)<CONF.QCKSNRCUT);
        nsnrcut=nnz(snrcut);
        snrcutmark(snrcut)=88; % 'X'
    end

    % WRITE SNRCUT INFO
    if(exist([CONF.OUTPUTBASENAME '.qcksnr'],'file'))
        delete([CONF.OUTPUTBASENAME '.qcksnr']);
    end
    dlmwrite([CONF.OUTPUTBASENAME '.qcksnr'],...
        ['snr cut = ' num2str(CONF.QCKSNRCUT)],'');
    dlmwrite([CONF.OUTPUTBASENAME '.qcksnr'],...
        ['number cut = ' num2str(nsnrcut)],'-append','delimiter','');
    dlmwrite([CONF.OUTPUTBASENAME '.qcksnr'],' ','-append','delimiter','');
    dlmwrite([CONF.OUTPUTBASENAME '.window'],...
        'FILENAME     SNR','-append','delimiter','');
    dlmwrite([CONF.OUTPUTBASENAME '.qcksnr'],...
        [char({data_init.name}') ones(nrecs,5)*32 ...
        num2str(gh(data_init,CONF.SNRFIELD)) ones(nrecs,5)*32 ...
        char(snrcutmark)],'-append','delimiter','');

    % CUT
    data_init(snrcut)=[];
end

% FILENAMES
names={data_init.name}.';
    
% TRIAL LOOP
while(1)
    data=data_init;
    
    % WINDOW
    if(CONF.WINDOW)
        disp('WINDOWING DATA')
    
        % USER REFINED WINDOW
        if(CONF.USERWIN)
            % INITIAL WINDOW (ALIGNS AND FOCUSES ON PHASE) 
            data=cutim(data,CONF.TTFIELD,CONF.INTWIN(1),CONF.TTFIELD,CONF.INTWIN(2),CONF.FILL,CONF.FILLER);
            
            % ADJUST TIMING TO BE RELATIVE TO PHASE
            [b,e,tt]=gh(data,'b','e',CONF.TTFIELD);
            data=ch(data,'b',b-tt,'e',e-tt);
        
            % REMOVE MEAN
            data=rmean(data);
        
            % USER WINDOW
            [data,CONF.SIGWIN]=userwindow(data,1);
        
            % UPDATE CONFIG
            if(isempty(CONF.SIGWIN)); CONF.SIGWIN=CONF.INTWIN; end % CHECK FOR NO WINDOW
        % AUTO WINDOW (FOR THE BRAVE)
        else
            data=cutim(data,CONF.TTFIELD,CONF.SIGWIN(1),CONF.TTFIELD,CONF.SIGWIN(2),CONF.FILL,CONF.FILLER);
        
            % ADJUST TIMING TO BE RELATIVE TO PHASE
            [b,e,tt]=gh(data,'b','e',CONF.TTFIELD);
            data=ch(data,'b',b-tt,'e',e-tt);
        end
    
        % WRITE WINDOW INFO
        if(exist([CONF.OUTPUTBASENAME '.window'],'file'))
            delete([CONF.OUTPUTBASENAME '.window']);
        end
        dlmwrite([CONF.OUTPUTBASENAME '.window'],...
            ['interactive = ' num2str(CONF.USERWIN)],'');
        dlmwrite([CONF.OUTPUTBASENAME '.window'],...
            ['length = ' num2str(CONF.SIGWIN(2)-CONF.SIGWIN(1))],...
            '-append','delimiter','');
        dlmwrite([CONF.OUTPUTBASENAME '.window'],' ','-append','delimiter','');
        dlmwrite([CONF.OUTPUTBASENAME '.window'],...
            'FILENAME     BEGIN     END    DELTA     NPTS',...
            '-append','delimiter','');
        nrecs=length(data); pad=ones(nrecs,5)*32;
        dlmwrite([CONF.OUTPUTBASENAME '.window'],...
            [char({data.name}') pad ...
            num2str(b) pad ...
            num2str(e) pad ...
            num2str(gh(data,'delta')) pad ...
            num2str(gh(data,'npts'))],'-append','delimiter','');
    end

    % TAPER
    if(CONF.TAPER)
        disp('TAPERING DATA')

        % USER REFINED TAPER
        if(CONF.USERTAPER)
            [data,CONF.TAPERTYPE,CONF.TAPERHW,CONF.TAPEROPT]=usertaper(data);
        % AUTO TAPER (FOR THE BRAVE)
        else
            data=taper(data,CONF.TAPERHW,CONF.TAPERTYPE,CONF.TAPEROPT);
        end

        % WRITE TAPER INFO
        if(exist([CONF.OUTPUTBASENAME '.taper'],'file'))
            delete([CONF.OUTPUTBASENAME '.taper']);
        end
        dlmwrite([CONF.OUTPUTBASENAME '.taper'],...
            ['interactive = ' num2str(CONF.USERTAPER)],'');
        dlmwrite([CONF.OUTPUTBASENAME '.taper'],...
            ['taper type = ' num2str(CONF.TAPERTYPE)],...
            '-append','delimiter','');
        dlmwrite([CONF.OUTPUTBASENAME '.taper'],...
            ['taper halfwidths = ' num2str(CONF.TAPERHW)],...
            '-append','delimiter','');
        dlmwrite([CONF.OUTPUTBASENAME '.taper'],...
            ['taper option = ' num2str(CONF.TAPEROPT)],...
            '-append','delimiter','');
    end
    error('user sepuku')
    % GET RELATIVE ARRIVALS
    summary=relarr_core(data,CONF);
    
    % ANOTHER TRIAL?
    response=menu('DO THE ALIGNMENT OVER?','YES - REDO','YES - USE THESE RESULTS','NO');
    if(response==3)
        % BREAK OUT OF THE LOOP
        break;
    elseif(response==2)
        % UPDATE TIMING FOR NON-NAN ENTRIES
        notnans=~isnan(summary(2,:)); arrivals=summary(2,:);
        data_init=ch(data_init(notnans),CONF.TTFIELD,arrivals(notnans));
    end
end

end