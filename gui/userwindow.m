function [data,win,h,h2]=userwindow(data,fill,rdrift,push)
%USERWINDOW    Interactively window SEIZMO records

% check nargin
msg=nargchk(1,4,nargin);
if(~isempty(msg)); error(msg); end

% check data structure
msg=seizmocheck(data,'dep');
if(~isempty(msg)); error(msg.identifier,msg.message); end

% defaults
if(nargin<4 || isempty(push)); push=0; end          % keep timing
if(nargin<3 || isempty(rdrift)); rdrift=1; end      % demean new window
if(nargin<2 || isempty(fill)); fill=0; end          % no zero fill

% push all records to start at zero
if(push) 
    [b,e]=gh(data,'b','e'); 
    data=ch(data,'b',0,'e',e-b);
    iwin=[0 max(e-b)];
else
    iwin=[min(gh(data,'b')) max(gh(data,'e'))];
end

% repeat until satisfied
h=-1; h2=-1;
while(1)
    % info
    mymenu={'Choose a new data window by clicking your';
            'mouse on the (soon to be) displayed figure.';
            '  ';
            'Usage:';
            '==================================';
            'Left clicking   -=> indicates the window start position';
            'Right clicking  -=> indicates the window end position';
            'Middle clicking -=> finalizes the window positions';
            '==================================';
            '  ';
            'You can rechoose window parameters as many';
            'times as you wish until you finalize them.';
            '  ';
            'When you finalize there will be a confirmation';
            'menu at the end so you will have a chance to redo';
            'the windowing if your not satisfied.';
            '  ';
            'Choose a plot type for analysis:'};
    i=menu(mymenu,'Overlay Plot','Evenly Spaced Plot','Distance Spaced Plot','DO NOT WINDOW','DIE!');
    
    % quick escape
    if (i==4)
        % NO WINDOW
        win=[];
        if(push)
            data=ch(data,'b',b,'e',e); 
        end
        return;
    elseif (i==5)
        % DEATH!
        error('user seppuku')
    
    % plot selectionand view window
    elseif (i==1)
        h=p2(data,'p2norm',true,'normmax',1);
    elseif (i==2)
        h=p0(data,'xlimits',iwin);
    else
        h=recsec(data,'xlimits',iwin);
    end
    
    % make window bars
    span=ylim;
    win=xlim;
    hold on
    goh(1)=plot([win(1) win(1)],span,'g','linewidth',4);
    goh(2)=plot([win(2) win(2)],span,'r','linewidth',4);
    hold off
    
    % loop until user is satisfied
    final=0;
    while (final==0)
        % focus plot and let user pick
        figure(h);
        [x,y,button]=ginput(1);
        
        % which mouse button
        if (button==1)
            % left - update window start
            win(1)=x;
            figure(h);
            delete(goh(1));
            hold on
            goh(1)=plot([win(1) win(1)],span,'g','linewidth',4);
            hold off
        elseif (button==3)
            % right - update window end
            win(2)=x;
            figure(h);
            delete(goh(2));
            hold on
            goh(2)=plot([win(2) win(2)],span,'r','linewidth',4);
            hold off
        elseif (button==2)
            % middle - satisfied so clean up and exit click loop
            if (win(1)>win(2))
                % start and end reversed - exchange
                win=fliplr(win);
                figure(h);
                delete(goh(1));
                delete(goh(2));
                hold on
                goh(1)=plot([win(1) win(1)],span,'g','linewidth',4);
                goh(2)=plot([win(2) win(2)],span,'r','linewidth',4);
                hold off
            elseif (win(1)==win(2))
                % start/end the same
                win=xlim;
                figure(h);
                delete(goh(1));
                delete(goh(2));
                hold on
                goh(1)=plot([win(1) win(1)],span,'g','linewidth',4);
                goh(2)=plot([win(2) win(2)],span,'r','linewidth',4);
                hold off
            end
            final=1;
        end
    end
    
    % cut
    data2=cut(data,'z',win(1),win(2),'fill',fill,'filler',0);
    
    % remove trends
    if(rdrift==2)
        data2=rtr(data2);
    elseif(rdrift==1)
        data2=rmean(data2);
    end
    
    % view window
    if (i==1)
        h2=p2(data2,'p2norm',true,'normmax',1);
    elseif (i==2)
        h2=p0(data2,'xlimits',win);
    else
        h2=recsec(data2,'xlimits',win);
    end
    
    % TRY AGAIN?
    k=menu('Keep window?','YES','NO - TRY AGAIN','NO - DIE!');
    if (k==1)
        % DONE
        data=data2;
        if(push)
            [b2,e2]=gh(data,'b','e');
            data=ch(data,'b',b2+b,'e',e2+b); 
        end
        return;
    elseif (k==2)
        % REPEAT
        close([h h2]);
        h=-1; h2=-1;
    elseif (k==3)
        % DEATH!
        error('user seppuku')
    end
end

end
