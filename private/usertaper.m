function [data,taptype,tap,tapopt,h,h2]=usertaper(data,skip,rdrift)
%USERTAPER    Interactively taper SAClab data

% check nargin
error(nargchk(1,3,nargin))

% check data structure
error(seischk(data,'dep'))

% defaults
if(nargin<3 || isempty(rdrift)); rdrift=1; end  % demean tapered data
if(nargin<2 || isempty(skip)); skip=1; end      % skip taper type/opt selection

% taperdefaults
taptype='hann';
tap=[0 0];
tapopt=[];

% length normalization
[b,e,npts,delta]=gh(data,'b','e','npts','delta');
data=ch(data,'b',0,'e',1,'delta',1./(npts-1));

% repeat until satisfied
h=-1; h2=-1;
while(1)
    % info
    mymenu={'Choose limits for the region you want untapered by clicking your';
            'mouse on the (soon to be) displayed figure.';
            '  ';
            'Usage:';
            '==================================';
            'Left clicking   -=> indicates the start of the untapered window';
            'Right clicking  -=> indicates the end of the untapered window';
            'Middle clicking -=> finalizes the tapers';
            '==================================';
            '  ';
            'You can rechoose taper parameters as many';
            'times as you wish until you finalize them.';
            '  ';
            'When you finalize there will be a confirmation';
            'menu at the end so you will have a chance to redo';
            'the tapering if your not satisfied.';
            '  ';
            '  ';
            'Choose a plot type for analysis:'};
    i=menu(mymenu,'Overlay Plot','Evenly Spaced Plot','Distance Spaced Plot','DO NOT TAPER','DIE!');
    
    % quick escape
    if (i==4)
        % NO TAPER
        data=ch(data,'b',b,'e',e,'delta',delta);
        return;
    elseif (i==5)
        % DEATH!
        error('user seppuku')
    
    % plot selection
    elseif (i==1)
        h=p2(data,'p2norm',true,'normmax',1);
    elseif (i==2)
        h=p0(data);
    else
        h=recsec(data);
    end
    
    % make taper bars
    span=ylim;
    tap=xlim;
    hold on
    goh(1)=plot([tap(1) tap(1)],span,'g','linewidth',4);
    goh(2)=plot([tap(2) tap(2)],span,'r','linewidth',4);
    hold off
    
    % loop until user is satisfied
    final=0;
    while (final==0)
        % focus plot and let user pick
        figure(h);
        [x,y,button]=ginput(1);
        
        % which mouse button
        if (button==1)
            % left - update untapered start
            tap(1)=x;
            figure(h);
            delete(goh(1));
            hold on
            goh(1)=plot([tap(1) tap(1)],span,'g','linewidth',4);
            hold off
        elseif (button==3)
            % right - update untapered end
            tap(2)=x;
            figure(h);
            delete(goh(2));
            hold on
            goh(2)=plot([tap(2) tap(2)],span,'r','linewidth',4);
            hold off
        elseif (button==2)
            % middle - finalize and click loop
            final=1;
        end
    end
    
    % fix tap(2)
    tap(2)=1-tap(2);
    
    % option to skip
    if(~skip)
        % taper types
        types={'barthannwin','bartlett','blackman','blackmanharris',...
            'bohmanwin','chebwin','flattopwin','gausswin','hamming',...
            'hann','kaiser','nuttallwin','parzenwin','rectwin',...
            'triang','tukeywin'};
        
        % choose taper type
        j=menu('CHOOSE A TAPER TYPE','DEFAULT','BARTHANN','BARTLETT',...
            'BLACKMAN','BLACKMAN-HARRIS','BOHMAN','CHEBYCHEV',...
            'FLAT TOP','GAUSSIAN','HAMMING','HANN','KAISER',...
            'NUTTALL','PARZEN','RECTANGULAR','TRIANGULAR','TUKEY');
        if(j>1)
            taptype=types{j-1};
        end
        
        % taper options
        if(j==7) % cheb
            tapopt=input('Chebychev - Stopband Attenuation? 0.1-200 [100]: ');
            if(isempty(tapopt))
                tapopt=100;
            end
        elseif(j==9) % gauss
            tapopt=input('Gaussian - Number of std dev? 0.1-10 [3.5]: ');
            if(isempty(tapopt))
                tapopt=3.5;
            end
        elseif(j==12) % kaiser
            tapopt=input('Kaiser - Stopband Attenuation Factor? 0-200 [7.5]: ');
            if(isempty(tapopt))
                tapopt=7.5;
            end
        elseif(j==17) % tukey
            tapopt=input('Tukey - Taper/Constant Ratio? 0-1 [1]: ');
            if(isempty(tapopt))
                tapopt=1;
            end
        end
    end
    
    % taper
    data2=taper(data,[tap(1) tap(2)],taptype,tapopt);
    
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
        h2=p0(data2);
    else
        h2=recsec(data2);
    end
    
    % TRY AGAIN?
    k=menu('Keep taper?','YES','NO - TRY AGAIN','NO - DIE!');
    if (k==1)
        % DONE
        data=ch(data2,'b',b,'e',e,'delta',delta);
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
