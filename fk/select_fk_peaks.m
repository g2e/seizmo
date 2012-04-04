function [peaks,ax]=select_fk_peaks(map,varargin)
%SELECT_FK_PEAKS    Interactively select peaks in a fk beamformer plot
%
%    Usage:    peaks=select_fk_peaks(map)
%              peaks=select_fk_peaks(map,...)
%
%    Description:
%     PEAKS=SELECT_FK_PEAKS(MAP) presents a fk beamformer plot that can be
%     interacted with by the user for selecting peaks.  Info on the peaks
%     is returned in PEAKS.  Interaction is fairly simple (mouse based):
%      left-click                 - select peak bounding box corner 1 or 2
%      right-click (ctrl-click)   - delete last bounding box corner
%      middle click (shift-click) - return
%     The output PEAKS is a struct with the following fields:
%      .xrng     - bounding box xlimits
%      .xtype    - xlimit type ('east' or 'baz')
%      .yrng     - bounding box ylimits
%      .ytype    - ylimit type ('north' or 'slow')
%      .db       - dB of peak
%      .backazi  - backazimuth of peak
%      .horzslow - slowness of peak
%      .meddb    - median dB of beamformer image (same for all peaks)
%     The user may bound as many peaks as they like (until they middle-
%     click the plot).  Only one peak is found per bounding box.  The PEAKS
%     output fields will expand in the number of rows for each new peak.
%
%     PEAKS=SELECT_FK_PEAKS(MAP,...) passes additional inputs to PLOTFKMAP.
%     This is useful for customizing the plotting.
%
%     [PEAKS,AX]=SELECT_FK_PEAKS(...) returns the handle to the axes that
%     the beamformer spectra map was plotted in as AX.
%
%    Notes:
%     - Closing the plot is the same as middle clicking.
%
%    Examples:
%     % Get the frequency-slowness of some array data and pick the peaks:
%     map=fkmap(data,50,201,[1/20 1/10],true);
%     peaks=select_fk_peaks(map,[0 6],'median');
%
%    See also: PLOTFKMAP, FKDBINFO, FKMAP, FKVOLUME, FKVOL2MAP

%     Version History:
%        Feb. 16, 2011 - initial version
%        Mar. 29, 2012 - ax output, minor struct change, added example,
%                        fix backazimuth bug for polar case
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 29, 2012 at 10:35 GMT

% todo:

% check nargin
error(nargchk(1,inf,nargin));

% create plot
ax=plotfkmap(map,varargin{:});

% polar plot?
polar=map.polar;
if(polar)
    xtype='baz';
    ytype='slow';
else % cartesian
    xtype='east';
    ytype='north';
end

% background color?
cmap=get(get(ax,'parent'),'colormap');
fgcolor=cmap(end,:);
bgcolor=cmap(1,:);

% get median info
[peak,med]=fkdbinfo(map);

% allocate output
peaks=struct('xrng',[],'xtype',xtype,'yrng',[],'ytype',ytype,...
    'db',[],'backazi',[],'horzslow',[],'meddb',med.db);

% check axes
if(~isscalar(ax) || ~ishandle(ax) || ~strcmp('axes',get(ax,'type')) ...
        || ~strcmp('fkmap',get(ax,'tag')))
    return;
end

% make selections
none=true; cnt=0;
peakh=[];
while(true)
    % get user input
    try
        % force axes to front first
        axes(ax);
        [x,y,button]=ginput(1);
    catch
        % plot closed so finish
        return;
    end
    
    % require ax
    if(gca~=ax); continue; end
    
    % act by button
    switch button
        case 1 % left or ctrl+shift
            % 1st or 2nd corner?
            if(none) % 1st = show point
                % save point
                bbx(1)=x;
                bby(1)=y;
                
                % draw + on image
                hold(ax,'on');
                pnth=plot(ax,x,y,'+','color',fgcolor);
                hold(ax,'off');
                
                % increment
                cnt=cnt+1;
                none=false;
            else % 2nd = draw box & get peak info
                % save point
                bbx(2)=x;
                bby(2)=y;
                
                % remove + sign
                delete(pnth);
                
                % draw box
                hold(ax,'on');
                if(polar)
                    % convert cartesian box corners to polar
                    [bbx,bby]=cart2pol(bby,bbx);
                    bbx=bbx*180/pi;
                    
                    % box must be < 180 degrees
                    if(diff(bbx)>180); bbx(2)=bbx(2)-360; end
                    if(diff(bbx)<-180); bbx(2)=bbx(2)+360; end
                    
                    % draw radial box lines
                    bbx2=[bbx(1:2); bbx(1:2)];
                    bby2=[bby([1 1]); bby([2 2])];
                    [bby2,bbx2]=pol2cart(bbx2*pi/180,bby2);
                    boxh(cnt,1:2)=plot(ax,bbx2,bby2,...
                        'color',fgcolor,'linewidth',2);
                    
                    % draw arc box lines
                    bbx2=[linspace(bbx(1),bbx(2),30)' ...
                        linspace(bbx(1),bbx(2),30)'];
                    bby2=bby([ones(30,1) 2*ones(30,1)]);
                    [bby2,bbx2]=pol2cart(bbx2*pi/180,bby2);
                    boxh(cnt,3:4)=plot(ax,bbx2,bby2,...
                        'color',fgcolor,'linewidth',2);
                else % cartesian
                    boxh(cnt,1:4)=plot(ax,...
                        [bbx([1 1 2 2]); bbx([1 2 2 1])],...
                        [bby([1 1 1 2]); bby([2 1 2 2])],...
                        fgcolor,'linewidth',2);
                end
                hold(ax,'off');
                
                % get peak info
                if(polar)
                    peak=fkdbinfo(map,[],sort(bbx),sort(bby));
                else
                    peak=fkdbinfo(map,[],[],[],sort(bbx),sort(bby));
                end
                peaks.db(cnt,1)=peak.db;
                peaks.backazi(cnt,1)=peak.backazi;
                peaks.horzslow(cnt,1)=peak.horzslow;
                peaks.xrng(cnt,:)=bbx;
                peaks.yrng(cnt,:)=bby;
                
                % print summary
                fprintf(['\nPEAK #' num2str(cnt) ' VALUES:\n']);
                fprintf('Back Azimuth:  %g\n',peak.backazi);
                fprintf('Horz Slowness: %g\n',peak.horzslow);
                fprintf('Abs Decibels:  %g\n',peak.db);
                fprintf('Rel Decibels:  %g\n',peak.db-med.db);
                
                % mark peak
                hold(ax,'on');
                [bby2,bbx2]=pol2cart(pi/180*peak.backazi,peak.horzslow);
                peakh(cnt)=plot(ax,bbx2,bby2,'x','color',bgcolor);
                hold(ax,'off');
                
                % no points
                none=true;
            end
        case 2 % middle or shift
            return;
        case 3 % right or ctrl
            % undo last action
            if(~cnt) % skip if nothing done
                continue;
            elseif(none) % delete last peak & bounding box
                % delete box & peak
                delete([boxh(cnt,:) peakh(cnt)]);
                
                % restore box corner 1
                bbx=peaks.xrng(cnt,1);
                bby=peaks.yrng(cnt,1);
                
                % fix box corner 1 if polar
                if(polar); [bby,bbx]=pol2cart(pi/180*bbx,bby); end
                
                % remove peak info
                peaks.xrng(cnt,:)=[];
                peaks.yrng(cnt,:)=[];
                peaks.db(cnt)=[];
                peaks.backazi(cnt)=[];
                peaks.horzslow(cnt)=[];
                
                % restore + on image for box corner 1
                hold(ax,'on');
                pnth=plot(ax,bbx,bby,'+','color',fgcolor);
                hold(ax,'off');
                
                % there is 1 pnt
                none=false;
            else % delete point & decrement
                delete(pnth);
                cnt=cnt-1;
                none=true;
            end
        otherwise
            key2zoompan(button,ax);
    end
end

end
