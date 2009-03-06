function ax=mkgetaxes(axhandle);
% mkgetaxes..........determine dimensions of all objects in axes
%
% call: ax=mkgetaxes(axhandle);
%
%       axhandle: axes handle as returned by gca
%
% result: ax: four elements vector as returned by axis, containing
%             min and max x and y values
%             (2D plots only!!)
%
% I have encountered some problems with axis auto settings. This routine
% determines the minimum axis extensions that are necessary to encompass
% all drawn objects.
%
% Martin Knapmeyer, 05.12.2003


%%% init result
ax=axis;


%%% identify all visible objects
handles=findobj(axhandle,'visible','on');
anz=length(handles); % so many objects


%%% loop over all objects
%maxlist=[];
minx=inf;
maxx=-inf;
miny=inf;
maxy=-inf;
for indy=1:anz
    current=get(handles(indy),'type');
    switch current
        case {'text'}
            textent=get(handles(indy),'extent');
            minx=min(minx,textent(1));
            maxx=max(maxx,textent(1)+textent(3));
            miny=min(miny,textent(2));
            maxy=max(maxy,textent(2)+textent(4));
        case {'line','patch','surface','image'}
            xdata=get(handles(indy),'xdata');
            xdata=xdata(find((~isinf(xdata))&(~isnan(xdata))));
            ydata=get(handles(indy),'ydata');
            ydata=ydata(find((~isinf(ydata))&(~isnan(ydata))));
            if (~isempty(xdata))&(~isempty(ydata))
               minx=min(minx,min(min(xdata)));
               maxx=max(maxx,max(max(xdata)));
               %maxlist=[maxlist max(max(ydata))];
               miny=min(miny,min(min(ydata)));
               maxy=max(maxy,max(max(ydata)));
            end; % if ~isempty
        case {'axes'}
            % don't do anything - its the axes we want to resize!
        otherwise
            error(['MKGETAXES: don''t know what to to with type ' upper(current)]);
    end; % switch current
end; % for indy
%error('stop')

ax=[minx-1 maxx+1 miny-1 maxy+1];