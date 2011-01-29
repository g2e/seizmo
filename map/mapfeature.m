function []=mapfeature(f,type,ax)
axis(ax);
hold(ax,'on');
cmap=hsv(numel(f));
switch lower(type)
    case 'point'
        for i=1:numel(f)
            plot(f(i).longitude,f(i).latitude,'.','color',cmap(i,:),'parent',ax);
            plot(f(i).longitude-360,f(i).latitude,'.','color',cmap(i,:),'parent',ax);
            plot(f(i).longitude+360,f(i).latitude,'.','color',cmap(i,:),'parent',ax);
        end
    case 'line'
        for i=1:numel(f)
            f(i).longitude=unwrap(f(i).longitude*pi/180)*180/pi;
            plot(f(i).longitude,f(i).latitude,'color',cmap(i,:),'parent',ax);
            plot(f(i).longitude-360,f(i).latitude,'color',cmap(i,:),'parent',ax);
            plot(f(i).longitude+360,f(i).latitude,'color',cmap(i,:),'parent',ax);
        end
    case 'patch'
        for i=1:numel(f)
            f(i).longitude=unwrap(f(i).longitude*pi/180)*180/pi;
            patch(f(i).longitude,f(i).latitude,cmap(i,:),'parent',ax);
            patch(f(i).longitude-360,f(i).latitude,cmap(i,:),'parent',ax);
            patch(f(i).longitude+360,f(i).latitude,cmap(i,:),'parent',ax);
        end
end
hold(ax,'off');
end
