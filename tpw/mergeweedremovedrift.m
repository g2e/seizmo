
dates=dir('data-hdrfix/*');
dates(strcmp({dates.name},'.') | strcmp({dates.name},'..'))=[];

for i=1:numel(dates)
    disp(num2str(i))
    mkdir(['data-mergeweedremovedrift/' dates(i).name]);
    cd(['data-hdrfix/' dates(i).name]);
    data=r('*');
    cd(['../../data-mergeweedremovedrift/' dates(i).name]);
    data=merge(data);
    [b,e]=gh(data,'b','e');
    short=(e-b)<(0.7*max(e-b));
    data=rtr(data(~short));
    w(data);
    cd('../../');
end
