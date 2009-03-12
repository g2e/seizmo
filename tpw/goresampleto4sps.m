
dates=dir('data-rinst/*');
dates(strcmp({dates.name},'.') | strcmp({dates.name},'..'))=[];

for i=1:numel(dates)
    disp(num2str(i))
    mkdir(['data-resample/' dates(i).name]);
    cd(['data-rinst/' dates(i).name]);
    data=r('*');
    cd(['../../data-resample/' dates(i).name]);
    data=syncrates(data,0.25);
    w(data);
    cd('../../');
end
