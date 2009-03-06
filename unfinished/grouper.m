function [pop,big,small]=grouper(T)
%GROUPER    Counts groups

% global
global CONF

% NUMBER OF GROUPS
ng=max(T);

% COUNTERS AND DEFAULTS
pop=zeros(ng,1);
big=[];
small=[];
c1=0;
c2=0;

% LOOPING THROUGH EACH GROUP
for i=1:ng
    % COUNTING GROUP POPULATION
    pop(i)=length(find(i==T));
    
    % BIG/SMALL/SINGULAR GROUPS
    if (pop(i)>=CONF.MINSIG) 
        c1=c1+1; 
        big(c1)=i;
    elseif (pop(i)>2)
        c2=c2+1;
        small(c2)=i;
    end
end

end