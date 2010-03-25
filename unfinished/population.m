function [pop]=population(grp)
%POPULATION    Returns cluster populations

pop=histc(grp.T,1:max(grp.T));

end
