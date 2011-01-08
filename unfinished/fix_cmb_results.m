function [results]=fix_cmb_results()

% get all results.mat files in this directory
list=onefilelist('*_outliers_results.mat');

% loop over each result, add correction, write
for i=1:numel(list)
    % read in result
    disp(list(i).name);
    tmp=load(list(i).name);
    
    % add corrections
    if(~isempty(tmp.useralign))
        tmp.corrections=cmb_corrections(tmp.phase,tmp.useralign.data);
    else
        tmp.corrections=[];
    end
    
    % save result
    save(list(i).name,'-struct','tmp');
    
    % add to results matrix
    results(i)=tmp;
end

end
