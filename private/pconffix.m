function [conf]=pconffix(conf)

% 
fields=fieldnames(conf.GENERAL).';
for i=fields
    for j=conf.GENERAL.(i{:})
        if(isempty(conf.(j{:}))); conf.(j{:})=conf.(i{:}); end
    end
end

end