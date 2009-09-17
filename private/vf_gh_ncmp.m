function [value]=vf_gh_ncmp(varargin)
%VF_GH_NCMP    Returns value for virtual field NCMP

% just return 1s as this is only for single-component filetypes
value=ones(size(varargin{2},2),1);

end
