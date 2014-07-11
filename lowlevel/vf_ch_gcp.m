function [head]=vf_ch_gcp(def,head,value)
%VF_CH_GCP    Sets virtual field GCP

% convert nan/inf to undef value
value(isnan(value) | isinf(value))=def.undef.ntype;

% get defined
good=(value~=def.undef.ntype);

% convert gcp to baz if defined
if(any(good))
    value(good)=mod(value(good)+180,360);
end

% find baz header field position
for i=1:numel(def.real)
    if(isfield(def.real(i).pos,'baz'))
        baz=def.real(i).pos.baz;
        break;
    end
end

% set header
head(baz,:)=value.';

end
