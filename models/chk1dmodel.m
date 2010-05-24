function [report]=chk1dmodel(model)
%CHK1DMODEL    Validate if a struct is as defined by PREM/AK135/IASP91
%
%    Usage:    msg=chk1dmodel(model)
%
%    Description: MSG=CHK1DMODEL(MODEL) checks that MODEL is a struct as
%     output from PREM/AK135/IASP91.  See those functions for details on
%     the struct layout.  MSG is an error structure if a problem is found
%     (otherwise it is empty).
%
%    Notes:
%
%    Examples:
%     Check that output from PREM is compatible:
%      mod=prem;
%      error(chk1dmodel(mod));
%
%    See also: PREM, AK135, IASP91

%     Version History:
%        May  24, 2010 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  24, 2010 at 15:40 GMT

% todo:

% check nargin
report=[];
error(nargchk(1,1,nargin));

% check model has the required fields
reqfields={'name' 'ocean' 'crust' 'isotropic' 'refperiod' 'flattened' ...
    'depth'};
if(~isstruct(model) || ~all(ismember(reqfields,fieldnames(model))))
    report.identifier='seizmo:chk1dmodel:dataNotStruct';
    report.message=['MODEL must be a struct with fields:\n' ...
        sprintf('%s ',reqfields{:})];
    return;
end

% check each required field
sz=size(model.name);
if(~ischar(model.name) || numel(sz)>2 || sz(1)~=1)
    report.identifier='seizmo:chk1dmodel:badName';
    report.message='The .name field of MODEL must be a string!';
    return;
elseif(~islogical(model.ocean) || ~isscalar(model.ocean))
    report.identifier='seizmo:chk1dmodel:badOcean';
    report.message='The .ocean field of MODEL must be TRUE or FALSE!';
    return;
elseif(~islogical(model.crust) || ~isscalar(model.crust))
    report.identifier='seizmo:chk1dmodel:badOcean';
    report.message='The .crust field of MODEL must be TRUE or FALSE!';
    return;
elseif(~islogical(model.isotropic) || ~isscalar(model.isotropic))
    report.identifier='seizmo:chk1dmodel:badOcean';
    report.message='The .isotropic field of MODEL must be TRUE or FALSE!';
    return;
elseif(~islogical(model.flattened) || ~isscalar(model.flattened))
    report.identifier='seizmo:chk1dmodel:badOcean';
    report.message='The .flattened field of MODEL must be TRUE or FALSE!';
    return;
elseif(~isreal(model.refperiod) || ~isscalar(model.refperiod) ...
        || model.refperiod<=0)
    report.identifier='seizmo:chk1dmodel:badRefPeriod';
    report.message='The .refperiod field of MODEL must be a scalar >0!';
    return;
end

% make sure all other fields are the same size as depth
sz=size(model.depth);
fields=setdiff(fieldnames(model),reqfields);
for i=1:numel(fields)
    if(~isequal(sz,size(model.(fields{i}))))
        report.identifier='seizmo:chk1dmodel:badPropSize';
        report.message=['The .' fields{i} ' field of MODEL must have ' ...
            'the same size as the .depth field!'];
        return;
    end
end

end
