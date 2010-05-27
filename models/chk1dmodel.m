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
%        May  27, 2010 - now allows arrays of models
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated May  27, 2010 at 15:05 GMT

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

% loop over each model
fields=setdiff(fieldnames(model),reqfields);
for i=1:numel(model)
    % check each required field
    sz=size(model(i).name);
    if(~ischar(model(i).name) || numel(sz)>2 || sz(1)~=1)
        report.identifier='seizmo:chk1dmodel:badName';
        report.message=['The .name field of MODEL ' ...
            num2str(i) ' must be a string!'];
        return;
    elseif(~islogical(model(i).ocean) || ~isscalar(model(i).ocean))
        report.identifier='seizmo:chk1dmodel:badOcean';
        report.message=['The .ocean field of MODEL ' ...
            num2str(i) ' must be TRUE or FALSE!'];
        return;
    elseif(~islogical(model(i).crust) || ~isscalar(model(i).crust))
        report.identifier='seizmo:chk1dmodel:badOcean';
        report.message=['The .crust field of MODEL ' ...
            num2str(i) ' must be TRUE or FALSE!'];
        return;
    elseif(~islogical(model(i).isotropic) || ~isscalar(model(i).isotropic))
        report.identifier='seizmo:chk1dmodel:badOcean';
        report.message=['The .isotropic field of MODEL ' ...
            num2str(i) ' must be TRUE or FALSE!'];
        return;
    elseif(~islogical(model(i).flattened) || ~isscalar(model(i).flattened))
        report.identifier='seizmo:chk1dmodel:badOcean';
        report.message=['The .flattened field of MODEL ' ...
            num2str(i) ' must be TRUE or FALSE!'];
        return;
    elseif(~isreal(model(i).refperiod) || ~isscalar(model(i).refperiod) ...
            || model(i).refperiod<=0)
        report.identifier='seizmo:chk1dmodel:badRefPeriod';
        report.message=['The .refperiod field of MODEL ' ...
            num2str(i) ' must be a scalar >0!'];
        return;
    end

    % make sure all other fields are the same size as depth
    sz=size(model(i).depth);
    for j=1:numel(fields)
        if(~isequal(sz,size(model(i).(fields{j}))))
            report.identifier='seizmo:chk1dmodel:badPropSize';
            report.message=['The .' fields{j} ' field of MODEL ' ...
                num2str(i) ' must be the same size as the .depth field!'];
            return;
        end
    end
end

end
