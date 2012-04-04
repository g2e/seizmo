function [report]=chk1dmodel(model)
%CHK1DMODEL    Validate if a struct is as defined by PREM/AK135/IASP91
%
%    Usage:    msg=chk1dmodel(model)
%
%    Description:
%     MSG=CHK1DMODEL(MODEL) checks that MODEL is a struct as output from
%     PREM/AK135/IASP91.  See those functions for details on the struct
%     layout.  MSG is an error structure if a problem is found (otherwise
%     it is empty).
%
%    Notes:
%
%    Examples:
%     % Check that output from PREM is compatible:
%     mod=prem;
%     error(chk1dmodel(mod));
%
%    See also: PREM, AK135, IASP91

%     Version History:
%        May  24, 2010 - initial version
%        May  27, 2010 - now allows arrays of models
%        Aug. 10, 2010 - require monotoniticity of .depth
%        Aug. 17, 2010 - require nonempty, fix error message
%        Sep. 19, 2010 - more .depth checks
%        Jan. 25, 2011 - require .crust is true if .ocean is true, fix some
%                        warning ids
%        Apr.  3, 2012 - minor doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr.  3, 2012 at 15:05 GMT

% todo:
% - models should have named discontinuities field
%   ('moho' '410' '660' 'cmb' 'icb')

% check nargin
report=[];
error(nargchk(1,1,nargin));

% check model has the required fields
reqfields={'name' 'ocean' 'crust' 'isotropic' 'refperiod' 'flattened' ...
    'depth'};
if(~isstruct(model) || ~all(ismember(reqfields,fieldnames(model))) ...
        || isempty(model))
    report.identifier='seizmo:chk1dmodel:dataNotStruct';
    report.message=[...
        sprintf('MODEL must be a non-empty struct with fields:\n') ...
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
        report.identifier='seizmo:chk1dmodel:badCrust';
        report.message=['The .crust field of MODEL ' ...
            num2str(i) ' must be TRUE or FALSE!'];
        return;
    elseif(model(i).ocean && ~model(i).crust)
        report.identifier='seizmo:chk1dmodel:badOceanCrust';
        report.message=['The .ocean field of MODEL ' ...
            num2str(i) ' is TRUE: .crust must be TRUE too!'];
        return;
    elseif(~islogical(model(i).isotropic) || ~isscalar(model(i).isotropic))
        report.identifier='seizmo:chk1dmodel:badIso';
        report.message=['The .isotropic field of MODEL ' ...
            num2str(i) ' must be TRUE or FALSE!'];
        return;
    elseif(~islogical(model(i).flattened) || ~isscalar(model(i).flattened))
        report.identifier='seizmo:chk1dmodel:badFlat';
        report.message=['The .flattened field of MODEL ' ...
            num2str(i) ' must be TRUE or FALSE!'];
        return;
    elseif(~isreal(model(i).refperiod) || ~isscalar(model(i).refperiod) ...
            || model(i).refperiod<=0)
        report.identifier='seizmo:chk1dmodel:badRefPeriod';
        report.message=['The .refperiod field of MODEL ' ...
            num2str(i) ' must be a scalar >0!'];
        return;
    elseif(isempty(model(i).depth) || ~isreal(model(i).depth) ...
            || any(isnan(model(i).depth)) ...
            || ~isvector(model(i).depth) || size(model(i).depth,2)~=1)
        report.identifier='seizmo:chk1dmodel:badDepth';
        report.message=['The .depth field of MODEL ' ...
            num2str(i) ' must be a non-empty real-valued column vector!'];
        return;
    elseif(any(diff(model(i).depth)<0))
        report.identifier='seizmo:chk1dmodel:badDepth';
        report.message=['The .depth field of MODEL ' ...
            num2str(i) ' must be monotonically non-decreasing!'];
        return;
    elseif(any(histc(model(i).depth,...
                    model(i).depth([find(diff(model(i).depth));end]))>3))
        report.identifier='seizmo:chk1dmodel:badDepth';
        report.message=['The .depth field of MODEL ' ...
            num2str(i) ' must not have values repeated 3+ times!'];
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
