function [report]=lgcchk(field,lgc)
%LGCCHK    Verifies that all elements of cellstr are 'true' or 'false'
%
%    Description:  lgcchk(LGC) will generate an appropriate error message
%     structure if all the elements of cellstr array LGC are not 'true' or
%     'false'.
%
%    Usage:  error(lgcchk('fieldname',logic_field))
%
%    Examples:
%     To assure record sample spacing logical leven is set:
%       leven=glgc(data,'leven');
%       error(lgcchk('leven',leven))
%
%    See also: glgc, seischk

% check logic
report=[];
tru=strcmp(lgc,'true');
fals=strcmp(lgc,'false');
if(~all(tru | fals))
    report.identifier='SAClab:lgcchk:logicBad';
    report.message=sprintf('logical field ''%s'' needs to be set to true of false',field); 
end

end
