function [report]=seischk(data,varargin)
%SEISCHK    Validate seislab data structure
%
%    Description: Returns an appropriate error message structure if the
%     input variable fails certain seislab data structure requirements.
%     Extra arguments must be strings corresponding to required fields in
%     the seislab data structure ('head' and 'version' are pre-included).
%     The structure contains the fields 'identifier' and 'message'.
%
%    Usage: error(seischk(data)) or warning(seischk(data))
%           error(seischk(data,'requiredfield1','requiredfield2',...))
%           warning(seischk(data,'requiredfield1','requiredfield2',...))
%
%    Examples:
%     Writing out seislab files requires the names and byte-orders of the
%     output files be set, so to include the fields 'name' and 'endian' in 
%     the check too:
%
%           error(seischk(data,'endian','name'))
%
%     Most functions require data records stored in the field 'x'.  This
%     will perform a regular check as well as assure the 'x' field exists:
%
%           error(seischk(data,'x')
%
%    See also: isseis, seishi

% check data structure
report=[];
if(~isstruct(data))
    report.identifier='seislab:seischk:dataNotStruct';
    report.message='seislab data must be a structure';
elseif(~isvector(data))
    report.identifier='seislab:seischk:dataNotVector';
    report.message='seislab data structure must be a vector';
else
    reqfields=[{'head' 'version'} varargin];
    lgc=isfield(data,reqfields);
    
    % only complain about the first field not found
    i=find(~lgc,1);
    if(~isempty(i))
        report.identifier='seislab:seischk:reqFieldNotFound';
        report.message=sprintf('seislab data structure must have field ''%s''',reqfields{i});
    end
end

end
