function [report]=seischk(data,varargin)
%SEISCHK    Validate SAClab data structure
%
%    Description: SEISCHK(DATA) returns an appropriate error message 
%     structure if the input variable fails certain SAClab data structure 
%     requirements (must be a nonempty 1D structure that has 'head' and 
%     'version' fields).  The output structure contains the fields 
%     'identifier' and 'message' following Matlab error report standards.
%
%     SEISCHK(DATA,FIELD1,...,FIELDN) allows extra fields to be required in
%     addition to the 'head' and 'version' fields.  FIELD must be a string.
%
%    Notes:
%     - see examples for extra field uses
%
%    System requirements: Matlab 7
%
%    Input/Output requirements: first arg can be anything, additional args
%     must be strings
%
%    Header changes: N/A
%
%    Usage: error(seischk(data))
%           warning(seischk(data))
%           error(seischk(data,'requiredfield1','requiredfield2',...))
%           warning(seischk(data,'requiredfield1','requiredfield2',...))
%
%    Examples:
%     Writing out SAClab files requires the names and byte-orders of the
%     output files be set, so to include the fields 'name' and 'endian' in 
%     the check:
%
%           error(seischk(data,'endian','name'))
%
%     Most functions require data records stored in the field 'x'.  This
%     will perform a regular check as well as assure the 'x' field exists:
%
%           error(seischk(data,'x')
%
%    See also: isseis, seisdef

%     Version History:
%        Feb. 28, 2008 - initial version
%        Mar.  2, 2008 - require nonempty data structure
%        Mar.  4, 2008 - fix error statement
%        Apr. 18, 2008 - fixed isfield to work with R14sp1
%        June 12, 2008 - doc update
%        Sep. 14, 2008 - doc update, input checks, return on first issue
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 14, 2008 at 17:45 GMT

% todo:

% check input
if(nargin<1)
    error('SAClab:seischk:notEnoughInputs','Not enough input arguments.');
elseif(nargin>1)
    for i=1:nargin-1
        if(~ischar(varargin{i}))
            error('SAClab:seischk:badInput',...
                'Additional argument FIELD%d must be a string!',i);
        end
    end
end

% check data structure
report=[];
if(~isstruct(data))
    report.identifier='SAClab:seischk:dataNotStruct';
    report.message='SAClab data must be a structure!';
    return;
elseif(isempty(data))
    report.identifier='SAClab:seischk:dataEmpty';
    report.message='SAClab data structure must not be empty!';
    return;
elseif(~isvector(data))
    report.identifier='SAClab:seischk:dataNotVector';
    report.message='SAClab data structure must be a vector!';
    return;
else
    reqfields=[{'head' 'version'} varargin];
    
    % this works with older isfield (check one at a time)
    for i=reqfields
        if(~isfield(data,i{:}))
            report.identifier='SAClab:seischk:reqFieldNotFound';
            report.message=sprintf('SAClab data structure must have field ''%s''!',i{:});
            return;
        end
    end
end

end
