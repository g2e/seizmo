function [report]=seischk(data,varargin)
%SEISCHK    Validate SAClab data structure
%
%    Description: SEISCHK(DATA) returns an appropriate error message 
%     structure if the input variable fails certain SAClab data structure 
%     requirements.  The output structure contains the fields 'identifier'
%     and 'message' following Matlab error report standards.
%
%     SEISCHK(DATA,FIELD1,...,FIELDN) allows extra fields to be required in
%     addition to the default ones.  FIELD must be a string.
%
%    Notes:
%     - Current SAClab Structure Requirements
%       - Fields: dir, name, filetype, version, endian, hasdata, head
%       - All default fields must be nonempty
%       - All default fields must be valid
%     - Non-default fields are not required to be nonempty or valid
%     - See examples for non-default field uses
%
%    System requirements: Matlab 7
%
%    Header changes: NONE
%
%    Usage: error(seischk(data))
%           error(seischk(data,'requiredfield1','requiredfield2',...))
%
%    Examples:
%     Most functions require records have data stored in the field 'dep'.
%     This will perform a regular check as well as assure the field exists:
%      error(seischk(data,'dep')
%
%    See also: isseis, seisdef

%     Version History:
%        Feb. 28, 2008 - initial version
%        Mar.  2, 2008 - require nonempty data structure
%        Mar.  4, 2008 - fix error statement
%        Apr. 18, 2008 - fixed isfield to work with R14sp1
%        June 12, 2008 - doc update
%        Sep. 14, 2008 - doc update, input checks, return on first issue
%        Sep. 25, 2008 - checks versions are valid
%        Oct. 15, 2008 - data no longer required to be vector; require the
%                        name, endian, and hasdata fields by default now; 
%                        require that fields name, endian, version, hasdata
%                        and head are not empty and are valid for each
%                        record
%        Oct. 17, 2008 - require new fields 'dir' and 'filetype'
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 17, 2008 at 00:00 GMT

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
else
    defreqfields={'dir' 'name' 'filetype' 'version' 'endian' 'hasdata' 'head'};
    reqfields=[defreqfields varargin];
    
    % this works with older isfield (check one at a time)
    for i=reqfields
        if(~isfield(data,i{:}))
            report.identifier='SAClab:seischk:reqFieldNotFound';
            report.message=sprintf('SAClab data structure must have field ''%s''!',i{:});
            return;
        end
    end
    
    for i=1:numel(data)
        % check for empties
        for j=defreqfields
            if(isempty(data(i).(j{:})))
                report.identifier='SAClab:seischk:reqFieldEmpty';
                report.message=sprintf('SAClab data structure field ''%s'' must not be empty!',j{:});
                return;
            end
        end
        
        % check dir is ok
        if(~isempty(data(i).dir) && ~ischar(data(i).dir))
            report.identifier='SAClab:seischk:nameBad';
            report.message='SAClab data records must have a valid directory!';
            return;
        end
        
        % check name is ok
        if(~ischar(data(i).name))
            report.identifier='SAClab:seischk:dirBad';
            report.message='SAClab data records must have a valid name!';
            return;
        end
        
        % check version/filetype is ok
        if(~isnumeric(data(i).version)...
                || ~isequal(union(data(i).version,...
                vvseis(data(i).filetype)),vvseis(data(i).filetype)))
            report.identifier='SAClab:seischk:versionBad';
            report.message='SAClab data records must have a valid version!';
            return;
        end
        
        % check endian is ok
        if(~any(strcmp(data(i).endian,{'ieee-be' 'ieee-le'})))
            report.identifier='SAClab:seischk:endianBad';
            report.message='SAClab data records must have a valid endian!';
            return;
        end
        
        % check hasdata is ok
        if(~islogical(data(i).hasdata))
            report.identifier='SAClab:seischk:hasdataBad';
            report.message='SAClab data records must have logical HASDATA field!';
            return;
        end
        
        % check header is ok
        if(~isnumeric(data(i).head) || ~isequal(size(data(i).head),[302 1]))
            report.identifier='SAClab:seischk:headerBad';
            report.message='SAClab data records must have a valid header!';
            return;
        end
    end
end

end
