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
%    Tested on: Matlab r2007b
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
%        Oct. 17, 2008 - require new fields DIR and FILETYPE
%        Oct. 27, 2008 - LOCATION field replaces DIR field, vectorized
%                        using cellfun, global SACLAB allows skipping check
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 27, 2008 at 04:55 GMT

% todo:

% check input
report=[];
if(nargin<1)
    error('SAClab:seischk:notEnoughInputs','Not enough input arguments.');
elseif(nargin>1)
    if(~iscellstr(varargin))
        error('SAClab:seischk:badInput',...
            'Additional arguments must be strings!');
    end
end

% check SACLAB global for quick exit
global SACLAB
if(isfield(SACLAB,'SEISCHK') && isfield(SACLAB.SEISCHK,'SKIP')...
        && islogical(SACLAB.SEISCHK.SKIP)...
        && isscalar(SACLAB.SEISCHK.SKIP) && SACLAB.SEISCHK.SKIP)
    return;
end

% check data structure
if(~isstruct(data))
    report.identifier='SAClab:seischk:dataNotStruct';
    report.message='SAClab data must be a structure!';
    return;
elseif(isempty(data))
    report.identifier='SAClab:seischk:dataEmpty';
    report.message='SAClab data structure must not be empty!';
    return;
else
    defreqfields={'location' 'name' 'filetype' 'version' 'endian' 'hasdata' 'head'};
    reqfields=sort([defreqfields varargin]);
    fields=sort(fieldnames(data).');
    
    % check that all fields are present
    if(~isequal(intersect(reqfields,fields),reqfields))
        i=setxor(intersect(reqfields,fields),reqfields);
        report.identifier='SAClab:seischk:reqFieldNotFound';
        report.message=sprintf('SAClab data structure must have field ''%s''!',i{1});
        return;
    end
    
    % compile into cell arrays
    locations={data.location};
    names={data.name};
    endians={data.endian};
    versions={data.version};
    filetypes={data.filetype};
    hasdatas={data.hasdata};
    headers={data.head};
    
    % check each using cellfun
    if(any(cellfun('isempty',locations)) || ~iscellstr(locations))
        report.identifier='SAClab:seischk:nameBad';
        report.message=['SAClab struct LOCATION field must be a '...
            'nonempty string!'];
    elseif(any(cellfun('isempty',names)) || ~iscellstr(names))
        report.identifier='SAClab:seischk:dirBad';
        report.message=['SAClab struct NAME field must be a '...
            'nonempty string!'];
    elseif(any(cellfun('isempty',endians)) || ~iscellstr(endians) ||...
            ~all(strcmpi(endians,'ieee-be') | strcmpi(endians,'ieee-le')))
        report.identifier='SAClab:seischk:endianBad';
        report.message=['SAClab struct ENDIAN field must be ''ieee-le'''...
            ' or ''ieee-be''!'];
    elseif(any(cellfun('isempty',hasdatas)) ||...
            ~all(cellfun('islogical',hasdatas)))
        report.identifier='SAClab:seischk:hasdataBad';
        report.message='SAClab struct HASDATA field must be a logical!';
    elseif(~iscellstr(filetypes)...
            || any(cellfun('prodofsize',versions)~=1)...
            || any(cellfun('isempty',cellfun(@(x,y)intersect(x,y),...
            cellfun(@(x)vvseis(x),filetypes,'UniformOutput',false),...
            versions,'UniformOutput',false))))
        report.identifier='SAClab:seischk:versionBad';
        report.message=['SAClab struct FILETYPE and VERSION fields '...
            'must be valid!'];
    elseif(any(cellfun('size',headers,1)~=302)...
            || any(cellfun('size',headers,2)~=1)...
            || any(cellfun('ndims',headers)~=2))
        report.identifier='SAClab:seischk:headerBad';
        report.message='SAClab struct HEAD field must be 302x1';
    end
end

end
