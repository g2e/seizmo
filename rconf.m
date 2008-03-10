function [conf]=rconf(conf_file_path,conf)
%RCONF    Read configuration file and assemble in a structure
%
%    Description:  Reads in a configuration file written in the format used
%     by the wconf function.  Takes a string that is the path to the
%     configuration file and optionally a structure (for updating).
%
%     Format of ascii file:
%       fieldname1 nrows spec val11 ... val1n
%                               .   ...   .
%                               .   ...   .
%                               .   ...   .
%                             valm1 ... valmn % comments ...
%       % comments ...
%       fieldname2 ...
%       ...
%
%    Notes:
%     - supports commenting with a '%' at the start of a line or at the end
%       of a value array (no commenting individual rows of a value array)
%     - does not support cell/cellstr arrays that have strings with spaces
%       but does support char arrays with spaces
%     - does not support substructures
%     - everything is whitespace delimited
%
%    Usage: CONF=rconf('path/to/file.conf')
%           CONF=rconf('path/to/file.conf',CONF)
%
%    Examples:
%
%    See also: wconf, makespecs

% open file for reading at its beginning
fid=fopen(conf_file_path);
if(fid<0)
    error('SAClab:rconf:badFID',...
        'Could not open config file %s',conf_file_path); 
end
fseek(fid,0,'bof');

% loop through field/value(s) pairs
while 1
    % read in field
    field=textscan(fid,'%s',1);
    
    % check if empty (no more input = eof)
    if(isempty(field{:})); break; end
    
    % uncell the field string
    field=field{:}{:};
    
    % make sure field is ok
    try
        conf.(field)=[];
        conf=rmfield(conf,field);
    catch
        % bad field name - skip rest of line
        fgetl(fid); continue;
    end
    
    % read in number of rows for field value array
    n=textscan(fid,'%d',1); n=n{:};
    
    % check if field value array is empty (nrows==0) - skip rest of line if so
    if(isempty(n)); fgetl(fid); continue; end
    if(n==0); conf.(field)=[]; fgetl(fid); continue; end
    
    % read in field's value format (indicates number of columns too)
    spec=textscan(fid,'%s',1); spec=spec{:}{:};
    
    % read in field's value array
    try
        values=textscan(fid,spec,n);
    catch
        warning(lasterr);
        fgetl(fid); continue;
    end
    
    % uncell char arrays
    if(length(spec)>2 ...
            && isscalar(regexp(spec,'%')) ...
            && strcmp(spec([1 end]),'%s'));
        values=values{:}; 
    end
    
    % assign field's value array to field
    conf.(field)=[values{:}];
end

% include input configuration file path/name
conf.INPUTCONFFILE=conf_file_path;

% close configuration file
fclose(fid);

end
