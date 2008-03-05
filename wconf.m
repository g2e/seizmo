function []=wconf(conf)
%WCONF    Write configuration structure to file
%
%    Description:  Writes a structure to an ascii file.  Requires the field
%     CONF.OUTPUTCONFFILE to be defined as a string representing the save
%     path/file.name.
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
%    Usage: wconf(CONF)
%
%    Examples:
%
%    See also: rconf, makespecs

% get info on structure
[fields,rows,cols,rspecs,wspecs]=makespecs(conf);

% open output file for writing
fid=fopen(conf.OUTPUTCONFFILE,'w');
if(fid<0)
    error('SAClab:wconf:badFID',...
        'Could not open file %s',conf_file_path); 
end
fseek(fid,0,'bof');

% loop through structure elements we have info on
for i=1:length(fields)
    % empty element shortcut
    if(rows(i)==0)
        fprintf(fid,'%s %d',fields{i},0);
    % special handling of cell/cellstr
    elseif(strncmp('%s',rspecs{i},2))
        % first line (without values yet)
        fprintf(fid,'%s %d %s',fields{i},rows(i),rspecs{i});
        % fprintf doesn't support cell input so loop through every element
        for j=1:rows(i)
            for k=1:cols(i)
                % tabbing subsequent rows to make easy reading
                if(k==1 && j>1); fprintf(fid,'\t\t\t'); end
                fprintf(fid,' %s',conf.(fields{i}){j,k});
            end
            % end of line
            fprintf(fid,'\n');
        end
    % char/numeric
    else
        % first line
        fprintf(fid,['%s %d %s' wspecs{i} '\n'],fields{i},rows(i),rspecs{i},conf.(fields{i})(1,:));
        % subsequent lines (rows) are indented
        for j=2:rows(i)
            fprintf(fid,['\t\t\t' wspecs{i} '\n'],conf.(fields{i})(j,:));
        end
    end
end

% close file
fclose(fid);

end
