function [varargout]=h1(func)
%H1    Returns H1 line for an mfile
%
%    Usage:    h1(func)
%              str=h1(func)
%
%    Description:
%     H1(FUNC) will print out the H1 line of the function or mfile FUNC.
%     FUNC should resolve to a known function on the Matlab path or be a
%     path to an mfile.  This should match the output of the LOOKFOR
%     function (as best as I can make it at least).
%
%     STR=H1(FUNC) returns the modified H1 line (that is the line that is
%     printed and not the actual line in the mfile).
%
%    Notes:
%
%    Examples:
%     % Get the H1 line of DIR:
%     str=h1('dir')
%
%     % The H1 of UNIQUE:
%     h1('unique')
%
%    See also: LOOKFOR, HELP, CREATE_CONTENTS_FILE, CLEAN_CONTENTS_FILE

%     Version History:
%        Jan.  3, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan.  3, 2011 at 23:00 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check input
% must be either a mfile or a built-in
type2=exist(func,'file');
type5=exist(func,'builtin');
if(~type2 && ~type5)
    error('seizmo:h1:notFile',...
        'FILE: %s\nDoes not exist!',func);
end

% get path to file
if(type2)
    path=which(func);
    if(isempty(path)); path=func; end
else % type5
    path=which(func);
    path=[path(11:end-1) '.m'];
end

% function name
[name,name]=fileparts(path);

% open
fid=fopen(path,'r');

% check that happened
if(fid<0)
    error('seizmo:h1:fileNotOpenable',...
        'FILE: %s\nCannot open file!',func);
end

% find 1st comment line *with non-whitespace text* in 1st comment block
% - must be before any execution line excepting the 'function' line
% - can be before or after the 'function' line
fseek(fid,0,'bof');
h1line='';
preblock1=true; inblock1=false; prefun=true; funcont=false;
while(true)
    % get next line
    try
        line=fgetl(fid);
    catch
        error('seizmo:h1:notText',...
            'FILE: %s\nFile appears not to be text or is malformed!',func);
    end
    
    % break if eof
    if(~ischar(line)); break; end
    
    % break into words
    words=getwords(line);
    
    % skip if nothing (but catch block comment & line continuation status)
    if(isempty(words))
        if(funcont); funcont=false; end
        if(inblock1); break; end
        continue;
    end
    
    % look for comment
    if(words{1}(1)=='%' && (numel(words)>1 || numel(words{1})>1))
        % remove 1 or 2 characters or first word
        if(numel(words{1})>1)
            if(strcmp(words{1}(1:2),'%%'))
                if(numel(words{1})>2)
                    words{1}=words{1}(3:end);
                else
                    words=words(2:end);
                end
            else
                words{1}=words{1}(2:end);
            end
        else
            words=words(2:end);
        end
        
        % success!
        if(~isempty(words))
            h1line=joinwords(words);
        end
        break;
    else
        % comment block start?
        if(words{1}(1)=='%' && (inblock1 || preblock1))
            preblock1=false;
            inblock1=true;
            continue;
        elseif(~preblock1)
            % no h1!
            break;
        else
            % 4 cases:
            % 1. function line (only allow 1 of these)
            % 2. continuation of function line (allow until "line" ends)
            % 3. free continuation line (ie starting with ...) (skip all)
            % 4. execution line (fail)
            if(prefun && strcmp(words{1},'function'))
                prefun=false;
                % check for continuation
                if(~isempty(regexp(line,'^\.{3,}|[^0-9]\.{3,}','once')))
                    funcont=true;
                end
            elseif(strcmp(words{1}(1:3),'...'))
                continue;
            elseif(funcont)
                % check for continuation
                if(~isempty(regexp(line,'^\.{3,}|[^0-9]\.{3,}','once')))
                    funcont=true;
                else
                    funcont=false;
                end
            else
                break;
            end
        end
    end
end

% format the line appropriately
%XXXXX - yyyy yyy yyyyyyyy
maxchar=30;
wlen=numel(name);
if(isempty(h1line))
    h1line=['<a href="matlab:help ' name '">' name '</a>' ...
        char(32*ones(1,min(1,maxchar-wlen))) '- '];
else
    % remove leading function name if there
    if(numel(words{1})>=wlen)
        if(strcmpi(words{1}(1:wlen),name))
            if(numel(words{1})==wlen)
                words=words(2:end);
            else
                words{1}=words{1}(wlen+1:end);
            end
        end
    end
    % remove leading dash if there
    if(strcmpi(words{1}(1),'-'))
        if(numel(words{1})>1)
            words{1}=words{1}(2:end);
        else
            words=words(2:end);
        end
    end
    h1line=['<a href="matlab:help ' name '">' name '</a>' ...
        char(32*ones(1,max(1,1+maxchar-wlen))) '- ' joinwords(words)];
end

% close
fclose(fid);

% print line if no output
if(nargout)
    % output
    varargout{1}=h1line;
else
    % print
    fprintf('%s\n',h1line);
end

end
