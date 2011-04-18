function []=clean_contents_file(file,skip,dbl)
%CLEAN_CONTENTS_FILE    Cleans up formatting of a Contents.m file
%
%    Usage:    clean_contents_file(file,[hlines tlines])
%              clean_contents_file(file,[hlines tlines],dbl)
%
%    Description:
%     CLEAN_CONTENTS_FILE(FILE,[HLINES TLINES]) edits the mfile given by
%     FILE that conforms to typical Matlab Help system.  The editing will
%     add HREF links around the first word of each line (should be the
%     function name) and space out the description that follows to be
%     indented consistently with the rest of the function descriptions.
%     The number of header and trailing lines to skip may be set with the
%     second arg [HLINES TLINES].
%
%     CLEAN_CONTENTS_FILE(FILE,[HLINES TLINES],DBL) passes the flag DBL
%     (must be TRUE or FALSE -- default is FALSE) indicating if the
%     Contents file is double indented.  The default (FALSE) is single
%     indent, which is the usual case.  Double indenting is for listing
%     shortcut function names by their longer counterparts.
%
%    Notes:
%
%    Examples:
%     % Copy all the H1 lines of the functions in a directory into an
%     % m-file, then call CLEAN_CONTENTS_FILE on that m-file to make
%     % it look slick.  You can continue to add new H1 lines to that file
%     % as new functions are added.  Just call CLEAN_CONTENTS_FILE on the
%     % m-file again to update it.
%     ... make a m-file of H1-lines ...
%     clean_contents_file(my_file)
%     ... add some more H1-lines to the same file ...
%     clean_contents_file(my_file)
%
%    See also: HELP

%     Version History:
%        Jan.  2, 2011 - initial version
%        Jan. 23, 2011 - allow files on path, fix path bug
%        Jan. 26, 2011 - allow help function to not match string, 3rd arg
%                        allows for double indent content files
%        Apr. 16, 2011 - better handling of dbl
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Apr. 16, 2011 at 23:00 GMT

% todo:

% check nargin
error(nargchk(0,3,nargin));

% default inputs
if(nargin<1); file=[]; end
if(nargin<2); skip=[]; end
if(nargin<3); dbl=false; end

% check double indent flag
if(~isscalar(dbl) || ~islogical(dbl))
    error('seizmo:clean_contents:file:badDoubleIndentFlag',...
        'DBL must be TRUE or FALSE!');
end

% filespec
global SEIZMO
SEIZMO.ONEFILELIST.FILTERSPEC={'*.m;*.M' 'M Files (*.m,*.M)'};

% get files
list=onefilelist(file);
nfiles=numel(list);

% if nothing listed, search on path
if(~nfiles)
    list=which(file);
    
    % error if none
    if(isempty(list))
        error('seizmo:clean_contents_file:noFiles',...
            'Unknown File: %',file);
    end
    
    % now replace with xdir output
    list=xdir(list);
    nfiles=1;
end

% check skip
if(isempty(skip))
    skip=zeros(nfiles,2);
elseif(any(size(skip,1)==[1 nfiles]) && any(size(skip,2)==[1 2]) ...
        && all(skip==fix(skip)) && all(skip>=0))
    if(size(skip,2)==1); skip(:,2)=0; end
    if(size(skip,1)==1); skip=skip(ones(nfiles,1),:); end
else
    error('seizmo:clean_contents_file:badInput',...
        'SKIP must be formatted as [HLINES TLINES]!');
end

% loop over files
for a=1:nfiles
    % make path/file
    p2f=fullfile(list(a).path,list(a).name);
    
    % read in file
    lines=getwords(readtxt(p2f),sprintf('\n'));
    nlines=numel(lines);
    
    % figure out maximum characters for indent
    maxchar=0; dblmax=0;
    for b=1:nlines
        % look only if in range
        if(b>skip(a,1) && b<nlines+1-skip(a,2))
            % split into words
            words=getwords(lines{b});
            
            % skip non-comment (%)
            if(strcmp(words{1}(1),'%'))
                % strip comment char(s)
                while(~isempty(words) && strcmp(words{1}(1),'%'))
                    % strip first character or first word?
                    if(numel(words{1})>1)
                        words{1}=words{1}(2:end);
                    else
                        words=words(2:end);
                    end
                end
                
                % skip if nothing left
                if(~numel(words)); continue; end
                
                % is it a href line or a raw one?
                if(strcmp(words{1}(1),'<'))
                    % href line looks like follows:
                    % <a href="matlab:help xxxx">zzzz</a>  - yyy yyyyy
                    %
                    % We want zzzz length, which is length of the 3rd
                    % "word" divided by 2 minus 3 if xxxx=zzzz.
                    %if((numel(words{3})/2-3)>maxchar)
                    %    maxchar=numel(words{3})/2-3;
                    %end
                    % This works when xxxx~=zzzz
                    tmpwords=getwords(words{3},'>');
                    if((numel(tmpwords{2})-3)>maxchar)
                        maxchar=numel(tmpwords{2})-3;
                    end
                    
                    % handle double indent
                    % getting yyy from above
                    if(dbl && numel(words{5})>dblmax)
                        dblmax=numel(words{5});
                    end
                else
                    % update maxchar if needed
                    if(numel(words{1})>maxchar)
                        maxchar=numel(words{1});
                    end
                    
                    % handle double indent
                    % XXXX - YYYY - zzz zzzzz zz
                    if(dbl && numel(words{3})>dblmax)
                        dblmax=numel(words{3});
                    end
                end
            end
        end
    end
    
    % open for editing
    fid=fopen(p2f,'w');
    if(fid<0)
        error('seizmo:clean_contents_file:cannotOpenFile',...
            'File: %s\nNot Openable!',list(a).name);
    end
    fseek(fid,0,'bof');
    
    % edit loop
    for b=1:nlines
        % edit only if in range
        if(b>skip(a,1) && b<nlines+1-skip(a,2))
            % split into words
            words=getwords(lines{b});
            
            % skip non-comment (%)
            if(strcmp(words{1}(1),'%'))
                % remove comment char(s)
                while(~isempty(words) && strcmp(words{1}(1),'%'))
                    % remove first character or first word
                    if(numel(words{1})>1)
                        words{1}=words{1}(2:end);
                    else
                        words=words(2:end);
                    end
                end
                
                % write comment line if nothing left
                if(~numel(words))
                    fprintf(fid,'%%\n');
                    continue;
                end
                
                % is it a href line or a raw one?
                if(strcmp(words{1}(1),'<'))
                    % double indented?
                    if(dbl) % double indent
                        % href line looks like follows:
                        % <a href="matlab:help xxx">zzz</a> - XXX - yyy yyy
                        % wlen=numel(words{3})/2-3; % old way (xxxx==zzzz)
                        tmpwords=getwords(words{3},'>');
                        wlen=numel(tmpwords{2})-3;
                        dbllen=numel(words{5});
                        lines{b}=['%' joinwords(words(1:3)) ...
                            char(32*ones(1,1+maxchar-wlen)) ...
                            joinwords(words(4:5)) ...
                            char(32*ones(1,1+dblmax-dbllen)) ...
                            joinwords(words(6:end))];
                    else % single indent
                        % href line looks like follows:
                        % <a href="matlab:help xxxx">zzzz</a>  - yyy yyyyy
                        % wlen=numel(words{3})/2-3; % old way (xxxx==zzzz)
                        tmpwords=getwords(words{3},'>');
                        wlen=numel(tmpwords{2})-3;
                        lines{b}=['%' joinwords(words(1:3)) ...
                            char(32*ones(1,1+maxchar-wlen)) ...
                            joinwords(words(4:end))];
                    end
                else
                    % strip separator if there is one
                    if(strcmp(words{2},'-')); words(2)=[]; end
                    
                    % double indent?
                    if(dbl) % double indent
                        % make new line
                        wlen=numel(words{1});
                        dbllen=numel(words{2});
                        lines{b}=['%<a href="matlab:help ' words{1} ...
                            '">' words{1} '</a>' ...
                            char(32*ones(1,1+maxchar-wlen)) '- ' ...
                            words(2) ...
                            char(32*ones(1,1+dblmax-dbllen)) ...
                            joinwords(words(3:end))];
                    else % single indent
                        % make new line
                        wlen=numel(words{1});
                        lines{b}=['%<a href="matlab:help ' words{1} ...
                            '">' words{1} '</a>' ...
                            char(32*ones(1,1+maxchar-wlen)) '- ' ...
                            joinwords(words(2:end))];
                    end
                end
            end
        end
        
        % write to file
        fprintf(fid,'%s\n',lines{b});
    end
    
    % close file
    fclose(fid);
end

end
