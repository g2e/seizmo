function [ok]=uninstall_seizmo()
%UNINSTALL_SEIZMO    Removes SEIZMO components
%
%    Usage:    ok=uninstall_seizmo
%
%    Description:
%     OK=UNINSTALL_SEIZMO will search out and remove all SEIZMO components
%     from Matlab's path and classpath.txt.  The path and classpath.txt are
%     saved with those components removed.  The path is saved to the
%     pathdef.m file that the path was loaded from at startup.  OK is TRUE
%     if the uninstall is successful.
%
%    Notes:
%     - This will uninstall ALL SEIZMO components (ie it will remove more
%       than just the SEIZMO that this function belongs to).
%
%    Examples:
%
%    See also: ABOUT_SEIZMO, SEIZMO, INSTALL_SEIZMO, SAVEPATH_SEIZMO

%     Version History:
%        Jan.  1, 2011 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan.  1, 2011 at 15:25 GMT

% todo:

% check nargin
error(nargchk(0,0,nargin));

% find classpath.txt
sjcp=which('classpath.txt');

% quick return if no classpath.txt (ie Octave)
if(isempty(sjcp)); return; end

% check that classpath is synced
% - this fails b/c of comments & arch-specific loads
%s1=javaclasspath('-static');
s2=textread(sjcp,'%s','delimiter','\n','whitespace','');
%if(~isequal(s1,s2))
%    warning('seizmo:uninstall_seizmo:restartMatlab',...
%        ['The static javaclasspath & classpath.txt are out of sync!\n' ...
%        'Please restart Matlab before uninstalling MatTauP!']);
%    return;
%end

% detect offending class path lines
yn=strfind(s2,'seizmo');
yn=~cellfun('isempty',yn);

% inform user about which lines are to be removed
disp('Removing the following lines:');
fprintf('%s\n',s2{yn});
fprintf('\nfrom:\n%s\n\n',sjcp);

% only remove if necessary
ok=false;
if(sum(yn))
    % the hard part (remove offending lines from classpath.txt
    fid=fopen(sjcp,'w');
    if(fid<0)
        warning('seizmo:uninstall_seizmo:failedToOpen',...
            ['You must have Root/Administrator privileges to edit\n' ...
            'the classpath.txt file.  To fully uninstall MatTauP\n' ...
            'you need to remove the following line(s):\n' ...
            strrep(sprintf('%s\n',s2{yn}),'\','\\') '\n' ...
            'from your Matlab''s classpath.txt located here:\n' ...
            strrep(sjcp,'\','\\') '\n\n' ...
            'This may be done by contacting your System Admin if you\n' ...
            'do not have Root/Administrator privileges.  Afterwards,\n' ...
            'please restart Matlab to complete the uninstallation!']);
    else
        fseek(fid,0,'bof');
        for i=find(~yn)'
            fprintf(fid,'%s\n',s2{i});
        end
        fclose(fid);
        fprintf(['\nSEIZMO java components will remain in running\n' ...
            'instances of Matlab.  Please restart Matlab to clear\n' ...
            'these out too.\n']);
        ok=true;
    end
else
    ok=true;
end

% find seizmo entries on path
% - this will detect all seizmo entries
p=getwords(path,':').';
yn=strfind(p,'seizmo');

% remove detections
disp('Removing the following directories from the path:');
for i=1:numel(yn)
    if(~isempty(yn{i}))
        disp(p{i});
        rmpath(p{i});
    end
end
fprintf('\n');

% save cleaned path
bad=savepath;
if(bad)
    warning('seizmo:uninstall_seizmo:failedSavingPath',...
        ['Could not save the path after uninstalling SEIZMO!\n' ...
        'Please contact your system administrator if you cannot\n' ...
        'edit and save your Matlab path!']);
end
ok=ok & ~bad;

end

function [words]=getwords(str,delimiter,collapse)
%GETWORDS    Returns a cell array of words from a string
%
%    Usage:    words=getwords('str')
%              words=getwords('str',delimiter)
%              words=getwords('str',delimiter,collapse)
%
%    Description: WORDS=GETWORDS('STR') extracts words in STR and returns
%     them separated into a cellstr array WORDS without any whitespace.
%
%     WORDS=GETWORDS('STR',DELIMITER) separates words in STR using the
%     single character DELIMITER.
%
%     WORDS=GETWORDS('STR',DELIMITER,COLLAPSE) toggles treating multiple
%     delimiters as a single delimiter.  Setting COLLAPSE to TRUE treats
%     multiple delimiters between words as a single delimiter.  Setting
%     COLLAPSE to FALSE will always return the word between a delimiter
%     pair, even if the word is '' (ie no characters).  The default is
%     TRUE.
%
%    Notes:
%
%    Examples:
%     Break up a sentence:
%      getwords('This example is pretty dumb!')
%       ans = 
%       'This'    'example'    'is'    'pretty'    'dumb!'
%
%     Turn off multi-delimiter collapsing to allow handling empty words:
%      getwords('the answer is  !',[],false)
%       ans = 
%       'the'    'answer'    'is'    [1x0 char]    '!'
%
%    See also: JOINWORDS, STRTOK, ISSPACE

%     Version History:
%        June 11, 2009 - initial version
%        Sep. 13, 2009 - minor doc update, added input check
%        Sep. 16, 2009 - add delimiter option
%        Nov. 20, 2009 - make multi-delimiter collapse optional
%        July 30, 2010 - nargchk fix
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 30, 2010 at 13:00 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% check str
if(~ischar(str) || ~isvector(str))
    error('seizmo:getwords:badInput','STR must be a char array!');
end

% check collapse
if(nargin<3 || isempty(collapse)); collapse=true; end
if(~islogical(collapse) || ~isscalar(collapse))
    error('seizmo:getwords:badInput','COLLAPSE must be a logical!');
end

% force str to row vector
str=str(:).';

% highlight word boundaries
if(nargin>1 && ~isempty(delimiter))
    % check delimiter
    if(~ischar(delimiter) || ~isscalar(delimiter))
        error('seizmo:getwords:badInput','DELIMITER must be a char!');
    end
    idx=[true str==delimiter true];
else
    idx=[true isspace(str) true];
end

% get word boundaries
if(collapse)
    idx=diff(idx);
    s=find(idx==-1);
    e=find(idx==1)-1;
else
    s=find(idx(1:end-1));
    e=find(idx(2:end))-1;
end

% number of words
nw=numel(s);

% get words
words=cell(1,nw);
for i=1:nw; words{i}=str(s(i):e(i)); end

end

