function [varargout]=install_seizmo()
%INSTALL_SEIZMO    Installs the SEIZMO Toolbox in Matlab or Octave
%
%    Usage:    install_seizmo
%
%    Description:
%     INSTALL_SEIZMO installs the SEIZMO toolbox in Matlab or Octave.  This
%     mainly involves editing the Matlab/Octave path so you get access to
%     all the functions that make SEIZMO work.  This will uninstall
%     previous SEIZMO installs (thus this is a way to "update").  You can
%     do an uninstallation yourself by calling UNINSTALL_SEIZMO (its hidden
%     in the 'uninstall' directory).
%
%    Notes:
%     - INSTALL_SEIZMO must be called from *WITHIN* a running session of
%       Matlab or Octave.  This is NOT a shell script or the Windows
%       equivalent, so do NOT "run" install_seizmo.m from your OS shell.
%       Start up Matlab/Octave and in the command window type
%       "install_seizmo" without the quotes and press enter.  Read what it
%       says and with any luck you will be ready to go after it finishes!
%
%     - Websites:
%        SEIZMO    - http://epsc.wustl.edu/~ggeuler/codes/m/seizmo/
%        MatTauP   - http://www.ess.washington.edu/SEIS/FMI/matTaup.htm
%        M_Map     - http://www.eos.ubc.ca/~rich/map.html
%        njTBX     - http://sourceforge.net/apps/trac/njtbx
%        GSHHS     - http://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html
%        GlobalCMT - http://www.globalcmt.org/
%
%    Examples:
%     % Amazingly, every step of installing SEIZMO can be done *WITHIN*
%     % Matlab/Octave.  All you need is an internet connection so you
%     % can grab the SEIZMO package, extract the contents and install:
%     file=gunzip('https://github.com/g2e/seizmo/tarball/master');
%     untar(file{:});
%     cd seizmo/
%     install_seizmo
%
%    See also: ABOUT_SEIZMO, SEIZMO, UNINSTALL_SEIZMO, WEBINSTALL_NJTBX,
%              WEBINSTALL_MMAP, WEBINSTALL_GSHHS

%     Version History:
%        Dec. 30, 2010 - initial version
%        Jan.  1, 2011 - added msgs, detects more issues
%        Jan.  4, 2011 - improved docs (examples needs more work)
%        Jan. 14, 2011 - example improved to get current release
%        Jan. 19, 2011 - updated to fixed examples
%        Apr.  6, 2011 - include verLessThan for pre-7.4 matlab
%        June 16, 2011 - doc update
%        June 24, 2011 - fix for octave warning in verLessThan, better
%                        uninstall output
%        Feb. 15, 2012 - use version_compare, clean up messages, drop
%                        separate installers except for external
%                        components, update cmt db, flip savepath logic,
%                        only use javaaddpath or edit classpath as needed
%        Feb. 16, 2012 - export_fig is externally managed
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 16, 2012 at 15:25 GMT

% todo:

% check nargin
error(nargchk(0,1,nargin));

% check application & version
ok=true;
disp('##################################################################');
disp('################## STARTING SEIZMO INSTALLATION ##################');
disp('##################################################################');
disp('Checking what application we are installing SEIZMO in...')
[application,version]=getapplication;
disp(['Application:  ' application]);
disp(['Version    :  ' version]);
switch lower(application)
    case 'matlab'
        if(version_compare(version,'7.1')<=0)
            warning('seizmo:install_seizmo:versionBad',...
                ['Matlab version too old for SEIZMO!\n' ...
                'Full function of toolbox is unlikely.']);
        end
        if(~license('checkout','signal_toolbox'))
            warning('seizmo:install_seizmo:noSigProcTbx',...
                ['Your Matlab does not have the Signal Processing\n' ...
                'Toolbox installed!  SEIZMO needs the Signal\n' ...
                'Processing Toolbox to be fully functional!']);
        elseif(~license('checkout','statistics_toolbox'))
            warning('seizmo:install_seizmo:noStatsTbx',...
                ['Your Matlab does not have the Statistics\n' ...
                'Toolbox installed!  SEIZMO needs the Statistics\n' ...
                'Toolbox to be fully functional!']);
        end
    case 'octave'
        warning('seizmo:install_seizmo:octaveIssues',...
            'Octave compatibility for SEIZMO is a work in progress!');
    otherwise
        warning('seizmo:install_seizmo:noClueWhatIsRunning',...
            'Installing SEIZMO on an UNKNOWN application!');
end

% remove old seizmo installation(s)
info=ver('seizmo');
while(~isempty(info))
    disp(['Uninstalling previous SEIZMO version: ' info.Version]);
    ok=uninstall_seizmo;
    if(~ok)
        warning('seizmo:install_seizmo:failedUninstall',...
            'Uninstalling previous SEIZMO failed!');
        if(nargout); varargout{1}=ok; end
        return;
    end
    info=ver('seizmo');
end

% where am i?
path=fileparts(mfilename('fullpath'));
disp(['SEIZMO install path:  ' path]);

% install new seizmo
fs=filesep;
addpath(path,...
    [path fs 'lowlevel'],...
    [path fs 'uninstall'],...
    [path fs 'behavior'],...
    [path fs 'toc'],...
    [path fs 'rw'],...
    [path fs 'hdr'],...
    [path fs 'sz'],...
    [path fs 'misc'],...
    [path fs 'time'],...
    [path fs 'position'],...
    [path fs 'audio'],...
    [path fs 'cmap'],...
    [path fs 'cmb'],...
    [path fs 'cmt'],...
    [path fs 'decon'],...
    [path fs 'event'],...
    [path fs 'filtering'],...
    [path fs 'fixes'],...
    [path fs 'fk'],...
    [path fs 'ftran'],...
    [path fs 'gui'],...
    [path fs 'invert'],...
    [path fs 'mapping'],...
    [path fs 'models'],...
    [path fs 'multi'],...
    [path fs 'noise'],...
    [path fs 'pick'],...
    [path fs 'plotting'],...
    [path fs 'resampling'],...
    [path fs 'response'],...
    [path fs 'shortnames'],...
    [path fs 'solo'],...
    [path fs 'sphpoly'],...
    [path fs 'synth'],...
    [path fs 'tomo'],...
    [path fs 'topo'],...
    [path fs 'tpw'],...
    [path fs 'ttcorrect'],...
    [path fs 'win'],...
    [path fs 'ww3'],...
    [path fs 'xcalign'],...
    [path fs 'mattaup']);
ok=ok & ~savepath;
if(~ok)
    warning('seizmo:install_seizmo:noPermission',...
        'Could not edit path!');
end

% check that classpath exists (Octave fails here)
mattaupjar=fullfile(path,'mattaup','lib','matTaup.jar');
sjcp=which('classpath.txt');
if(isempty(sjcp))
    %warning('seizmo:webinstall_njtbx:noJavaClassPath',...
    %    'Octave has no classpath.txt to save .jar files!');
    
    % no classpath.txt so add to dynamic path
    if(~ismember(mattaupjar,javaclasspath)); javaaddpath(mattaupjar); end
else % install matTaup.jar to classpath
    % read classpath.txt
    s2=textread(sjcp,'%s','delimiter','\n','whitespace','');
    
    % detect offending classpath.txt lines
    yn=~cellfun('isempty',strfind(s2,mattaupjar));
    
    % only add if not there already
    if(~sum(yn))
        fid=fopen(sjcp,'a+');
        if(fid<0)
            warning('seizmo:webinstall_njtbx:noWriteClasspath',...
                ['Cannot edit classpath.txt! Adding ' ...
                'matTaup.jar to dynamic java class path!']);
            if(~ismember(mattaupjar,javaclasspath))
                javaaddpath(mattaupjar);
            end
        else
            fseek(fid,0,'eof');
            fprintf(fid,'%s\n',mattaupjar);
            fclose(fid);
        end
    end
end

% ask to install external components
% - eventually: taup
reply=input('Install njTBX (30MB)? Y/N [Y]: ','s');
if(isempty(reply) || strncmpi(reply,'y',1))
    ok=ok & webinstall_njtbx;
end
reply=input('Install M_Map (<1MB)? Y/N [Y]: ','s');
if(isempty(reply) || strncmpi(reply,'y',1))
    ok=ok & webinstall_mmap;
    reply=input('Install GSHHS (120MB!)? Y/N [Y]: ','s');
    if(isempty(reply) || strncmpi(reply,'y',1))
        ok=ok & webinstall_gshhs;
    end
end
reply=input('Install export_fig (25KB)? Y/N [Y]: ','s');
if(isempty(reply) || strncmpi(reply,'y',1))
    ok=ok & webinstall_exportfig;
end

% update dbs
% - eventually: iris sacpzdb
reply=input('Update GlobalCMT database? Y/N [Y]: ','s');
if(isempty(reply) || strncmpi(reply,'y',1)); globalcmt_update; end

% download models/features
% - eventually: 3D mantle models & features in mapping

disp('DON''T FORGET TO RESTART!');
disp('##################################################################');
disp('################## FINISHED SEIZMO INSTALLATION ##################');
disp('##################################################################');
fprintf('\n\n\n');

% display helpful info
about_seizmo;
help seizmo;

% output
if(nargout); varargout{1}=ok; end

end

function [varargout]=getapplication()
%GETAPPLICATION    Returns application running this script and its version
%
%    Usage:    [application,version]=getapplication()
%
%    Description:
%     [APPLICATION,VERSION]=GETAPPLICATION() will determine and return the
%     name and version of the application running this script (obviously
%     only if the application can run this script in the first place).
%     Both APPLICATION and VERSION are strings.
%
%    Notes:
%     - returns 'UNKNOWN' if it cannot figure out the application
%
%    Examples:
%     % Matlab and Octave still behave quite differently for a number of
%     % different functions so it is best in some cases to use different
%     % function calls depending on which we are running:
%     [app,ver]=getapplication;
%     if(strcmp(app,'MATLAB'))
%       % do something via matlab routines
%     else
%       % do something via octave routines
%     end
%
%    See also: NATIVEBYTEORDER, VER

%     Version History:
%        Nov. 13, 2008 - initial version
%        Mar.  3, 2009 - minor doc cleaning
%        Apr. 23, 2009 - move usage up
%        Sep.  8, 2009 - minor doc update
%        Feb. 15, 2012 - minor doc update, cache values using persistent
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 15, 2012 at 19:55 GMT

% todo:

% cache as persistent variables
persistent application version
if(~isempty(application) && ~isempty(version))
    varargout={application version};
    return;
end

% checking for Matlab will throw an error in Octave
try
    % first check if we are in Matlab
    a=ver('matlab');
    
    % we are in Matlab
    application=a.Name;
    version=a.Version;
catch
    % check if we are in Octave
    if(exist('OCTAVE_VERSION','builtin')==5)
        application='OCTAVE';
        version=OCTAVE_VERSION;
    % ok I have no clue what is running
    else
        application='UNKNOWN';
        version='UNKNOWN';
    end
end

% output
varargout={application version};

end

function [cmp]=version_compare(ver1,ver2)
%VERSION_COMPARE    Compares versions strings given in XX.XX.XX.... format
%
%    Usage:    cmp=version_compare(ver1,ver2)
%
%    Description:
%     CMP=VERSION_COMPARE(VER1,VER2) compares the versions in VER1 & VER2
%     returning CMP=1 if VER1>VER2, CMP=-1 if VER1<VER2, or CMP=0 if
%     VER1=VER2.  VER1 & VER2 must be strings.  Versions are split based on
%     the '.' delimiter (typically versions are formatted as major.minor
%     and so on).  Alphabetical characters have higher values than numbers
%     so that 'a'>'9'.
%
%    Notes:
%
%    Examples:
%     % Compare some simple versions:
%     cmp=version_compare('1','3')
%
%     % A case from wgrib2:
%     cmp=version_compare('v0.1.5f','v0.1.9.4')
%
%    See also: VERLESSTHAN, VER, PARSE_ALPHANUMERIC

%     Version History:
%        Feb. 13, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 13, 2012 at 13:00 GMT

% todo:

% check nargin
error(nargchk(2,2,nargin));

% require strings
if(~(ischar(ver1) || iscellstr(ver1)) ...
        || ~(ischar(ver2) || iscellstr(ver2)))
    error('seizmo:version_compare:badInput',...
        'VER must be a string!');
end

% now convert to cell strings
ver1=cellstr(ver1);
ver2=cellstr(ver2);

% expand scalars
if(isscalar(ver1)); ver1=ver1(ones(numel(ver2),1),1); end
if(isscalar(ver2)); ver2=ver2(ones(numel(ver1),1),1); end
ncmp=numel(ver1);

% loop over each row in ver1/ver2
cmp=nan(ncmp,1);
for i=1:ncmp
    % require single row
    if(size(ver1{i},1)>1 || size(ver2{i},1)>1)
        error('seizmo:version_compare:badInput',...
            'VER string must be a row vector!');
    end
    
    % get major.minor.revision...
    f1=getwords(ver1{i},'.');
    f2=getwords(ver2{i},'.');
    
    % loop over comparible fields until there is a difference
    nf1=numel(f1); nf2=numel(f2);
    for j=1:min(nf1,nf2)
        % try parsing as numbers
        d1=str2double(f1{j});
        d2=str2double(f2{j});
        
        if(~isnan(d1) && ~isnan(d2))
            % both are numbers
            cmp(i)=sign(d1-d2);
        else
            % parse as alphanumeric
            [an1,isnum1]=parse_alphanumeric(f1{j});
            [an2,isnum2]=parse_alphanumeric(f2{j});
            nbit1=numel(isnum1);
            nbit2=numel(isnum2);
            
            % loop over parsed bits
            for k=1:min(nbit1,nbit2)
                % n vs n
                if(isnum1(k) && isnum2(k))
                    cmp(i)=sign(an1{k}-an2{k});
                % a vs a
                elseif(~isnum1(k) && ~isnum2(k))
                    last=max(numel(an1{k}),numel(an2{k}));
                    for l=1:last
                        cmp(i)=sign(an1{k}(l)-an2{k}(l));
                        if(cmp(i)); break; end
                    end
                    % no diff so use number of alpha
                    if(~cmp(i))
                        cmp(i)=sign(numel(an1{k})-numel(an2{k}));
                    end
                % n vs a
                elseif(isnum1(k) && ~isnum2(k))
                    % alpha wins (0-9, a-z)
                    cmp(i)=-1;
                % a vs n
                elseif(~isnum1(k) && isnum2(k))
                    % alpha wins (0-9, a-z)
                    cmp(i)=1;
                end
                if(cmp(i)); break; end
            end
        end
        
        % continue until 1 or -1
        if(cmp(i)); break; end
    end
    
    % failed to find a difference in comparible fields
    % so just use the number of fields
    if(isnan(cmp(i)) || ~cmp(i)); cmp(i)=sign(nf1-nf2); end
end

end

function [an,isnum]=parse_alphanumeric(str)
%PARSE_ALPHANUMERIC    Split alphanumeric string into words & numbers
%
%    Usage:    [an,isnum]=parse_alphanumeric(str)
%
%    Description:
%     [AN,ISNUM]=PARSE_ALPHANUMERIC(STR) parses out alphabet and digit
%     sequences from character string STR as a cell array of "words" and
%     numbers (converted to double) in AN.  ISNUM indicates the elements in
%     AN that are numeric.
%
%    Notes:
%
%    Examples:
%     % Split a date string:
%     parse_alphanumeric('2000may03')
%
%    See also: ISSTRPROP, GETWORDS, READTXT, VERSION_COMPARE

%     Version History:
%        Feb. 13, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 13, 2012 at 13:00 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% require character input
if(~ischar(str))
    error('seizmo:parse_alphanumeric:badInput',...
        'STR must be a character string!');
end

% identify alphanumeric bits
num=isstrprop(str,'digit');
alf=isstrprop(str,'alpha');

% loop over bits getting words & numbers
cnt=1; in=false; idx=nan(0,2); isnum=false(0,1);
for i=1:numel(str)
    if(num(i))
        % number
        if(in(cnt) && isnum(cnt))
            % continuing number
            continue;
        elseif(in(cnt) && ~isnum(cnt))
            % terminate alpha, begin number
            idx(cnt,2)=i-1;
            cnt=cnt+1;
            idx(cnt,1)=i;
            isnum(cnt)=true;
            in(cnt)=true;
        else
            % new number
            idx(cnt,1)=i;
            isnum(cnt)=true;
            in(cnt)=true;
        end
    elseif(alf(i))
        % word
        if(in(cnt) && isnum(cnt))
            % terminate number, begin word
            idx(cnt,2)=i-1;
            cnt=cnt+1;
            idx(cnt,1)=i;
            isnum(cnt)=false;
            in(cnt)=true;
        elseif(in(cnt) && ~isnum(cnt))
            % continuing word
            continue;
        else
            % new word
            idx(cnt,1)=i;
            isnum(cnt)=false;
            in(cnt)=true;
        end
    else
        % not either
        if(in(cnt) && isnum(cnt))
            % terminate number
            idx(cnt,2)=i-1;
            cnt=cnt+1;
            in(cnt)=false;
        elseif(in(cnt) && ~isnum(cnt))
            % terminate word
            idx(cnt,2)=i-1;
            cnt=cnt+1;
            in(cnt)=false;
        end
    end
end

% close final word/number
if(in(cnt)); idx(cnt,2)=i; end

% make cellstr
an=cell(1,numel(isnum));
for i=1:numel(isnum)
    if(isnum(i))
        an{i}=str2double(str(idx(i,1):idx(i,2)));
    else
        an{i}=str(idx(i,1):idx(i,2));
    end
end

end

function [words]=getwords(str,delimiter,collapse)
%GETWORDS    Returns a cell array of words from a string
%
%    Usage:    words=getwords('str')
%              words=getwords('str',delimiter)
%              words=getwords('str',delimiter,collapse)
%
%    Description:
%     WORDS=GETWORDS('STR') extracts words in STR and returns them
%     separated into a cellstr array WORDS without any whitespace.
%
%     WORDS=GETWORDS('STR',DELIMITER) separates words in STR using the
%     single character DELIMITER. The default is '' or [] which indicates
%     any whitespace character.
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
%     % Break up a sentence:
%     getwords('This example is pretty dumb!')
%       ans = 
%       'This'    'example'    'is'    'pretty'    'dumb!'
%
%     % Turn off multi-delimiter collapsing to allow handling empty words:
%     getwords('the answer is  !',[],false)
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
%        Jan.  3, 2011 - use isstring
%        Nov.  1, 2011 - doc update
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Nov.  1, 2011 at 13:00 GMT

% todo:

% check nargin
error(nargchk(1,3,nargin));

% check str
if(~isstring(str))
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

function [lgc]=isstring(str)
%ISSTRING    True for a string (row vector) of characters
%
%    Usage:    lgc=isstring(str)
%
%    Description:
%     LGC=ISSTRING(STR) returns TRUE if STR is a string (ie row vector) of
%     characters.  This means ISCHAR(STR) must be TRUE and SIZE(STR,1)==1.
%
%    Notes:
%
%    Examples:
%     % A 2x2 character array will return FALSE:
%     isstring(repmat('a',2,2))
%
%    See also: ISCHAR

%     Version History:
%        Sep. 13, 2010 - initial version (added docs)
%        Oct. 11, 2010 - allow empty string
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Oct. 11, 2010 at 11:00 GMT

lgc=ischar(str) && ndims(str==2) ...
    && (size(str,1)==1 || isequal(size(str),[0 0]));

end
