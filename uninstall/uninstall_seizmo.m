function [ok]=uninstall_seizmo()
%UNINSTALL_SEIZMO    Uninstalls the currently installed SEIZMO
%
%    Usage:    ok=uninstall_seizmo
%
%    Description:
%     OK=UNINSTALL_SEIZMO removes the directories and jar-files associated
%     with the SEIZMO toolbox from the path and tries to save the edited
%     path.  OK is TRUE if the uninstall succeeded.
%
%    Notes:
%     - Uses the location of the 'seizmodef' function to find the toolbox.
%     - This function is not in the top level directory as that would
%       lead to running 'uninstall_seizmo' from a not-installed SEIZMO
%       directory when trying to install that SEIZMO directory.
%
%    Examples:
%     % Update SEIZMO
%     install_seizmo  % calls uninstall_seizmo
%
%    See also: ABOUT_SEIZMO, SEIZMO, INSTALL_SEIZMO, UNINSTALL_NJTBX,
%              UNINSTALL_MMAP, UNINSTALL_GSHHS, UNINSTALL_EXPORTFIG

%     Version History:
%        Jan.  1, 2011 - initial version
%        June 24, 2011 - octave bugfix: uninstall seizmo from path before
%                        exiting when no classpath file found
%        Feb. 15, 2012 - doc update, cleaner code, call external component
%                        uninstallers, only use javarmpath when needed,
%                        don't force failure for octave
%        Feb. 16, 2012 - export_fig is externally managed
%        Feb. 27, 2012 - multi-jar mattaup update
%        Mar.  8, 2012 - fix mattaup multi-jar breakage
%        Mar. 28, 2013 - add ocean/xc directories
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 28, 2013 at 15:25 GMT

% todo:

% ask to install external components
% - has to be done before removing 'seizmo/uninstall' directory
ok=true;
reply=input('Uninstall njTBX? Y/N [Y]: ','s');
if(isempty(reply) || strncmpi(reply,'y',1))
    ok=ok & uninstall_njtbx;
end
reply=input('Uninstall M_Map? Y/N [Y]: ','s');
if(isempty(reply) || strncmpi(reply,'y',1))
    ok=ok & uninstall_mmap;
    reply=input('Uninstall GSHHS? Y/N [Y]: ','s');
    if(isempty(reply) || strncmpi(reply,'y',1))
        ok=ok & uninstall_gshhs;
    end
end
reply=input('Uninstall export_fig? Y/N [Y]: ','s');
if(isempty(reply) || strncmpi(reply,'y',1))
    ok=ok & uninstall_exportfig;
end

% does seizmodef exist?
% - this is the "kernel" of seizmo
fs=filesep;
if(exist('seizmodef','file'))
    path=fileparts(fileparts(which('seizmodef'))); % root directory
    rmpath(...
        path,...
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
        [path fs 'ocean'],...
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
        [path fs 'xc'],...
        [path fs 'xcalign'],...
        [path fs 'mattaup']);
    ok=ok & ~savepath;
    if(~ok)
        warning('seizmo:uninstall_seizmo:noPermission',...
            'Could not edit path!');
    end
else
    % not found, so toolbox not installed...
    return;
end

% clean out mattaup jars from dynamic java path
jar=dir(fullfile(path,'mattaup','lib','*.jar'));
for i=1:numel(jar)
    if(ismember(fullfile(path,'mattaup','lib',...
            jar(i).name),javaclasspath))
        javarmpath(fullfile(path,'mattaup','lib',jar(i).name));
    end
end

% find classpath.txt
sjcp=which('classpath.txt');
if(isempty(sjcp)); return; end

% read classpath.txt
s2=textread(sjcp,'%s','delimiter','\n','whitespace','');

% detect offending classpath.txt lines
injcp=~cellfun('isempty',strfind(s2,fullfile(path,'mattaup','lib','')));

% only remove if necessary
if(sum(injcp)>0)
    % inform user about which lines are to be removed
    fprintf(' Removing the following lines:\n');
    fprintf('  %s\n',s2{injcp});
    fprintf(' from:\n  %s\n',sjcp);
    
    % the hard part (remove offending lines from classpath.txt)
    fid=fopen(sjcp,'w');
    if(fid<0)
        warning('seizmo:uninstall_seizmo:failedToOpen',...
            'Cannot edit classpath.txt!');
        disp('#######################################################');
        disp('You must have Root/Administrator privileges to edit');
        disp('the classpath.txt file.  To fully uninstall SEIZMO');
        disp('you need to remove the following line(s):');
        disp(strrep(sprintf('%s\n',s2{injcp}),'\','\\'));
        disp(' ');
        disp('from your Matlab''s classpath.txt located here:');
        disp(strrep(sjcp,'\','\\'));
        disp(' ');
        disp('This may be done by contacting your System Admin if you');
        disp('do not have Root/Administrator privileges.  Afterwards,');
        disp('please restart Matlab to complete the uninstallation!');
        disp('#######################################################');
        ok=false;
    else
        fseek(fid,0,'bof');
        fprintf(fid,'%s\n',s2{~injcp});
        fclose(fid);
    end
end

end
