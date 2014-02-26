function [ok]=uninstall_njtbx()
%UNINSTALL_NJTBX    Uninstalls the currently installed njTBX
%
%    Usage:    ok=uninstall_njtbx
%
%    Description:
%     OK=UNINSTALL_NJTBX removes the directories and jar-files associated
%     with the njTBX toolbox from the path and tries to save the edited
%     path.  OK is TRUE if the uninstall succeeded.
%
%    Notes:
%     - Uses the location of function 'nj_time' to detect the toolbox.
%
%    Examples:
%     % Update njtbx:
%     uninstall_njtbx & webinstall_njtbx
%
%    See also: WEBINSTALL_NJTBX, UNINSTALL_MMAP, WEBINSTALL_MMAP,
%              UNINSTALL_GSHHG, WEBINSTALL_GSHHG, UNINSTALL_EXPORTFIG,
%              WEBINSTALL_EXPORTFIG, UNINSTALL_EXTRAS, WEBINSTALL_EXTRAS,
%              UNINSTALL_IRISWS, WEBINSTALL_IRISWS, UNINSTALL_TAUP,
%              WEBINSTALL_TAUP, UNINSTALL_SEIZMO, INSTALL_SEIZMO

%     Version History:
%        Feb. 14, 2012 - initial version
%        Feb. 15, 2012 - handle not installed, flip logic from savepath,
%                        doc update, only use javarmpath when needed,
%                        don't force failure for octave
%        Mar.  8, 2012 - make code changes for clarity
%        Apr. 25, 2012 - fix classpath.txt jar removal
%        Jan. 15, 2014 - updated See also list
%        Feb. 20, 2014 - fixed dynamic java path uninstall, no fullfile
%                        use, updated see also list
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 20, 2014 at 15:25 GMT

% todo:

% directory separator
fs=filesep;

% does nj_time exist?
if(exist('nj_time','file'))
    path=fileparts(fileparts(which('nj_time'))); % root directory
    rmpath([path fs 'njTBX-2.0' fs 'Utilities']);
    rmpath([path fs 'njTBX-2.0']);
    rmpath([path fs 'njFunc']);
    rmpath([path fs 'examples']);
    rmpath(path);
    ok=~savepath;
else
    % not found, so toolbox not installed...
    ok=true;
    return;
end

% clear the dynamic java path
jarpath=fileparts(path);
jars=dir([jarpath fs '*.jar']);
for i=1:numel(jars)
    if(ismember([jarpath fs jars(i).name],javaclasspath))
        javarmpath([jarpath fs jars(i).name]);
    end
end

% find classpath.txt
sjcp=which('classpath.txt');
if(isempty(sjcp)); return; end

% read classpath.txt
s2=textread(sjcp,'%s','delimiter','\n','whitespace','');

% detect offending classpath.txt lines
injcp=~cellfun('isempty',strfind(s2,jarpath));

% only remove if necessary
if(sum(injcp)>0)
    % inform user about which lines are to be removed
    fprintf(' Removing the following lines:\n');
    fprintf('  %s\n',s2{injcp});
    fprintf(' from:\n  %s\n',sjcp);
    
    % the hard part (remove offending lines from classpath.txt)
    fid=fopen(sjcp,'w');
    if(fid<0)
        warning('seizmo:uninstall_njtbx:failedToOpen',...
            'Cannot edit classpath.txt!');
        disp('#######################################################');
        disp('You must have Root/Administrator privileges to edit');
        disp('the classpath.txt file.  To fully uninstall njTBX');
        disp('you need to remove the following line(s):');
        disp(strrep(sprintf('%s\n',s2{injcp}),'\','\\'));
        disp(' ');
        disp('from your Matlab''s classpath.txt located here:');
        disp(strrep(sjcp,'\','\\'));
        disp(' ');
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
