function [ok]=uninstall_njtbx()
%UNINSTALL_NJTBX    Uninstalls the currently installed njtbx
%
%    Usage:    ok=uninstall_njtbx
%
%    Description:
%     OK=UNINSTALL_NJTBX removes the directories and jar-files associated
%     with the njtbx toolbox from the path and tries to save the edited
%     path.  OK is the savepath exit status.
%
%    Notes:
%     - Uses the function 'nj_time' to detect the toolbox.
%
%    Examples:
%     % Update njtbx:
%     uninstall_njtbx & webinstall_njtbx
%
%    See also: WEBINSTALL_NJTBX, UNINSTALL_MMAP, WEBINSTALL_MMAP,
%              UNINSTALL_GSHHS, WEBINSTALL_GSHHS

%     Version History:
%        Feb. 14, 2012 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 14, 2012 at 15:25 GMT

% todo:

% does nj_time exist?
if(exist('nj_time','file'))
    path=fileparts(fileparts(which('nj_time'))); % root directory
    rmpath(fullfile(path,'njTBX-2.0','Utilities'));
    rmpath(fullfile(path,'njTBX-2.0'));
    rmpath(fullfile(path,'njFunc'));
    rmpath(fullfile(path,'examples'));
    rmpath(path);
    ok=savepath;
end

% clear the dynamic java path
jars=dir([path '*.jar']);
javarmpath(fullfile(path,jars(1).name));
javarmpath(fullfile(path,jars(2).name));

% find classpath.txt
sjcp=which('classpath.txt');
if(isempty(sjcp))
    warning('seizmo:uninstall_njtbx:noJavaClassPath',...
        'Octave has no classpath.txt to remove .jar files from!');
    disp('Skipping remainder of njtbx uninstall!');
    ok=false;
    return;
end

% read classpath.txt
s2=textread(sjcp,'%s','delimiter','\n','whitespace','');

% detect offending classpath.txt lines
yn=~cellfun('isempty',strfind(s2,path));

% inform user about which lines are to be removed
disp('Removing the following lines:');
fprintf('%s\n',s2{yn});
fprintf('\nfrom:\n%s\n\n',sjcp);

% only remove if necessary
if(sum(yn))
    % the hard part (remove offending lines from classpath.txt)
    fid=fopen(sjcp,'w');
    if(fid<0)
        warning('seizmo:uninstall_njtbx:failedToOpen',...
            'Cannot edit classpath.txt!');
        disp('You must have Root/Administrator privileges to edit');
        disp('the classpath.txt file.  To fully uninstall njtbx');
        disp('you need to remove the following line(s):');
        disp(strrep(sprintf('%s\n',s2{yn}),'\','\\'));
        disp(' ');
        disp('from your Matlab''s classpath.txt located here:');
        disp(strrep(sjcp,'\','\\'));
        disp(' ');
        disp(' ');
        disp('This may be done by contacting your System Admin if you');
        disp('do not have Root/Administrator privileges.  Afterwards,');
        disp('please restart Matlab to complete the uninstallation!');
        ok=false;
    else
        fseek(fid,0,'bof');
        for i=find(~yn)'
            fprintf(fid,'%s\n',s2{i});
        end
        fclose(fid);
        ok=ok & true;
    end
else
    ok=ok & true;
end

end
