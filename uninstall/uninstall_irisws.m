function [ok]=uninstall_irisws()
%UNINSTALL_IRISWS    Uninstalls the currently installed IRIS web services
%
%    Usage:    ok=uninstall_irisws
%
%    Description:
%     OK=UNINSTALL_IRISWS removes the jar-file associated with IRIS web
%     services from the path and tries to save the edited path.  OK is TRUE
%     if the uninstall succeeded.
%
%    Notes:
%     - Uses the java class TraceData to detect the toolbox.
%
%    Examples:
%     % Reinstall IRIS web services:
%     uninstall_irisws & webinstall_irisws
%
%    See also: WEBINSTALL_IRISWS, UNINSTALL_NJTBX, WEBINSTALL_NJTBX,
%              UNINSTALL_MMAP, WEBINSTALL_MMAP, UNINSTALL_GSHHG,
%              WEBINSTALL_GSHHG, UNINSTALL_EXPORTFIG, WEBINSTALL_EXPORTFIG,
%              UNINSTALL_EXTRAS, WEBINSTALL_EXTRAS, UNINSTALL_TAUP,
%              WEBINSTALL_TAUP, UNINSTALL_SEIZMO, INSTALL_SEIZMO

%     Version History:
%        Feb. 20, 2014 - initial version
%        Mar.  1, 2014 - only remove specific jars, java detection
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2014 at 15:25 GMT

% todo:

% directory separator
fs=filesep;

% does IRISWS-*.jar exist on javaclasspath?
ok=true;
if(exist('edu.iris.dmc.extensions.fetch.TraceData','class'))
    path=fileparts(which('irisFetch')); % root directory
else
    % not found, so toolbox not installed...
    return;
end

% check that java pkg is installed
java_in_octave=true;
if(exist('OCTAVE_VERSION','builtin')==5 && isempty(ver('java')))
    java_in_octave=false;
end

% clear the dynamic java path
jars{1}=[path fs 'IRIS-WS-2.0.6.jar'];
for i=1:numel(jars)
    if(java_in_octave && ismember(jars{i},javaclasspath))
        javarmpath(jars{i});
    end
end

% find classpath.txt
sjcp=which('classpath.txt');
if(isempty(sjcp)); return; end

% read classpath.txt
s2=textread(sjcp,'%s','delimiter','\n','whitespace','');

% detect offending classpath.txt lines
injcp=~cellfun('isempty',strfind(s2,jars{1}));

% only remove if necessary
if(sum(injcp)>0)
    % inform user about which lines are to be removed
    fprintf(' Removing the following lines:\n');
    fprintf('  %s\n',s2{injcp});
    fprintf(' from:\n  %s\n',sjcp);
    
    % the hard part (remove offending lines from classpath.txt)
    fid=fopen(sjcp,'w');
    if(fid<0)
        warning('seizmo:uninstall_irisws:failedToOpen',...
            'Cannot edit classpath.txt!');
        disp('#######################################################');
        disp('You must have Root/Administrator privileges to edit');
        disp('the classpath.txt file.  To fully uninstall IRISWS');
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
