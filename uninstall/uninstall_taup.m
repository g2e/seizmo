function [ok]=uninstall_taup()
%UNINSTALL_TAUP    Uninstalls the currently installed IRIS web services
%
%    Usage:    ok=uninstall_taup
%
%    Description:
%     OK=UNINSTALL_TAUP removes the jar-file associated with MatTauP from
%     the path and tries to save the edited path.  OK is TRUE if the
%     uninstall succeeded.
%
%    Notes:
%     - Uses the java class TauP_Time to detect the toolbox.
%
%    Examples:
%     % Reinstall IRIS web services:
%     uninstall_taup & webinstall_irisws
%
%    See also: WEBINSTALL_TAUP, UNINSTALL_NJTBX, WEBINSTALL_NJTBX,
%              UNINSTALL_MMAP, WEBINSTALL_MMAP, UNINSTALL_GSHHG,
%              WEBINSTALL_GSHHG, UNINSTALL_EXPORTFIG, WEBINSTALL_EXPORTFIG,
%              UNINSTALL_EXTRAS, WEBINSTALL_EXTRAS, UNINSTALL_IRISWS,
%              WEBINSTALL_IRISWS, UNINSTALL_SEIZMO, INSTALL_SEIZMO

%     Version History:
%        Feb. 20, 2014 - initial version
%        Mar.  1, 2014 - only remove specific jars, java detection
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2014 at 15:25 GMT

% todo:

% directory separator
fs=filesep;

% does TauP-*.jar exist on javaclasspath?
ok=true;
if(exist('edu.sc.seis.TauP.TauP_Time','class'))
    path=[fileparts(which('tauptime')) fs 'lib']; % jar directory
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
jars{1}=[path fs 'TauP-2.1.1.jar'];
jars{2}=[path fs 'seisFile-1.5.1.jar'];
jars{3}=[path fs 'MatTauP-2.1.1.jar'];
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
injcp=~cellfun('isempty',strfind(s2,jars{1})) ...
    | ~cellfun('isempty',strfind(s2,jars{2})) ...
    | ~cellfun('isempty',strfind(s2,jars{3}));

% only remove if necessary
if(sum(injcp)>0)
    % inform user about which lines are to be removed
    fprintf(' Removing the following lines:\n');
    fprintf('  %s\n',s2{injcp});
    fprintf(' from:\n  %s\n',sjcp);
    
    % the hard part (remove offending lines from classpath.txt)
    fid=fopen(sjcp,'w');
    if(fid<0)
        warning('seizmo:uninstall_taup:failedToOpen',...
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
