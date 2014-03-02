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
%     - Uses the 'validseizmo' function location to find the toolbox.
%     - This function is not in the top level directory as that would
%       lead to running 'uninstall_seizmo' from a not-installed SEIZMO
%       directory when trying to install that SEIZMO directory.
%
%    Examples:
%     % Update SEIZMO & extras (must be in the main seizmo directory):
%     delete('*.zip');  % forces redownload of extra packages
%     install_seizmo    % calls uninstall_seizmo
%
%    See also: ABOUT_SEIZMO, SEIZMO, INSTALL_SEIZMO, UNINSTALL_NJTBX,
%              UNINSTALL_MMAP, UNINSTALL_GSHHG, UNINSTALL_EXPORTFIG,
%              UNINSTALL_TAUP, UNINSTALL_IRISWS, UNINSTALL_EXTRAS

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
%        Jan. 15, 2014 - update for gshhs to gshhg rename, use validseizmo
%                        to locate the toolbox (octave workaround), update
%                        example for updating extra stuff too
%        Feb. 25, 2014 - uninstall_taup support, uninstall_irisws support
%        Feb. 27, 2014 - uninstall_extras support
%        Mar.  1, 2014 - updated See also list, savepath only called if
%                        needed
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar.  1, 2014 at 15:25 GMT

% todo:

% ask to install external components
% - has to be done before removing 'seizmo/uninstall' directory
ok=true;
reply=input('Uninstall TauP? Y/N [Y]: ','s');
if(isempty(reply) || strncmpi(reply,'y',1))
    ok=ok & uninstall_taup;
end
reply=input('Uninstall IRISWS? Y/N [Y]: ','s');
if(isempty(reply) || strncmpi(reply,'y',1))
    ok=ok & uninstall_irisws;
end
reply=input('Uninstall njTBX? Y/N [Y]: ','s');
if(isempty(reply) || strncmpi(reply,'y',1))
    ok=ok & uninstall_njtbx;
end
reply=input('Uninstall M_Map? Y/N [Y]: ','s');
if(isempty(reply) || strncmpi(reply,'y',1))
    ok=ok & uninstall_mmap;
    reply=input('Uninstall GSHHG? Y/N [Y]: ','s');
    if(isempty(reply) || strncmpi(reply,'y',1))
        ok=ok & uninstall_gshhg;
    end
end
reply=input('Uninstall export_fig? Y/N [Y]: ','s');
if(isempty(reply) || strncmpi(reply,'y',1))
    ok=ok & uninstall_exportfig;
end
reply=input('Uninstall Extras? Y/N [Y]: ','s');
if(isempty(reply) || strncmpi(reply,'y',1))
    ok=ok & uninstall_extras;
end

% does validseizmo exist?
fs=filesep;
if(exist('validseizmo','file'))
    path=fileparts(fileparts(which('validseizmo'))); % root directory
    seizmo_path={path ...
        [path fs 'lowlevel'] ...
        [path fs 'uninstall'] ...
        [path fs 'behavior'] ...
        [path fs 'toc'] ...
        [path fs 'rw'] ...
        [path fs 'hdr'] ...
        [path fs 'sz'] ...
        [path fs 'misc'] ...
        [path fs 'time'] ...
        [path fs 'position'] ...
        [path fs 'audio'] ...
        [path fs 'cmap'] ...
        [path fs 'cmb'] ...
        [path fs 'cmt'] ...
        [path fs 'decon'] ...
        [path fs 'event'] ...
        [path fs 'filtering'] ...
        [path fs 'fixes'] ...
        [path fs 'fk'] ...
        [path fs 'ftran'] ...
        [path fs 'gui'] ...
        [path fs 'invert'] ...
        [path fs 'mapping'] ...
        [path fs 'models'] ...
        [path fs 'multi'] ...
        [path fs 'noise'] ...
        [path fs 'ocean'] ...
        [path fs 'pick'] ...
        [path fs 'plotting'] ...
        [path fs 'resampling'] ...
        [path fs 'response'] ...
        [path fs 'shortnames'] ...
        [path fs 'solo'] ...
        [path fs 'sphpoly'] ...
        [path fs 'synth'] ...
        [path fs 'tomo'] ...
        [path fs 'topo'] ...
        [path fs 'tpw'] ...
        [path fs 'ttcorrect'] ...
        [path fs 'win'] ...
        [path fs 'ww3'] ...
        [path fs 'ws'] ...
        [path fs 'xc'] ...
        [path fs 'xcalign'] ...
        [path fs 'mattaup']};
    rmpath(seizmo_path{:});
    if(is_on_static_path(seizmo_path{:}))
        ok=ok & ~savepath;
    end
else
    % not found, so toolbox not installed...
    return;
end

end


function [lgc]=is_on_static_path(varargin)
% find pathdef.m
spd=which('pathdef.m');

% read pathdef.m
s=textread(spd,'%s','delimiter','\n','whitespace','');

% detect offending pathdef.m lines
for i=1:nargin
    lgc=any(~cellfun('isempty',strfind(s,varargin{i})));
    if(lgc); return; end
end
end

