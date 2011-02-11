function [ok]=install_seizmo_mmap(mypath,varargin)
%INSTALL_SEIZMO_MMAP    Check & install M_Map for SEIZMO
%
%    Usage:    ok=install_seizmo_mmap(mypath)
%              ok=install_seizmo_mmap(mypath,type)
%
%    Description:
%     OK=INSTALL_SEIZMO_MMAP(MYPATH) installs the M_Map files associated
%     with the current SEIZMO on the path.  The path is saved to the
%     pathdef.m file that the path was loaded from at startup.  GSHHS
%     coastline/border/river files are needed to complete the installation.
%     See the SEIZMO_GSHHS_WEBINSTALL for details on retrieving them.
%
%     OK=INSTALL_SEIZMO_MMAP(MYPATH,TYPE) allows editing the pathdef.m
%     save preference.  See SAVEPATH_SEIZMO for details.  The default is
%     no input to SAVEPATH_SEIZMO.
%
%    Notes:
%
%    Examples:
%
%    See also: INSTALL_SEIZMO, INSTALL_SEIZMO_MATTAUP, INSTALL_SEIZMO_WW3,
%              INSTALL_SEIZMO_CORE, INSTALL_SEIZMO_OPTIONAL,
%              SEIZMO_GSHHS_WEBINSTALL

%     Version History:
%        Dec. 30, 2010 - initial version
%        Jan.  1, 2011 - added msg about gshhs files
%        Jan.  4, 2011 - improved docs about GSHHS files
%        Feb.  5, 2011 - updates related to seizmo_gshhs_webinstall
%        Feb. 10, 2011 - add seizmo_gshhs_webinstall to See also section,
%                        remove notes & example
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 10, 2011 at 15:25 GMT

% todo:

% check nargin
error(nargchk(0,2,nargin));

% default path input
if(nargin<1 || isempty(mypath)); mypath='.'; end

% check path
fs=filesep;
if(~exist(mypath,'dir'))
    error('seizmo:install_seizmo_mmap:badPath',...
        ['SEIZMO directory (' mypath ') does not exist!']);
end

% check for m_map
ok=false;
if(~isempty(ver('m_map')))
    % complain about non-SEIZMO M_Map
    mmappath=fileparts(which('m_coast'));
    warning('seizmo:install_seizmo_mmap:previousInstall',...
        ['M_Map has been previously installed outside of SEIZMO!\n' ...
        'Offending M_Map Path:\n' ...
        '    ' mmappath '\n\n' ...
        'SEIZMO''s M_Map will still be installed and will take\n' ...
        'precedence over the previous installation of M_Map!']);
end

% add seizmo m_map mfiles to path
addpath([mypath fs 'm_map']);
bad=savepath_seizmo(varargin{:});
if(~bad); ok=true; end

% final message about coastlines
fprintf(...
    ['\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n' ...
    'To complete the M_Map installation, GMT coastline, rivers,\n' ...
    '& border files must be installed on the path. See \n' ...
    '<a href="matlab:help seizmo_gshhs_webinstall">' ...
    'seizmo_gshhs_webinstall</a> for more info\n' ...
    '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n']);

end
