function [ok]=install_seizmo_mmap(mypath,varargin)
%INSTALL_SEIZMO_MMAP    Check & install M_Map for SEIZMO
%
%    Usage:    ok=install_seizmo_mmap(mypath)
%              ok=install_seizmo_mmap(mypath,type)
%
%    Description:
%     OK=INSTALL_SEIZMO_MMAP(MYPATH) installs the M_Map files
%     associated with the current SEIZMO on the path.  The path is saved to
%     the pathdef.m file that the path was loaded from at startup.  GSHHS
%     coastline/border/river files are needed to complete the installation.
%     See the Examples section below to download these too (WARNING: 106
%     megabyte download).
%
%     OK=INSTALL_SEIZMO_MMAP(MYPATH,TYPE) allows editing the pathdef.m
%     save preference.  See SAVEPATH_SEIZMO for details.  The default is
%     no input to SAVEPATH_SEIZMO.
%
%    Notes:
%     - To complete the installation GMT coastline, rivers, & border files
%       must be installed on the path.  They are available from here:
%        http://www.ngdc.noaa.gov/mgg/shorelines/data/gshhs/oldversions/
%       I would suggest grabbing gshhs_1.10.zip under the 1.10 version for
%       now.  Support for 2+ in M_Map is coming but not currently
%       supported.  Extract the *.b files on to the Matlab path somewhere.
%       If you do not put them on the path, you will need to edit the
%       m_gshhs*.m files!  See the examples section below for help doing
%       this in Matlab/Octave.
%
%    Examples:
%     % Download the GSHHS package & extract the files onto path:
%     % WARNING! This is a 100+ MEGABYTE download...
%     url='http://www.ngdc.noaa.gov/mgg/shorelines/data/gshhs/';
%     url=[url 'oldversions/version1.10/gshhs_1.10.zip'];
%     unzip(url,'m_map'); % extracting into the m_map folder
%
%    See also: INSTALL_SEIZMO, INSTALL_SEIZMO_MATTAUP, INSTALL_SEIZMO_WW3,
%              INSTALL_SEIZMO_CORE, INSTALL_SEIZMO_OPTIONAL

%     Version History:
%        Dec. 30, 2010 - initial version
%        Jan.  1, 2011 - added msg about gshhs files
%        Jan.  4, 2011 - improved docs about GSHHS files
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Jan.  4, 2011 at 15:25 GMT

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
    '& border files must be installed on the path.  They are\n' ...
    'available from here:\n' ...
    ' http://www.ngdc.noaa.gov/mgg/shorelines/data/gshhs/oldversions/\n'...
    'I would suggest grabbing gshhs_1.10.zip under the 1.10 version\n' ...
    'for now.  Support for 2+ in M_Map is coming but not currently\n' ...
    'supported.  Extract the *.b files on to the Matlab path\n' ...
    'anywhere.  If you do not put them on the path, you will need to\n' ...
    'edit M_Map''s m_gshhs*.m files to point to them!\n' ...
    '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n']);

end
