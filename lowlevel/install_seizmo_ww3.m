function [ok]=install_seizmo_ww3()
%INSTALL_SEIZMO_WW3    Checks that njtbx is installed
%
%    Usage:    ok=install_seizmo_ww3()
%
%    Description:
%     OK=INSTALL_SEIZMO_WW3() checks that njtbx is installed.  If not
%     a warning is given indicating where to go to get it.  OK is TRUE if
%     njtbx is found and FALSE otherwise.
%
%    Notes:
%
%    Examples:
%
%    See also: INSTALL_SEIZMO, INSTALL_SEIZMO_MMAP, INSTALL_SEIZMO_MATTAUP,
%              INSTALL_SEIZMO_CORE, INSTALL_SEIZMO_OPTIONAL

%     Version History:
%        Dec. 30, 2010 - initial version
%        Feb. 14, 2012 - switch from read_grib to njtbx
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb. 14, 2012 at 15:25 GMT

% todo:

% check that necessary njtbx functions are available
if(~exist('mDataset','file') || ~exist('getVars','file') ...
        || ~exist('nj_time','file'))
    warning('seizmo:install_seizmo_ww3:noNJTBX',...
        ['NJTBX is not installed!  You can get it here:\n' ...
        'http://sourceforge.net/apps/trac/njtbx/' ...
        'wiki/DownloadNjtbx-current']);
    ok=false;
else
    ok=true;
end

end
