function [events]=read_usgs_fm(file,hdrlines,flag)
%READ_USGS_FM    Reads FM format text from the USGS website
%
%    Usage:    events=read_usgs_fm(file,headerlines)
%              events=read_usgs_fm(string,headerlines,true)
%
%    Description:
%     !!! WARNING !!!
%     UNFORTANATELY THE USGS HAVE ELIMINATED THE SOPAR SITE IN FAVOR OF
%     BURYING SOURCE PARAMETER INFO INTO INDIVIDUAL EVENT PAGES (SO IT IS
%     NO LONGER ACCESSIBLE VIA SEARCH) WHILE EXPANDING THE SEARCH
%     CATEGORIES FOR "COOLER" EARTHQUAKE INFORMATION LIKE SOCIAL MEDIA
%     STATISTICS.  ENJOY YOUR "OMG EARTHQUAKE!" TWITTER WAVES USGS.
%     I WILL KEEP THIS FUNCTION AROUND FOR THOSE THAT HAVE SAVED SOPAR
%     OUTPUT AS A FILE AND IN CASE THE USGS REGROWS A BRAIN.  FOR NOW IF
%     YOU WANT SOURCE PARAMETERS USE THE NEW FUNCTION PARSE_ISC_FM.
%     !!! WARNING !!!
%
%     EVENTS=READ_USGS_FM(FILE,HEADERLINES) reads in an ascii file in USGS
%     FM format containing moment tensor & fault plane solutions.  FILE
%     may be omitted to allow the user to graphically select the file.
%     HEADERLINES allows skipping lines at the start of the file and is 0
%     by default.  Note that only lines that are 146 characters in length
%     are parsed by this routine so you probably don't need to set the
%     HEADERLINES input.  Entries without a depth or scalar moment are
%     omitted.
%
%     EVENTS=READ_USGS_FM(STRING,HEADERLINES,TRUE) allows directly passing
%     the FM formatted text into the routine for parsing.  STRING should be
%     a vector.  This might be useful for reading directly from a website
%     via URLREAD.
%
%    Notes:
%     - Exponent in Nm is converted to dyne-cm (SI to CGS standard) to
%       match GlobalCMT catalog units.
%     - Centroid fields are set to be equivalent to epicenter since there
%       is no centroid determination for USGS moment tensors.
%     - mb, ms are set to 0
%     - Source function half-durations are determined with MO2HD which
%       uses the Harvard empirical scaling.
%
%    Examples:
%     % Read in the entire USGS moment tensor catalog (minus GlobalCMT
%     % entries) via their website (note the string flag must be set!):
%     mts=read_usgs_mt(urlread( ...
%         ['http://neic.usgs.gov/cgi-bin/sopar/sopar.cgi' ...
%          '?GS=1&OTHER=4&FILEFORMAT=1']),[],true);
%
%    See also: READ_USGS_MT, READNDK, FINDCMT, FINDCMTS, PARSE_ISC_FM

%     Version History:
%        June 14, 2011 - initial version
%        Jan. 27, 2014 - abs path exist fix
%        Feb.  9, 2014 - use readtxt, noted uselessness of usgs
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  9, 2014 at 13:50 GMT

% todo:

% check nargin
error(nargchk(0,3,nargin));

% default flag
if(nargin<2 || isempty(hdrlines)); hdrlines=0; end
if(nargin<3 || isempty(flag)); flag=false; end

% directory separator
fs=filesep;

% skip if string input
if(~flag)
    % read in ndk file
    if(nargin<1); file=[]; end
    txt=readtxt(file,{'*.fm;*.FM' 'USGS FM Files (*.fm,*.FM)';
            '*.*' 'All Files (*.*)'});
else
    % just copy file to txt
    if(nargin<1 || isempty(file) || ~isstring(file))
        error('seizmo:read_usgs_fm:emptyStr',...
            'STRING must be non-empty!');
    else
        txt=file;
    end
end

% delete carriage return characters
txt(txt==13)=[];

% split lines
lines=getwords(txt,sprintf('\n'));

% skip header lines
lines=lines(hdrlines+1:end);

% skip all lines that are not 146 characters long
lines=lines(cellfun('prodofsize',lines)==146);
if(isempty(lines))
    error('seizmo:read_usgs_fm:malformedFile',...
        ['USGS FM file: %s\n'...
        'Looks like no lines are 146 characters long!'],file);
end

% now go back to a char array so we can extract certain column sections
txt=char(lines);

% remove bad (using eigval1 field)
good=~isnan(str2double(cellstr(txt(1:end,65:68))));

% actual parsing of text (slooow)
events.year=str2num(txt(good,2:5));
events.month=str2num(txt(good,7:8));
events.day=str2num(txt(good,10:11));
events.hour=str2num(txt(good,13:14));
events.minute=str2num(txt(good,16:17));
events.seconds=str2num(txt(good,19:23));
events.latitude=str2num(txt(good,25:31));
events.longitude=str2num(txt(good,33:40));
events.catalog=cellstr(txt(good,42:44)); % hypocenter catalog
events.depth=str2double(cellstr(txt(good,46:50))); % missing values
events.centroidtime=zeros(size(events.day)); % setting to 0 (no centroid)
events.centroidlat=events.latitude;  % no centroid determination
events.centroidlon=events.longitude; % so we set this to the same
events.centroiddep=events.depth;     % as the pde locations
events.mw=str2double(cellstr(txt(good,52:54))); % missing values
events.mb=zeros(size(events.mw)); % not
events.ms=zeros(size(events.mw)); % given
events.scalarmoment=str2double(cellstr(txt(good,56:58))); % missing values
events.exponent=str2double(cellstr(txt(good,61:62)))+7; % missing values
% using harvard empirical relationship to get halfdur
events.srcfuncdur=mo2hd(events.scalarmoment.*10.^events.exponent);
events.eigval1=str2double(cellstr(txt(good,65:68))); % missing values
events.plunge1=str2num(txt(good,70:71));
events.azimuth1=str2num(txt(good,73:75));
events.eigval2=str2double(cellstr(txt(good,77:81))); % missing values
events.plunge2=str2double(cellstr(txt(good,83:84))); % missing values
events.azimuth2=str2double(cellstr(txt(good,86:88))); % missing values
events.eigval3=str2double(cellstr(txt(good,90:94))); % missing values
events.plunge3=str2num(txt(good,96:97));
events.azimuth3=str2num(txt(good,99:101));
events.strike1=str2num(txt(good,103:105));
events.dip1=str2num(txt(good,107:108));
events.rake1=str2num(txt(good,110:113));
events.strike2=str2num(txt(good,115:117));
events.dip2=str2num(txt(good,119:120));
events.rake2=str2num(txt(good,122:125));
events.srcinvtype=cellstr(txt(good,131:134)); % tensor source catalog

% get harvard moment tensor components
[events.mrr,events.mtt,events.mpp,events.mrt,events.mrp,events.mtp]=...
    mt_g2c(tpb2mt([events.eigval1 events.plunge1 events.azimuth1],...
    [events.eigval3 events.plunge3 events.azimuth3],...
    [events.eigval2 events.plunge2 events.azimuth2]));

end
