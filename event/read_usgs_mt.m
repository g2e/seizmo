function [events]=read_usgs_mt(file,hdrlines,flag)
%READ_USGS_MT    Reads MT format text from the USGS website
%
%    Usage:    events=read_usgs_mt(file,headerlines)
%              events=read_usgs_mt(string,headerlines,true)
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
%     EVENTS=READ_USGS_MT(FILE,HEADERLINES) reads in an ascii file in USGS
%     MT formatting containing moment tensor solutions.  FILE may be
%     omitted to allow the user to graphically select the file.
%     HEADERLINES allows skipping lines at the start of the file and is 0
%     by default.  Note that only lines that are 146 characters in length
%     are parsed by this routine so you probably don't need to set the
%     HEADERLINES input.
%
%     EVENTS=READ_USGS_MT(STRING,HEADERLINES,TRUE) allows directly passing
%     the ascii text into the routine for parsing.  STRING should be a
%     vector.  This might be useful for reading directly from a website.
%
%    Notes:
%     - Exponent in Nm is converted to dyne-cm (SI to CGS standard) to
%       match GlobalCMT catalog units.
%     - Centroid fields are set to be equivalent to hypocenter since there
%       is no centroid determination for USGS moment tensors
%     - Magnitudes not given are set to 0
%     - Source function half-durations are determined with MO2HD which
%       uses the Harvard empirical scaling
%     - Principal axes are determined using moment tensor components
%
%    Examples:
%     % Read in usgs.momten file and find moment tensors from 1990:
%     mts=findcmts('st',[1990 1],'et',[1990 365],...
%         'catalog',read_usgs_mt(which('usgs.mt')));
%
%     % Read in the entire USGS moment tensor catalog (minus GlobalCMT
%     % entries) via their website (note the string flag must be set!):
%     mts=read_usgs_mt(urlread( ...
%         ['http://neic.usgs.gov/cgi-bin/sopar/sopar.cgi' ...
%          '?GS=1&OTHER=4&FILEFORMAT=3']),[],true);
%
%    See also: READ_USGS_FM, READNDK, FINDCMT, FINDCMTS, PARSE_ISC_FM

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

% skip if string input
if(~flag)
    % read in ndk file
    if(nargin<1); file=[]; end
    txt=readtxt(file,{'*.mt;*.MT' 'USGS MT Files (*.mt,*.MT)';
        '*.*' 'All Files (*.*)'});
else
    % just copy file to txt
    if(nargin<1 || isempty(file) || ~isstring(file))
        error('seizmo:read_usgs_mt:emptyStr',...
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
    error('seizmo:read_usgs_mt:malformedFile',...
        ['USGS MT file: %s\n'...
        'Looks like no lines are 146 characters long!'],file);
end

% now go back to a char array so we can extract certain column sections
txt=char(lines);

% actual parsing of text (slooow)
events.year=str2num(txt(1:end,2:5));
events.month=str2num(txt(1:end,7:8));
events.day=str2num(txt(1:end,10:11));
events.hour=str2num(txt(1:end,13:14));
events.minute=str2num(txt(1:end,16:17));
events.seconds=str2num(txt(1:end,19:23));
events.centroidtime=zeros(size(events.day)); % setting to 0 (no centroid)
events.latitude=str2num(txt(1:end,25:31));
events.longitude=str2num(txt(1:end,33:40));
events.catalog=cellstr(txt(1:end,42:44)); % hypocenter catalog
events.depth=str2num(txt(1:end,46:50));
events.centroidlat=events.latitude;  % no centroid determination
events.centroidlon=events.longitude; % so we set this to the same
events.centroiddep=events.depth;     % as the pde locations
events.mw=str2num(txt(1:end,52:54));
events.mb=zeros(size(events.mw)); % not
events.ms=zeros(size(events.mw)); % given
events.scalarmoment=str2num(txt(1:end,56:58));
events.exponent=str2num(txt(1:end,61:62))+7; % 1Nm => 1e7 dyne-cm
events.srcfuncdur=str2double(cellstr(txt(1:end,64:67))); % missing values
% using harvard empirical relationship to fill-in usgs mt w/o halfdur
bad=isnan(events.srcfuncdur);
events.srcfuncdur(bad)=mo2hd(...
    events.scalarmoment(bad).*10.^events.exponent(bad));
events.srcinvtype=cellstr(txt(1:end,69:72)); % tensor source catalog
%events.exponent=str2num(txt(1:end,74:75)); % repeated exponent
events.mrr=str2num(txt(1:end,77:81));
events.mtt=str2num(txt(1:end,83:87));
events.mpp=str2num(txt(1:end,89:93));
events.mrt=str2num(txt(1:end,95:99));
events.mrp=str2num(txt(1:end,101:105));
events.mtp=str2num(txt(1:end,107:111));

% calculate principal axes
[t,p,b]=mt2tpb(mt_c2v(...
    [events.mrr events.mtt events.mpp events.mrt events.mrp events.mtp]));
events.eigval1=t(:,1);
events.eigval2=b(:,1);
events.eigval3=p(:,1);
events.plunge1=t(:,2);
events.plunge2=b(:,2);
events.plunge3=p(:,2);
events.azimuth1=t(:,3);
events.azimuth2=b(:,3);
events.azimuth3=p(:,3);

end
