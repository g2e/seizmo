function [isc]=parse_isc_origin(file,flag)
%PARSE_ISC_ORIGIN    Parses hypocenter output from the ISC bulletin
%
%    Usage:    isc=parse_isc_origin(file)
%              isc=parse_isc_origin(string,true)
%
%    Description:
%     ISC=PARSE_ISC_ORIGIN(FILE) reads an ascii file containing origin
%     output from the ISC bulletin.  To create your own file go here:
%       http://www.isc.ac.uk/iscbulletin/search/bulletin/
%     and make sure the output format is set to ISF Bulletin and to select
%     Prime Hypocenter output and unselect the output for Headers,
%     Comments, Links, Magnitudes and Phases.  Save the page as an .html
%     file.  FILE should be the filename (and path if necessary) of the
%     saved .html file.  The output (ISC) is a scalar struct with fields
%     for all of the ISC origin datums.  Each field will be NORIGINSx1 and
%     will either be a double array or a cellstr array and the function
%     SSIDX is provided to pull info from scalar structs.
%
%     ISC=PARSE_ISC_ORIGIN(STRING,TRUE) will take the first input to be a
%     character string of the contents of an ISC hypocenter .html file
%     (rather than the filename) if the second input is set to the logical
%     TRUE.  The string should be the same as if the file was read with
%     READTXT (a single row char vector with linefeeds included).  This is
%     also useful when combined with URLREAD to directly grab info from the
%     web rather than having to save it to your filesystem first.
%
%    Notes:
%     - Empty numeric fields have a NaN placeholder and empty string fields
%       are left as empty strings.
%
%    Examples:
%     % Passing no arguments brings up a graphical menu to select a file:
%     isc=parse_isc_origin();
%
%     % Do a test run:
%     link=['http://isc-mirror.iris.washington.edu/cgi-bin/web-db-v4' ...
%           '?request=COMPREHENSIVE&out_format=ISF&bot_lat=&top_lat=' ...
%           '&left_lon=&right_lon=&ctr_lat=&ctr_lon=&radius=' ...
%           '&max_dist_units=deg&searchshape=GLOBAL&srn=&grn=' ...
%           '&start_year=2011&start_month=7&start_day=31' ...
%           '&start_time=00%3A00%3A00&end_year=2011&end_month=8' ...
%           '&end_day=01&end_time=00%3A00%3A00&min_dep=&max_dep=' ...
%           '&min_mag=&max_mag=&req_mag_type=Any&req_mag_agcy=Any' ...
%           '&min_def=&max_def=&prime_only=on'];
%     isc=parse_isc_origin(urlread(link),true)
%
%     % Extract the 3rd origin from the returned isc origin struct:
%     origin=ssidx(isc,3)
%
%    See also: PARSE_ISC_FM, READNDK, READSODEVENTCSV, READ_USGS_MT,
%              READ_USGS_FM, SSIDX

%     Version History:
%        July 26, 2010 - initial version by Erica, added code to allow for
%                        graphical file selection and handling blank fields
%        July 28, 2010 - added documentation
%        Mar.  7, 2011 - mention ssidx
%        Jan. 27, 2014 - abs path exist fix
%        Feb.  9, 2014 - no more header line option, direct string (html)
%                        input allowed with flagged second option, parses
%                        html page directly, updated docs to reflect ISC
%                        webpage changes, altered examples
%
%     Written by Erica Emry (ericae at wustl dot edu)
%                Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  9, 2014 at 17:25 GMT

% todo:

% check nargin
error(nargchk(0,2,nargin));

% default flag
if(nargin<2 || isempty(flag)); flag=false; end

% skip if string input
if(~flag)
    % read in isc origin file
    if(nargin<1); file=[]; end
    txt=readtxt(file,{'*.html;*.HTML' 'HTML Files (*.html,*.HTML)';
        '*.txt;*.TXT' 'TXT Files (*.txt,*.TXT)';
        '*.*' 'All Files (*.*)'});
else
    % just copy file to txt
    if(nargin<1 || isempty(file) || ~isstring(file))
        error('seizmo:parse_isc_origin:emptyStr',...
            'STRING must be non-empty!');
    else
        txt=file;
    end
end

% delete carriage return characters
txt(txt==13)=[];

% separate lines
lines=getwords(txt,sprintf('\n'));

% find data lines
b=find(strcmp('ISC Bulletin',lines));
e=find(strcmp('STOP',lines));

% convert to character array (to extract fixed column sections)
txt0=char(lines(b+1:e-1));

% Parse ISC origin output -
%   Input file should be: 
%       Bulletin Type - Comprehensive
%       Format - IMS1.0
%       No headers, comments, links, secondaries, magnitudes, or phases
%           selected
%
% Record 	Position 	Format 	Description
% Origin Sub-block 		
% 1 	4-7 	a4 	Date
% (header) 	15-18 	a4 	Time
% 	27-29 	a3 	Err
% 	33-35 	a3 	RMS
% 	37-44 	a8 	Latitude
% 	46-54 	a9 	Longitude
% 	57-60 	a4 	Smaj
% 	63-66 	a4 	Smin
% 	69-70 	a2 	Az
% 	72-76 	a5 	Depth
% 	80-82 	a3 	Err
% 	84-87 	a4 	Ndef
% 	89-92 	a4 	Nst
% 	94-96 	a3 	Gap
% 	99-103 	a5 	mdist
% 	106-110 	a5 	Mdist
% 	112-115 	a4 	Qual
% 	19-124 	a6 	Author
% 	131-136 	a6 	OrigID
% Record 	Position 	Format 	Description
% 	1-10 	i4,a1,i2,a1,i2 	epicenter date (yyyy/mm/dd)
% 	12-22 	i2,a1,i2,a1,f5.2 	epicenter time (hh:mm:ss.ss)
% 	23 	a1 	fixed flag (f = fixed origin time solution, blank if not a fixed origin time)
% 	25-29 	f5.2 	origin time error (seconds; blank if fixed origin time)
% 	31-35 	f5.2 	root mean square of time residuals (seconds)
% 	37-44 	f8.4 	latitude (negative for South)
% 	46-54 	f9.4 	longitude (negative for West)
% 	55 	a1 	fixed flag (f = fixed epicenter solution, blank if not a fixed epicenter solution)
% 	56-60 	f5.1 	semi-major axis of 90% ellipse or its estimate (km, blank if fixed epicenter)
% 	62-66 	f5.1 	semi-minor axis of 90% ellipse or its estimate (km, blank if fixed epicenter)
% 	68-70 	i3 	strike (0 <= x <= 360) of error ellipse clock-wise from North (degrees)
% 	72-76 	f5.1 	depth (km)
% 	77 	a1 	fixed flag (f = fixed depth station, d = depth phases, blank if not a fixed depth)
% 	79-82 	f4.1 	depth error 90% (km; blank if fixed depth)
% 	84-87 	i4 	number of defining phases
% 	89-92 	i4 	number of defining stations
% 	94-96 	i3 	gap in azimuth coverage (degrees)
% 	98-103 	f6.2 	distance to closest station (degrees)
% 	105-110 	f6.2 	distance to furthest station (degrees)
% 	112 	a1 	analysis type: (a = automatic, m = manual, g = guess)
% 	114 	a1 	location method: (i = inversion, p = pattern recognition, g = ground truth, o = other)
% 	116-117 	a2 	event type:
% 			uk = unknown
% 			de = damaging earthquake ( Not standard IMS )
% 			fe = felt earthquake ( Not standard IMS )
% 			ke = known earthquake
% 			se = suspected earthquake
% 			kr = known rockburst
% 			sr = suspected rockburst
% 			ki = known induced event
% 			si = suspected induced event
% 			km = known mine expl.
% 			sm = suspected mine expl.
% 			kh = known chemical expl. ( Not standard IMS )
% 			sh = suspected chemical expl. ( Not standard IMS )
% 			kx = known experimental expl.
% 			sx = suspected experimental expl.
% 			kn = known nuclear expl.
% 			sn = suspected nuclear explosion
% 			ls = landslide
% 	119-127 	a9 	author of the origin
% 	129-136 	a8 	origin identification 

% extract fields
% - note we use str2double+cellstr to set empty numeric entries with nans
% - blank string fields are left as blank
isc.year=str2double(cellstr(txt0(:,1:4)));
isc.month=str2double(cellstr(txt0(:,6:7)));
isc.day=str2double(cellstr(txt0(:,9:10)));
isc.hour=str2double(cellstr(txt0(:,12:13)));
isc.minute=str2double(cellstr(txt0(:,15:16)));
isc.seconds=str2double(cellstr(txt0(:,18:22)));
isc.fixed_origin_time=cellstr(txt0(:,23));
isc.origin_time_error=str2double(cellstr(txt0(:,25:29)));
isc.rms_time_residuals=str2double(cellstr(txt0(:,31:35)));
isc.latitude=str2double(cellstr(txt0(:,37:44)));
isc.longitude=str2double(cellstr(txt0(:,46:54)));
isc.fixed_epicenter=cellstr(txt0(:,55));
isc.semimajor=str2double(cellstr(txt0(:,56:60)));
isc.semiminor=str2double(cellstr(txt0(:,62:66)));
isc.strike=str2double(cellstr(txt0(:,68:70)));
isc.depth=str2double(cellstr(txt0(:,72:76)));
isc.fixed_depth=cellstr(txt0(:,77));
isc.depth_error=str2double(cellstr(txt0(:,79:82)));
isc.nphase=str2double(cellstr(txt0(:,84:87)));
isc.nsta=str2double(cellstr(txt0(:,89:92)));
isc.gap=str2double(cellstr(txt0(:,94:96)));
isc.mindist=str2double(cellstr(txt0(:,98:103)));
isc.maxdist=str2double(cellstr(txt0(:,105:110)));
isc.analysis_type=cellstr(txt0(:,112));
isc.location_method=cellstr(txt0(:,114));
isc.event_type=cellstr(txt0(:,116:117));
isc.author=cellstr(txt0(:,119:127));
isc.origin_id=str2double(cellstr(txt0(:,129:136)));

end
