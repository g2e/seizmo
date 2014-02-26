function [isc]=parse_isc_fm(file,flag)
%PARSE_ISC_FM    Parses focal mechanism & moment tensor output from the ISC
%
%    Usage:    isc=parse_isc_fm(file)
%              isc=parse_isc_fm(string,true)
%
%    Description:
%     ISC=PARSE_ISC_FM(FILE) reads an ascii file containing focal mechanism
%     output from the ISC bulletin.  To create your own file go here:
%       http://www.isc.ac.uk/iscbulletin/search/fmechanisms/
%     and make sure the output format is set to CSV formatted focal
%     mechanisms.  Save the page as an .html file.  FILE should be the
%     filename (and path if necessary) of the saved .html file.  The output
%     of this function (ISC) is a scalar struct with fields for all of the
%     ISC origin datums.  Each field will be NFMECHx1 and will either be a
%     double array or a cellstr array and the function SSIDX is provided to
%     pull info from scalar structs.
%
%     ISC=PARSE_ISC_FM(STRING,TRUE) will take the first input to be a
%     character string of the contents of an ISC focal mechanism .html file
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
%     isc=parse_isc_fm();
%
%     % Do a test run:
%     link=[];
%     isc=parse_isc_fm(urlread(link),true)
%
%     % Extract the 3rd focal mechanism from
%     % the returned isc focal mechanism struct:
%     fm=ssidx(isc,3)
%
%    See also: PARSE_ISC_ORIGIN, READNDK, READSODEVENTCSV, READ_USGS_MT,
%              READ_USGS_FM, SSIDX

%     Version History:
%        Feb.  9, 2014 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Feb.  9, 2014 at 17:25 GMT

% todo:

% check nargin
error(nargchk(0,2,nargin));

% default flag
if(nargin<2 || isempty(flag)); flag=false; end

% skip if string input
if(~flag)
    % read in ndk file
    if(nargin<1); file=[]; end
    txt=readtxt(file,{'*.html;*.HTML' 'HTML Files (*.html,*.HTML)';
        '*.txt;*.TXT' 'TXT Files (*.txt,*.TXT)';
        '*.*' 'All Files (*.*)'});
else
    % just copy file to txt
    if(nargin<1 || isempty(file) || ~isstring(file))
        error('seizmo:parse_isc_fm:emptyStr',...
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

% replace header because we need to fix some names
fixhdr=['EVENT_ID,ORIGIN_AUTHOR,DATE,TIME,LATITUDE,LONGITUDE,DEPTH,' ...
    'CENTROID,FM_AUTHOR,EXPONENT,SCALARMOMENT,MW,MXX_EXPONENT,MRR,MTT,' ...
    'MPP,MRT,MTP,MRP,STRIKE1,DIP1,RAKE1,STRIKE2,DIP2,RAKE2,' ...
    'PA_EXPONENT,T_VAL,T_PL,T_AZM,P_VAL,P_PL,P_AZM,N_VAL,N_PL,N_AZM'];
lines{b+2}=fixhdr;

% lowercase field names
lines{b+2}=lower(lines{b+2});

% parse comma separated values portion using readcsv
isc=readcsv(joinwords(lines(b+2:e-1),sprintf('\n')),[],true);

% convert some fields from cellstr to something more useful
isc.event_id=str2double(isc.event_id);
%isc.date
date=char(isc.date);
isc.year=str2double(cellstr(date(:,1:4)));
isc.month=str2double(cellstr(date(:,6:7)));
isc.day=str2double(cellstr(date(:,9:10)));
%isc.time
time=char(isc.time);
isc.hour=str2double(cellstr(time(:,1:2)));
isc.minute=str2double(cellstr(time(:,4:5)));
isc.seconds=str2double(cellstr(time(:,7:11)));
isc.latitude=str2double(isc.latitude);
isc.longitude=str2double(isc.longitude);
isc.depth=str2double(isc.depth);
isc.centroid=strcmp('TRUE',isc.centroid);
isc.exponent=str2double(isc.exponent);
isc.scalarmoment=str2double(isc.scalarmoment);
isc.mw=str2double(isc.mw);
isc.mxx_exponent=str2double(isc.mxx_exponent);
isc.mrr=str2double(isc.mrr);
isc.mtt=str2double(isc.mtt);
isc.mpp=str2double(isc.mpp);
isc.mrt=str2double(isc.mrt);
isc.mtp=str2double(isc.mtp);
isc.mrp=str2double(isc.mrp);
isc.strike1=str2double(isc.strike1);
isc.dip1=str2double(isc.dip1);
isc.rake1=str2double(isc.rake1);
isc.strike2=str2double(isc.strike2);
isc.dip2=str2double(isc.dip2);
isc.rake2=str2double(isc.rake2);
isc.pa_exponent=str2double(isc.pa_exponent);
isc.t_val=str2double(isc.t_val);
isc.t_pl=str2double(isc.t_pl);
isc.t_azm=str2double(isc.t_azm);
isc.p_val=str2double(isc.p_val);
isc.p_pl=str2double(isc.p_pl);
isc.p_azm=str2double(isc.p_azm);
isc.n_val=str2double(isc.n_val);
isc.n_pl=str2double(isc.n_pl);
isc.n_azm=str2double(isc.n_azm);

end
