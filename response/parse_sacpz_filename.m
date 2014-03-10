function [gd,knetwk,kstnm,kcmpnm,khole,b,e]=parse_sacpz_filename(filename)
%PARSE_SACPZ_FILENAME    Parses info in the SAC PoleZero filename
%
%    Usage:    [good,knetwk,kstnm,kcmpnm,khole,b,e]=...
%                  parse_sacpz_filename(filenames)
%
%    Description:
%     [GOOD,KNETWK,KSTNM,KCMPNM,KHOLE,B,E]=PARSE_SACPZ_FILENAME(FILENAMES)
%     parses SAC Polezero filenames for critical info on which records
%     and the time frame the PoleZero is valid.  Supported formats are
%     detailed in the notes section below.  FILENAMES must be a char or
%     cellstr array of filenames.  GOOD is a logical array indicating which
%     filenames were valid.  KNETWK, KSTNM, KCMPNM, KHOLE are column vector
%     cellstr arrays defining the network, station, component, and stream
%     names associated with the PoleZero files.  B and E are Nx5 numeric
%     arrays of [yr jday hr mn secs] giving the time range the PoleZeros
%     are valid.  All returned variables except GOOD are reduced (the
%     invalid file info is removed).
%
%    Notes:
%     - Supported SAC Polezero filename formats:
%        SAC_PZs_NT_STA_CMP_LL_YYYY.DDD.HH.MM.SS.FFF_YYYY.DDD.HH.MM.SS.FFF
%        SAC_PZs_NT_STA_CMP_LL_YYYY.DDD.HH.MM.SS.FFF
%       where
%        NT   = Network Name (ie IU, XB, etc) - at least 1 char
%        STA  = Station Name (ie T033, CMB, etc) - at least 1 char
%        CMP  = Component Name (ie BHZ, LH1, etc) - at least 1 char
%        LL   = Stream Name (ie __, 01, 02) - 0-2 chars
%        YYYY = Year - 1+ digits
%        DDD  = Julian day - 1+ digits
%        HH   = Hour - 1+ digits
%        MM   = Minute - 1+ digits
%        SS   = Second - 1+ digits
%        FFF  = Fractional seconds - 1+ digits
%
%    Examples:
%     % PARSE_SACPZ_FILENAME not only parses a list of filenames but also
%     % filters out the invalid names.  This reduces the burden of having
%     % to eliminate non-SAC PoleZero files from the file list.  For
%     % instance, to get all the filename info for valid SAC PoleZero files
%     % in the current directory:
%     files=dir('*');
%     filenames={files.name}.';
%     [good,knetwk,kstnm,kcmpnm,khole,b,e]=...
%         parse_sacpz_filename(filenames)
%
%    See also: READSACPZ, WRITESACPZ, GENSACPZNAME

%     Version History:
%        Sep. 20, 2009 - initial version
%        May  27, 2010 - seizmoverbose support
%        Feb. 11, 2011 - mass nargchk fix
%        Mar.  5, 2011 - improved empty khole detection
%        Feb.  3, 2012 - doc update
%        Mar. 10, 2014 - update See also section
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Mar. 10, 2014 at 15:05 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check filename
sz=size(filename);
if(numel(sz)>2)
    error('seizmo:parse_sacpz_filename:badInput',...
        'FILENAME must be 2D array!');
elseif(ischar(filename))
    filename=strtrim(cellstr(filename));
elseif(~iscellstr(filename))
    error('seizmo:parse_sacpz_filename:badInput',...
        'FILENAME must be a char or cellstr array!');
end

% verbosity
verbose=seizmoverbose;
if(verbose)
    disp('Parsing SAC Polezero Filename(s)');
end

% allocating output
npz=numel(filename);
knetwk=cell(npz,1); kstnm=knetwk; kcmpnm=knetwk; khole=knetwk;
b=[2999 365 23 59 59.99999]; b=b(ones(npz,1),:); e=b; gd=false(npz,1);

% skip last step if empty
if(npz==0); return; end

% loop over filenames
for i=1:npz
    % detail message
    if(verbose); print_time_left(i-1,npz); end
    
    % parse filename
    w=getwords(filename{i},'_');
    nw=numel(w);
    
    % require 6 or 7 or 8 fields
    if(~any(nw==[6 7 8])); continue; end
    
    % assign name fields
    [knetwk{i},kstnm{i},kcmpnm{i}]=deal(w{3:5});
    
    % check parsed name fields
    if(~strcmpi(w{1},'sac')); continue; end
    if(~strcmpi(w{2},'pzs')); continue; end
    if(isempty(knetwk{i})); continue; end
    if(isempty(kstnm{i})); continue; end
    if(isempty(kcmpnm{i})); continue; end
    
    % handle blank khole
    khole{i}=w{6}; skip=0;
    if(numel(khole{i})>5 && any(khole{i}(5:6)=='.')); skip=1; khole{i}=''; end
    
    % parse/check/convert begin time
    t=getwords(w{7-skip},'.');
    if(numel(t)~=6); continue; end
    t{5}=[t{5} '.' t{6}]; t(6)=[];
    b(i,:)=str2double(t);
    if(any(isnan(b(i,:)))); continue; end
    
    % parse/check/convert end time (if exists)
    if(nw==(8-skip))
        t=getwords(w{8-skip},'.');
        if(numel(t)~=6); continue; end
        t{5}=[t{5} '.' t{6}]; t(6)=[];
        e(i,:)=str2double(t);
        if(any(isnan(b(i,:)))); continue; end
    end
    gd(i)=true;
end

% detail message
if(verbose); print_time_left(i,npz); end

% remove unworthy
knetwk=knetwk(gd);
kstnm=kstnm(gd);
kcmpnm=kcmpnm(gd);
khole=khole(gd);
b=b(gd,:);
e=e(gd,:);

end
