function [events]=readndk(filename)
%READNDK    Reads a .ndk formatted text file into a struct
%
%    Description:
%     
%     Layout:
%      

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end;

% check filename is char
if(~ischar(filename)); error('FILENAME must be a string!'); end

% open file
fid=fopen(filename);

% check if file exists
if(fid<0); error('%s not openable!',filename); end

% just read in plain as matlab licks the nutsack when it comes to parsing
% strict formats (stop fucking deleting whitespace before/after asshats)
a=textscan(fid,'%s','delimiter','\n','whitespace',''); a=a{1}; a=char(a);

% close file
fclose(fid);

% allocate struct
n=size(a,1)/5;
events(1:n,1)=struct();

% extract each field separately (cause matlab sure as fuck will not)
% wrap them in cell arrays so we can easily push into a struct
temp=cellstr(a(1:5:end,1:4));
[events.catalog]=deal(temp{:});
temp=num2cell(str2num(a(1:5:end,6:9)));
[events.year]=deal(temp{:});
temp=num2cell(str2num(a(1:5:end,11:12)));
[events.month]=deal(temp{:});
temp=num2cell(str2num(a(1:5:end,14:15)));
[events.day]=deal(temp{:});
temp=num2cell(str2num(a(1:5:end,17:18)));
[events.hour]=deal(temp{:});
temp=num2cell(str2num(a(1:5:end,20:21)));
[events.minute]=deal(temp{:});
temp=num2cell(str2num(a(1:5:end,23:26)));
[events.seconds]=deal(temp{:});
temp=num2cell(str2num(a(1:5:end,28:33)));
[events.latitude]=deal(temp{:});
temp=num2cell(str2num(a(1:5:end,35:41)));
[events.longitude]=deal(temp{:});
temp=num2cell(str2num(a(1:5:end,43:47)));
[events.depth]=deal(temp{:});
temp=num2cell(str2num(a(1:5:end,49:51)));
[events.mb]=deal(temp{:});
temp=num2cell(str2num(a(1:5:end,53:55)));
[events.ms]=deal(temp{:});
temp=cellstr(a(1:5:end,57:80));
[events.location]=deal(temp{:});

temp=cellstr(a(2:5:end,1:16));
[events.name]=deal(temp{:});
temp=cellstr(a(2:5:end,18));
[events.data1type]=deal(temp{:});
temp=num2cell(str2num(a(2:5:end,20:22)));
[events.data1nstn]=deal(temp{:});
temp=num2cell(str2num(a(2:5:end,23:27)));
[events.data1ncmp]=deal(temp{:});
temp=num2cell(str2num(a(2:5:end,28:31)));
[events.data1sper]=deal(temp{:});
temp=cellstr(a(2:5:end,33));
[events.data2type]=deal(temp{:});
temp=num2cell(str2num(a(2:5:end,35:37)));
[events.data2nstn]=deal(temp{:});
temp=num2cell(str2num(a(2:5:end,38:42)));
[events.data2ncmp]=deal(temp{:});
temp=num2cell(str2num(a(2:5:end,43:46)));
[events.data2sper]=deal(temp{:});
temp=cellstr(a(2:5:end,48));
[events.data3type]=deal(temp{:});
temp=num2cell(str2num(a(2:5:end,50:52)));
[events.data3nstn]=deal(temp{:});
temp=num2cell(str2num(a(2:5:end,53:57)));
[events.data3ncmp]=deal(temp{:});
temp=num2cell(str2num(a(2:5:end,58:61)));
[events.data3sper]=deal(temp{:});
temp=cellstr(a(2:5:end,63:68));
[events.srcinvtype]=deal(temp{:});
temp=cellstr(a(2:5:end,70:74));
[events.srcfunctype]=deal(temp{:});
temp=num2cell(str2num(a(2:5:end,76:80)));
[events.srcfuncdur]=deal(temp{:});

temp=cellstr(a(3:5:end,1:8));
[events.centroidstring]=deal(temp{:});
temp=num2cell(str2num(a(3:5:end,13:18)));
[events.centroidtime]=deal(temp{:});
temp=num2cell(str2num(a(3:5:end,20:22)));
[events.centroidtimeerr]=deal(temp{:});
temp=num2cell(str2num(a(3:5:end,24:29)));
[events.centroidlat]=deal(temp{:});
temp=num2cell(str2num(a(3:5:end,31:34)));
[events.centroidlaterr]=deal(temp{:});
temp=num2cell(str2num(a(3:5:end,36:42)));
[events.centroidlon]=deal(temp{:});
temp=num2cell(str2num(a(3:5:end,44:47)));
[events.centroidlonerr]=deal(temp{:});
temp=num2cell(str2num(a(3:5:end,49:53)));
[events.centroiddep]=deal(temp{:});
temp=num2cell(str2num(a(3:5:end,55:58)));
[events.centroiddeperr]=deal(temp{:});
temp=cellstr(a(3:5:end,60:63));
[events.depthtype]=deal(temp{:});
temp=cellstr(a(3:5:end,65:80));
[events.timestamp]=deal(temp{:});

temp=num2cell(str2num(a(4:5:end,1:2)));
[events.exponent]=deal(temp{:});
temp=num2cell(str2num(a(4:5:end,4:9)));
[events.mrr]=deal(temp{:});
temp=num2cell(str2num(a(4:5:end,11:15)));
[events.mrrerr]=deal(temp{:});
temp=num2cell(str2num(a(4:5:end,17:22)));
[events.mtt]=deal(temp{:});
temp=num2cell(str2num(a(4:5:end,24:28)));
[events.mtterr]=deal(temp{:});
temp=num2cell(str2num(a(4:5:end,30:35)));
[events.mpp]=deal(temp{:});
temp=num2cell(str2num(a(4:5:end,37:41)));
[events.mpperr]=deal(temp{:});
temp=num2cell(str2num(a(4:5:end,43:48)));
[events.mrt]=deal(temp{:});
temp=num2cell(str2num(a(4:5:end,50:54)));
[events.mrterr]=deal(temp{:});
temp=num2cell(str2num(a(4:5:end,56:61)));
[events.mrp]=deal(temp{:});
temp=num2cell(str2num(a(4:5:end,63:67)));
[events.mrperr]=deal(temp{:});
temp=num2cell(str2num(a(4:5:end,69:74)));
[events.mtp]=deal(temp{:});
temp=num2cell(str2num(a(4:5:end,76:80)));
[events.mtperr]=deal(temp{:});

temp=cellstr(a(5:5:end,1:3));
[events.version]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,5:11)));
[events.eigval1]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,13:14)));
[events.plunge1]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,16:18)));
[events.azimuth1]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,20:26)));
[events.eigval2]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,28:29)));
[events.plunge2]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,31:33)));
[events.azimuth2]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,35:41)));
[events.eigval3]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,43:44)));
[events.plunge3]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,46:48)));
[events.azimuth3]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,50:56)));
[events.scalarmoment]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,58:60)));
[events.strike1]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,62:63)));
[events.dip1]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,65:68)));
[events.rake1]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,70:72)));
[events.strik2]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,74:75)));
[events.dip2]=deal(temp{:});
temp=num2cell(str2num(a(5:5:end,77:80)));
[events.rake2]=deal(temp{:});

end
