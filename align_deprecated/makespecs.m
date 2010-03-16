function [f,rows,cols,rspecs,wspecs]=makespecs(s)
%MAKESPECS    Makes format specs for fields in a structure
%
%    Description:  Returns info on each field in a structure so that they
%     can be written by wconf and subsequently read by rconf.  Specifically
%     writing is assumed to be done with fprintf and reading with textscan.
%
%    Notes:
%     - Substructures will throw a warning and will be deleted from the
%       output (they need to be delt with separately)
%
%    Usage: [fieldnames,nrows,ncols,readspecs,writespecs]=makespecs(CONF)
%
%    Examples:
%
%    See also: wconf, rconf

% get fields
f=fieldnames(s);

% allocation
nf=length(f);
rows=zeros(nf,1);
cols=rows;
rspecs=cell(nf,1);
wspecs=cell(nf,1);
destroy=false(nf,1);

% loop through fields
for i=1:nf
    % get size
    [rows(i),cols(i)]=size(s.(f{i}));
    % special handling for empty fields
    if(rows(i)*cols(i)==0)
        rows(i)=0;
    % numeric or logical
    elseif(isnumeric(s.(f{i})) || islogical(s.(f{i})))
        rspecs{i}=repmat('%f',1,cols(i));
        wspecs{i}=repmat(' %f',1,cols(i));
    % char
    elseif(ischar(s.(f{i})))
        rspecs{i}=['%' num2str(cols(i)) 'c'];
        wspecs{i}=[' %' num2str(cols(i)) 's'];
        cols(i)=1; % just one column
    % cell/cellstr
    elseif(iscellstr(s.(f{i})) || iscell(s.(f{i})))
        rspecs{i}=repmat('%s',1,cols(i));
        wspecs{i}=repmat('%s',1,cols(i));
    % otherwise warn and delete (substructure)
    else
        warning('SAClab:makespecs:badFieldType',...
            'Field %s is an unsupported output type',f{i})
        destroy(i)=true;
    end
end

% destroy bad fields
f(destroy)=[];
rows(destroy)=[];
cols(destroy)=[];
rspecs(destroy)=[];
wspecs(destroy)=[];

end
