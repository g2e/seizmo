function [data]=extract_plot_data(h)
%EXTRACT_PLOT_DATA    Extracts SEIZMO data struct from SEIZMO plot
%
%    Usage:    data=extract_plot_data(handle)
%
%    Description:
%     DATA=EXTRACT_PLOT_DATA(HANDLE) returns the SEIZMO dataset plotted in
%     the axes/figure given by HANDLE.  HANDLE must either be a single
%     figure handle or axes handle (multiple axes handles are allowed but
%     should come from the same plotting call).
%
%    Notes:
%     - for a multi-dataset PLOT2 figure data is returned as a cell array
%       with each dataset in a separate cell
%     - does not work with spectral records
%
%    Examples:
%     % Plot, extract, verify:
%     plot2(data);
%     data2=extract_plot_data(gcf);
%     isequal(data,data2)
%
%    See also: PLOT0, PLOT1, PLOT2, RECORDSECTION

%     Version History:
%        Apr. 19, 2011 - initial version
%        Mar. 13, 2012 - use getheader improvements, code reduction
%        Sep. 28, 2012 - fixed nasty eps usage bugs
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Sep. 28, 2012 at 23:00 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check handle
if(~isreal(h) || any(~ishandle(h)))
    error('seizmo:extract_plot_data:badInput',...
        'HANDLE must be a valid figure, axes or line handle!');
elseif(~isscalar(h) && any(strcmpi(get(h,'type'),'figure')))
    error('seizmo:extract_plot_data:badInput',...
        'Figure handle input must be scalar!');
end

% find records
rh=findobj(h,'tag','record');

% error if none found
if(isempty(rh))
    error('seizmo:extract_plot_data:badInput',...
        'No SEIZMO records found under the given graphics handle(s)!');
end

% extract metadata
data=get(rh,'userdata');
if(iscell(data)); data=[data{:}]'; end

% extract data
leven=getheader(data,'leven lgc');
for i=1:numel(data)
    switch leven{i}
        case 'true'  % evenly spaced
            data(i).dep=get(rh(i),'ydata').';
        case 'false' % unevenly spaced
            data(i).ind=get(rh(i),'xdata').';
            data(i).dep=get(rh(i),'ydata').';
        otherwise    % fail
            error('seizmo:extract_plot_data:corruptHeader',...
                'Record with handle %d has a bad LEVEN field!',rh(i));
    end
end

% sort by indexing
idx=cat(1,data.index);
[idx,didx]=sortrows(idx);
data=data(didx);

% break if any duplicates in idx
% as this indicates that records
% were plotted using multiple calls
% which we do not handle!
if(~isequal(unique(idx,'rows'),idx))
    error('seizmo:extract_plot_data:badInput',...
        'Plotted records must come from a single plotting call!');
end

% multiple datasets (plot2 only)
if(size(idx,2)==3)
    % 1st column indices
    [a,c,c]=unique(idx(:,1));
    cnt=nan(size(a));
    for i=1:numel(a); cnt(i)=sum(c==a(i)); end
    idx=mat2cell(idx(:,2:3),cnt);
    data=mat2cell(data,cnt);
else
    data={data};
    idx={idx};
end

% loop over datasets
for i=1:numel(data)
    % combine components
    if(any(idx{i}(:,2)>1))
        % get record indices
        [a,c,c]=unique(idx{i}(:,1));
        
        % find multicomponent
        delete=false(numel(data{i}),1);
        for j=1:numel(a)
            if(sum(a(j)==c)>1)
                % combine .dep
                midx=find(a(j)==c);
                data{i}(midx(1)).dep=[data{i}(midx(1:end)).dep];
                
                % set delete bits on
                delete(midx(2:end))=true;
            end
        end
        data{i}(delete)=[];
    end
    
    % remove index field
    data{i}=rmfield(data{i},'index');
    
    % fixing amplitudes (due to normalization/offset)
    [b,npts,delta,z6,leven,depmin,depmax]=getheader(data{i},...
        'b','npts','delta','z6','leven lgc','depmin','depmax');
    leven=~strcmpi(leven,'false');
    z6=datenum(cell2mat(z6));
    for j=1:numel(data{i})
        % fixing amplitudes (due to normalization/offset)
        dmax=max(data{i}(j).dep(:));
        dmin=min(data{i}(j).dep(:));
        if(dmax~=depmax(j) || dmin~=depmin(j))
            % rescale to original limits
            data{i}(j).dep=(depmax(j)-depmin(j))/(dmax-dmin)...
                *data{i}(j).dep;
            data{i}(j).dep=data{i}(j).dep+depmin(j)-min(data{i}(j).dep(:));
        end
        
        % fixing timing of uneven records plotted in abs time
        if(~leven && (abs(b(j)-data{i}(j).ind(1))>eps(single(10*b(j))) ...
                || ((data{i}(j).ind(end)-data{i}(j).ind(1))...
                /(npts(j)-1)-delta(j)>eps(single(10*delta(j))))))
            data{i}(j).ind=(data{i}(j).ind-z6(j))*86400;
        end
    end
end

% uncell single dataset
if(isscalar(data)); data=data{1}; end

end
