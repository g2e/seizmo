function [varargout]=coverage(data)
%COVERAGE    Convert SEIZMO data records to coverage
%
%    Usage:    spans=coverage(data)
%              [spans,gaps,overlaps]=coverage(data)
%              coverage(data)
%
%    Description:
%     SPANS=COVERAGE(DATA) converts SEIZMO records in DATA to just contain
%     the start and end times for that particular component.  Components
%     are grouped based on the KNAME group field.  This is for making
%     coverage plots (see Examples section below).
%
%     [SPANS,GAPS,OVERLAPS]=COVERAGE(DATA) also returns gaps and overlaps
%     start and end times using the SEIZMO record format.  See the Examples
%     section below to show how display this information.
%
%     COVERAGE(DATA) with no outputs will plot the spans, gaps & overlaps
%     in 3 figures like in the Examples section below.
%
%    Notes:
%
%    Examples:
%     % Look at coverage, gaps and overlaps:
%     [c,g,o]=coverage(data);
%     plot0(c,'abs',true,'namesonyaxis',true,'recwidth',2,'title','spans');
%     plot0(g,'abs',true,'namesonyaxis',true,'recwidth',2,'title','gaps');
%     plot0(o,'abs',true,'namesonyaxis',true,'recwidth',2,...
%             'title','overlaps');
%
%    See also: PLOT0, MELD

%     Version History:
%        Mar.  7, 2012 - initial version
%        July 25, 2014 - fix no overlaps bug, fix spectral records
%                        requiring data, bugfix: data/delta needed
%                        reordering
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated July 25, 2014 at 23:00 GMT

% todo:

% check nargin
error(nargchk(1,1,nargin));

% check struct
error(seizmocheck(data));

% detect spectral
iftype=getheader(data,'iftype id');
spec=strcmpi(iftype,'irlim') | strcmpi(iftype,'iamph');

% group by kname
[idx,componentcode]=getcomponentidx(data);
ncmp=max(idx);

% get necessary header fields
[butc,eutc,z,delta,sb,nspts,sdelta]=getheader(data,'b utc','e utc','z',...
    'delta','sb utc','nspts','sdelta');
butc=cell2mat(butc); eutc=cell2mat(eutc); z=cell2mat(z); sb=cell2mat(sb);

% fix b/e for spectral
if(any(spec))
    delta(spec)=sdelta(spec);
    butc(spec,:)=sb(spec,:);
    se=sb;
    se(:,end)=se(:,end)+(nspts-1).*sdelta;
    eutc(spec,:)=fixtimes(se(spec,:),'utc');
end

% combine start/end times into one record
[keep,npts]=deal(nan(ncmp,1));
for i=1:ncmp
    % keep/use first record
    cmp=idx==i;
    j=find(cmp,1);
    keep(i)=j;
    
    % create coverage record for this component
    data(j).ind=[timediff(z(j,:),butc(cmp,:)) ...
        timediff(z(j,:),eutc(cmp,:)) nan(sum(cmp),1)];
    data(j).ind=sortrows(data(j).ind).';
    data(j).ind=data(j).ind(:);
    data(j).dep=zeros(size(data(j).ind));
    npts(i)=numel(data(j).ind);
    
    % delta check (set delta to 0 for variable sample rate channels)
    if(~isscalar(unique(delta(cmp)))); delta(j)=0; end
end

% reorder & remove unneeded
data=data(keep);
delta=delta(keep);

% adjust fields as needed
[data.hasdata]=deal(true);
[data.name]=deal(componentcode{:});
data=changeheader(data,'npts',npts,'leven',false,'delta',delta,'dep',0);
if(nargout==1); varargout={data}; return; end

% get gaps & overlaps
gap=data;
if(any(nargout==[0 3])); overlap=data; end
gnpts=npts; onpts=npts;
for i=1:ncmp
    % create gap record for this component
    b=data(i).ind(1:3:end); e=data(i).ind(2:3:end);
    gap(i).ind=[]; gap(i).dep=[];
    for j=1:numel(e)-1
        % is this end time within some other record?
        if(any(e(j)+delta(i)>=b & e(j)<e)) % yes
            % do nothing since there is no data gap here
        else % nope, so there is possibly a gap
            % this includes duplicates or those ending at the same time
            % so we need to remove all but the first of those for gaps
            if(any(e(j)==e(1:j-1))); continue; end
            
            % first begin time after current end time
            gt=b-e(j);
            gt=min(gt(gt>0));
            if(isempty(gt)); continue; end
            
            % gap found, so add it in
            gap(i).ind=[gap(i).ind; e(j)+delta(i) e(j)+gt nan];
        end
    end
    
    % clean up
    gap(i).ind=sortrows(gap(i).ind).';
    gap(i).ind=gap(i).ind(:);
    gap(i).dep=zeros(size(gap(i).ind));
    gnpts(i)=numel(gap(i).ind);
    if(nargout==2); continue; end
    
    % create overlap record for this component
    overlap(i).ind=[]; overlap(i).dep=[];
    for j=1:numel(b)-1
        for k=j+1:numel(b)
            % get overlap
            if(b(j)==b(k) && e(j)==e(k))
                % exact duplicate
                overlap(i).ind=[overlap(i).ind; b(j) e(j)+delta(i) nan];
            elseif(b(j)==b(k))
                % duplicate start
                overlap(i).ind=[overlap(i).ind; ...
                    b(j) min(e(j),e(k))+delta(i) nan];
            elseif(e(j)==e(k))
                % duplicate end
                overlap(i).ind=[overlap(i).ind; ...
                    max(b(j),b(k)) e(j)+delta(i) nan];
            elseif(b(j)>b(k) && e(j)<e(k))
                % within
                overlap(i).ind=[overlap(i).ind; b(j) e(j)+delta(i) nan];
            elseif(b(j)<b(k) && e(j)>e(k))
                % surround
                overlap(i).ind=[overlap(i).ind; b(k) e(k)+delta(i) nan];
            elseif(b(j)>b(k) && b(j)<e(k)+delta(i))
                % front overlap
                overlap(i).ind=[overlap(i).ind; b(j) e(k)+delta(i) nan];
            elseif(e(j)<e(k) && e(j)+delta(i)>b(k))
                % end overlap
                overlap(i).ind=[overlap(i).ind; b(k) e(j)+delta(i) nan];
            end
        end
    end
    
    % reduce overlaps to minimum set
    nov=size(overlap(i).ind,1);
    if(nov>1)
        b=overlap(i).ind(:,1);
        e=overlap(i).ind(:,2);
        j=1;
        keep=true(nov,1);
        while j<=nov
            for k=j+1:nov
                % skip deleted
                if(~keep); continue; end
                
                % compare overlap timings
                if(b(j)>=b(k) && e(j)<=e(k))
                    % within another so delete me and go to next overlap
                    keep(j)=false;
                    break;
                elseif(b(j)<=b(k) && e(j)>=e(k))
                    % surrounds another so delete them
                    keep(k)=false;
                elseif(b(j)>b(k) && b(j)<e(k))
                    % extend b, delete them, redo
                    overlap(i).ind(j,1)=overlap(i).ind(k,1);
                    keep(k)=false;
                    j=j-1; % balance increment below
                    break;
                elseif(e(j)<e(k) && e(j)>b(k))
                    % extend e, delete them, redo
                    overlap(i).ind(j,2)=overlap(i).ind(k,2);
                    keep(k)=false;
                    j=j-1; % balance increment below
                    break;
                end
            end
            j=j+1;
        end
        overlap(i).ind(~keep,:)=[];
    end
    
    % clean up
    overlap(i).ind=sortrows(overlap(i).ind).';
    overlap(i).ind=overlap(i).ind(:);
    overlap(i).dep=zeros(size(overlap(i).ind));
    onpts(i)=numel(overlap(i).ind);
end

% fix headers and return
gap=changeheader(gap,'npts',gnpts);
if(nargout==2); varargout={data gap}; return; end
overlap=changeheader(overlap,'npts',onpts);
if(nargout==3); varargout={data gap overlap}; return; end

% no outputs so plot results
plot0(data,'abs',true,'namesonyaxis',true,'recwidth',2,'title','spans');
plot0(gap,'abs',true,'namesonyaxis',true,'recwidth',2,'title','gaps');
plot0(overlap,'abs',true,'namesonyaxis',true,'recwidth',2,...
        'title','overlaps');

end
