function [slowness]=add_2sta_slowness_table(db)
%ADD_2STA_SLOWNESS_TABLE    Adds 2-station slowness table to the db
%
%    Usage:
%
%    Description:
%
%    Notes:
%
%    Examples:
%
%    See also:

%     Version History:
%        Aug. 25, 2009 - initial version
%
%     Written by Garrett Euler (ggeuler at wustl dot edu)
%     Last Updated Aug. 25, 2009 at 21:15 GMT

% todo:
% - output table
% - dispersion plots
%   - bank comparison
%   - corrected vs uncorrected
%   - phase vs group
%   - 2sta for all event

% check nargin
msg=nargchk(1,1,nargin);
if(~isempty(msg)); error(msg); end

% check db
if(~isstruct(db))
    error('seizmo:add_2sta_slowness_table:badInput',...
        'Input DB must be a struct!');
end
fields=fieldnames(db);
reqfields={'event' 'phase' 'bank' 'band' 'station' 'align' 'correct'};
if(any(~ismember(reqfields,fields)))
    error('seizmo:add_2sta_slowness_table:badInput',...
        ['Input DB must have the following fields:\n'...
        sprintf('%s ',reqfields{:})]);
end

% ok we need to loop over event, phase, filter bank and band, and stations
for i=1:size(db.event,1)
    % lets get event to station geometry first
    [gcarc,az,baz]=sphericalinv(db.event(i,7),db.event(i,8),db.station(:,2),db.station(:,3));
    
    % extract align records for this event
    evrec=db.align(db.align(:,1)==db.event(i,1),:);
    
    % jump if empty
    if(isempty(evrec)); continue; end
    
    % loop over phase
    for j=1:size(db.phase,1)
        % extract align records for this event and phase
        phrec=evrec(evrec(:,5)==db.phase(j,1),:);
        
        % jump if empty
        if(isempty(phrec)); continue; end
        
        % loop over station1
        for k=1:size(db.station,1)
            % extract align records for this event, phase, sta1
            strec1=phrec(phrec(:,2)==db.station(k,1));
            
            % jump if empty
            if(isempty(strec1)); continue; end
            
            % get travel time corrections
            cor1=db.correct(db.correct(:,1)==db.event(i,1) ...
                & db.correct(:,2)==db.phase(j,1) ...
                & db.correct(:,3)==db.station(k,1),:);
            [crust1,elev1,ellip1,mantle1]=deal(cor1(4),cor1(5),cor1(6),cor1(7));
            
            % loop over station2
            for l=k+1:size(db.station,1)
                % extract align records for this event, phase, sta2
                strec2=phrec(phrec(:,2)==db.station(l,1),:);
                
                % jump if empty
                if(isempty(strec2)); continue; end
                
                % get travel time corrections
                cor2=db.correct(db.correct(:,1)==db.event(i,1) ...
                    & db.correct(:,2)==db.phase(j,1) ...
                    & db.correct(:,3)==db.station(l,1),:);
                [crust2,elev2,ellip2,mantle2]=deal(cor2(4),cor2(5),cor2(6),cor2(7));
                
                % loop over bank
                h=[]; % set empty figure handle
                for m=1:size(db.bank,1)
                    % loop over band
                    c=1; % initialize counter
                    for n=find(db.band(:,1)'==db.bank(m,1))
                        % extract align records for this event, phase, bank, band
                        bdrec1=strec1(strec1(:,3)==db.band(n,1) & strec1(:,4)==db.band(n,2),:);
                        bdrec2=strec2(strec2(:,3)==db.band(n,1) & strec2(:,4)==db.band(n,2),:);
                        
                        % jump if not 2 recs
                        if(isempty(bdrec1) || isempty(bdrec2)); continue; end
                        
                        % get uncorrected slowness and error
                        [ucslow,uccovslow]=wlinem(gcarc([k l]),...
                            [bdrec1(11) bdrec2(11)],...
                            diag([bdrec1(12) bdrec2(12)].^2),...
                            diag([bdrec1(10) bdrec2(10)]));
                        ucslowerr(c)=2*sqrt(uccovslow(2,2));
                        ucdisp(c)=ucslow(2);
                        
                        % get corrected travel times
                        ctt1=bdrec1(11)-crust1-elev1-ellip1-mantle1;
                        ctt2=bdrec2(11)-crust2-elev2-ellip2-mantle2;
                        
                        % get slowness and error
                        [cslow,ccovslow]=wlinem(gcarc([k l]),[ctt1 ctt2],...
                            diag([bdrec1(12) bdrec2(12)].^2),...
                            diag([bdrec1(10) bdrec2(10)]));
                        cslowerr(c)=2*sqrt(ccovslow(2,2));
                        cdisp(c)=cslow(2);
                        c=c+1;
                    end
                    if(c==1); continue; end
                    if(isempty(h)); h=figure; end
                    plot(disp)
                    hold on
                end
                if(isempty(h)); continue; end
                hold off
                die
            end
        end
    end
end

end
