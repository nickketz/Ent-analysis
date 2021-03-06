function trl = ent_trialfun_mff(cfg)

% operates using Net Station evt files and event structs

%cfg.eventinfo.offsetMSTCP = 0; %default offset if no DIN found

% convert single string into cell-array, otherwise intersection does not
% work as intended
if ischar(cfg.trialdef.eventvalue)
    cfg.trialdef.eventvalue = {cfg.trialdef.eventvalue};
end

% get the header and event information
fprintf('Reading flags from EEG file using FieldTrip...');
ft_hdr = ft_read_header(cfg.dataset);
[pathstr,name] = fileparts(cfg.dataset);
ftEventsFile = fullfile(pathstr,sprintf('%s_ftEvents.mat',name));
if exist(ftEventsFile,'file')
    ft_event = load(ftEventsFile);
    if isfield(ft_event,'date_string')
        warning('Using pre-saved FT events from this date: %s!',ft_event.date_string);
    else
        warning('Using pre-saved FT events from an unknown date!');
    end
    ft_event = ft_event.ft_event;
else
    tic
    ft_event = ft_read_event(cfg.dataset);
    toc
    date_string = datestr(now);
    fprintf('Saving FT events from MFF (current time: %s): %s...',date_string,ftEventsFile);
    save(ftEventsFile,'ft_event','date_string');
    fprintf('Done.\n');
end
fprintf('Done.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in external data, if wanted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if cfg.eventinfo.useMetadata
    md = cfg.eventinfo.metadata;
    
    if ismember('eventStruct', md.types)
        % read the events file
        eventsFile = fullfile(md.dirs.dataroot,md.dirs.behDir,md.subject,'events','events.mat');
        if exist(eventsFile,'file')
            fprintf('Loading events file: %s...',eventsFile);
            events_all = load(eventsFile,'events');
            events_all = events_all.events;
            fprintf('Done.\n');
        else
            error('Cannot find events file: %s\n',eventsFile)
        end
    end
    
    if ismember('expParam', md.types)
        % read the experiment parameters file
        expParamFile = fullfile(md.dirs.dataroot,md.dirs.behDir,md.subject,'experimentParams.mat');
        if exist(expParamFile,'file')
            fprintf('Loading experiment parameters file: %s...',expParamFile);
            load(expParamFile,'expParam');
            fprintf('Done.\n');
        else
            error('Cannot find experiment parameters file: %s\n',expParamFile)
        end
    end
end

if cfg.eventinfo.usePhotodiodeDIN
    triggers = {'FLK2', cfg.eventinfo.photodiodeDIN_str};
else
    triggers = {'FLK2'};
end

% initialize the trl matrix
trl = [];

% all trls need to have the same length
maxTrlCols = -Inf;
fn_trl_ord = fieldnames(cfg.eventinfo.trl_order);
for fn = 1:length(fn_trl_ord)
    if ismember(fn_trl_ord{fn},cfg.eventinfo.eventValues)
        if length(cfg.eventinfo.trl_order.(fn_trl_ord{fn})) > maxTrlCols
            maxTrlCols = length(cfg.eventinfo.trl_order.(fn_trl_ord{fn}));
        end
    end
end
if maxTrlCols == -Inf
    fprintf('Did not set maximum number of trialinfo columns!\n');
    keyboard
end
timeCols = 3;
eventNumCols = 1;
trl_ini = -1 * ones(1, timeCols + eventNumCols + maxTrlCols);


% only keep the ft events with triggers
ft_event = ft_event(ismember({ft_event.type},{cfg.trialdef.eventtype, 'DIN_1'}));
tmp = {ft_event.value};
ft_event = ft_event(ismember(tmp,triggers));


if cfg.eventinfo.usePhotodiodeDIN
    hasDIN = sum(ismember({ft_event.value},cfg.eventinfo.photodiodeDIN_str))>0;
    if hasDIN
        photodiodeDIN_toleranceMS = cfg.eventinfo.photodiodeDIN_toleranceMS;
        photodiodeDIN_toleranceSamp = ceil((photodiodeDIN_toleranceMS / 1000) * ft_hdr.Fs);
    else
        warning('No DIN events found, using standard offset for all found triggers');
    end
end

offsetSampAA = ceil((cfg.eventinfo.offsetMS / 1000) * ft_hdr.Fs);
%if ~isfield(cfg.eventinfo, 'offsetMSTCP') error('missing standard offset incase of missing DIN'); end
%offsetSampTCP = ceil((cfg.eventinfo.offsetMSTCP / 1000) * ft_hdr.Fs);




%% go through events and add metadata to trl matrix

ses = cfg.eventinfo.sessionNum;
sesName = cfg.eventinfo.sessionNames{ses};
sesType = find(ismember(cfg.eventinfo.sessionNames,cfg.eventinfo.sessionNames{ses}));
% sesType = ismember(cfg.eventinfo.sessionNames,cfg.eventinfo.sessionNames{ses});

fprintf('FT event count of NS flags (out of %d): %s',length(ft_event),repmat(' ',1,length(num2str(length(ft_event)))));

for i = 1:length(ft_event)
    fprintf(1,[repmat('\b',1,length(num2str(i))),'%d'],i);
    
    if strcmp(ft_event(i).type,cfg.trialdef.eventtype)
        % found an EEG event that we might want to process
        
        switch ft_event(i).value
            case 'FLK2'
                nKeys = length(ft_event(i).orig.keys);                
                defOffset = ceil((15/1000) * ft_hdr.Fs); %TCP offset if DIN triggers not used, determined by photodiode tests
                
                % set column types because Net Station evt files can vary
                ns_evt_cols = {};
                for ns = 1:nKeys
                    ns_evt_cols = cat(1,ns_evt_cols,ft_event(i).orig.keys(ns).key.keyCode);
                end
                cols.pfrq = find(strcmp(ns_evt_cols,'pfrq'));
                if isempty(cols.pfrq)
                    keyboard
                end
                
                pfrq = ft_event(i).orig.keys(cols.pfrq).key.data.data;
                evVal = sprintf('flckr%i',str2num(pfrq));
                
                % get the order of trl columns for this phase and event type
                trl_order = cfg.eventinfo.trl_order.(evVal);
                
                % find where this event type occurs in the list
                eventNumber = find(ismember(cfg.trialdef.eventvalue,evVal));
                if isempty(eventNumber)
                    eventNumber = -1;
                end
                
                if length(eventNumber) == 1 && eventNumber ~= -1
                    % set the times we need to segment before and after the
                    % trigger
                    prestimSec = abs(cfg.eventinfo.prepost{1}(eventNumber,1));
                    poststimSec = cfg.eventinfo.prepost{1}(eventNumber,2);
                    
                    % prestimulus period should be negative
                    prestimSamp = -round(prestimSec * ft_hdr.Fs);
                    poststimSamp = round(poststimSec * ft_hdr.Fs);
                else
                    fprintf('event number not found for %s!\n',evVal);
                    keyboard
                end
                
                % add it to the trial definition
                this_trl = trl_ini;
                
                % get the time of this event
                this_sample = ft_event(i).sample;

                %set offset
                offsetSamp = offsetSampAA; %anti-aliasing offset
%                offsetSamp = offsetSampAA + offsetSampTCP;
                
                % if we're using the photodiode DIN and we find one
                % within the threshold, replace the current sample time
                % with that of the DIN
                if cfg.eventinfo.usePhotodiodeDIN      
%                     
%                     if i~=length(ft_event) && ... %is it not the last event?
%                             strcmp(ft_event(i+1).value,cfg.eventinfo.photodiodeDIN_str) && ...%is next event a DIN?
%                             abs(ft_event(i+1).sample - this_sample) <= photodiodeDIN_toleranceSamp %is it in the tolerance range?
%                         
%                         this_sample = ft_event(i+1).sample;
%                         offsetSamp = offsetSampAA;
%                         
%                     else
%                         warning('No DIN found for %s event %i; using standard TCP offset of %ims', ft_event(i).value, i, cfg.eventinfo.offsetMSTCP);
%                         offsetSamp = offsetSampAA + offsetSampTCP;
%                     end   
                    
                     %matt's DIN detection... seems to be grabbing dins on
                     %either side? this doesn't make sense to me, DINS
                     %should always follow the TCP event
                     try
                         if strcmp(ft_event(i+1).value,cfg.eventinfo.photodiodeDIN_str) && strcmp(ft_event(i-1).value,cfg.eventinfo.photodiodeDIN_str)
                             
                             % if there is a DIN before and after the stim, pick
                             % the closer one
                             preDiff = (ft_event(i-1).sample - this_sample);
                             postDiff = (ft_event(i+1).sample - this_sample);
                             
                             if preDiff < 0 && abs(preDiff) <= photodiodeDIN_toleranceSamp
                                 preFlag = true;
                             else
                                 preFlag = false;
                             end
                             if postDiff <= photodiodeDIN_toleranceSamp
                                 postFlag = true;
                             else
                                 postFlag = false;
                             end
                             
                             if preFlag && ~postFlag
                                 % only the pre-DIN makes sense
                                 this_sample = ft_event(i-1).sample;
                             elseif ~preFlag && postFlag
                                 % only the post-DIN makes sense
                                 this_sample = ft_event(i+1).sample;
                             elseif preFlag && postFlag
                                 % choose the smaller one
                                 if abs(preDiff) < abs(postDiff)
                                     this_sample = ft_event(i-1).sample;
                                 elseif abs(preDiff) > abs(postDiff)
                                     this_sample = ft_event(i+1).sample;
                                 elseif abs(preDiff) == abs(postDiff)
                                     keyboard
                                 end
                             end
                         elseif strcmp(ft_event(i+1).value,cfg.eventinfo.photodiodeDIN_str) && ~strcmp(ft_event(i-1).value,cfg.eventinfo.photodiodeDIN_str) && (ft_event(i+1).sample - this_sample) <= photodiodeDIN_toleranceSamp
                             this_sample = ft_event(i+1).sample;
                         %elseif strcmp(ft_event(i-1).value,cfg.eventinfo.photodiodeDIN_str) && (ft_event(i-1).sample - this_sample) < 0 && abs(ft_event(i-1).sample - this_sample) <= photodiodeDIN_toleranceSamp
                             % apparently the ethernet tags can be delayed
                             % enough that the DIN shows up first
                         %    this_sample = ft_event(i-1).sample;
                             
                         else
                             error('DIN triggers not found, using default offset of %ims\n',defOffset);
                         end
                     catch
                         warning('DIN triggers not used, using default offset of %i samples\n',defOffset);
                         offsetSamp = offsetSamp + defOffset;
                     end
                else
                    offsetSamp = offsetSamp + defOffset; %use average offset from photodiode                      
                end
                
                % prestimulus sample
                this_trl(1) = this_sample + prestimSamp + offsetSamp;
                % poststimulus sample
                this_trl(2) = this_sample + poststimSamp + offsetSamp;
                % offset in samples
                this_trl(3) = prestimSamp;
                
                %event number
                this_trl(4) = eventNumber;
                
                                
                %add trial info
                %{'subn', 'blkn', 'phsn', 'trln', 'ordr', 'stmt', 'pfrq', 'itid', 'corr'};
                for to = 1:length(trl_order)
                    thisInd = find(ismember(trl_order,trl_order{to}));
                    keyInd = find(ismember(ns_evt_cols, trl_order{to}));
                            
                    if ~isempty(thisInd) & ~isempty(keyInd)
                        evtvar = [];
                        switch trl_order{to}
                            case {'subn','blkn','phsn','trln','ordr','pfrq','itid','corr'}
                                evtvar = str2num(ft_event(i).orig.keys(keyInd).key.data.data);
                                if isempty(evtvar)
                                    evtvar = nan;
                                end
                                
                            case {'stmt'}
                                stimtypes = {'Study','Flick','Old','New'};
                                evtvar = find(ismember(stimtypes,ft_event(i).orig.keys(keyInd).key.data.data));
                        end
                        
                        this_trl(timeCols + eventNumCols + thisInd) = evtvar;
                    end
                    
                end
                
                
                % put all the trials together
                trl = cat(1,trl,double(this_trl));
                
                
                
                %       case 'RESP'
                %
                %       case 'FIXT'
                %
                %       case 'PROM'
                %
                %       case 'REST'
                %
                %       case 'REND'
                %
                %       case cfg.eventinfo.photodiodeDIN_str
        end % switch
    end
end
%add event number to trl_oder
fprintf('\n');
