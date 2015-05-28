function trl = ent_trialfun_mff(cfg)

% operates using Net Station evt files and event structs

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
    triggers = {'STIM', 'RESP', 'FIXT', 'PROM', 'REST', 'REND', cfg.eventinfo.photodiodeDIN_str};
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

if cfg.eventinfo.usePhotodiodeDIN
    photodiodeDIN_toleranceMS = cfg.eventinfo.photodiodeDIN_toleranceMS;
    photodiodeDIN_toleranceSamp = ceil((photodiodeDIN_toleranceMS / 1000) * ft_hdr.Fs);
end

offsetSamp = ceil((cfg.eventinfo.offsetMS / 1000) * ft_hdr.Fs);

% ft_event_ind = true(1,length(ft_event));
% for i = 1:length(ft_event)
%   if isempty(ft_event(i).value)
%     ft_event_ind(i) = false;
%   end
% end
% ft_event = ft_event(ft_event_ind);

ft_event = ft_event(ismember({ft_event.type},cfg.trialdef.eventtype));

% only keep the ft events with triggers
ft_event = ft_event(ismember({ft_event.value},triggers));


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
                
                % if we're using the photodiode DIN and we find one
                % within the threshold, replace the current sample time
                % with that of the DIN
                if cfg.eventinfo.usePhotodiodeDIN
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
                    elseif strcmp(ft_event(i-1).value,cfg.eventinfo.photodiodeDIN_str) && (ft_event(i-1).sample - this_sample) < 0 && abs(ft_event(i-1).sample - this_sample) <= photodiodeDIN_toleranceSamp
                        % apparently the ethernet tags can be delayed
                        % enough that the DIN shows up first
                        this_sample = ft_event(i-1).sample;
                    end
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
