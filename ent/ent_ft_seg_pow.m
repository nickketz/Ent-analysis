%segment netstation mff files and create time-frequency data for each
%subject

% initialize the analysis structs
exper = struct;
files = struct;
dirs = struct;
ana = struct;

%% Experiment-specific setup

exper.name = 'ENT';

exper.sampleRate = 250;

% type of NS file for FieldTrip to read; raw or sbin must be put in
% dirs.dataroot/ns_raw; egis must be put in dirs.dataroot/ns_egis
% exper.eegFileExt = 'egis';
% exper.eegFileExt = 'raw';
exper.eegFileExt = 'mff';

% types of events to find in the NS file; these must be the same as the
% events in the NS files; or space_trialfun.m must be set up to find the
% corrct events
exper.eventValues = {'flckr0','flckr6','flckr10','flckr20'};

% pre- and post-stimulus times to read, in seconds (pre is negative).
% Construct as a cell with one Nx2 matrix per session where N is
% length(exper.eventValues{ses}) Order must correspond to the event order
% in exper.eventValues.
exper.prepost = {[-1.0 2.5; -1.0 2.5; -1.0 2.5; -1.0 2.5]};

exper.subjects = {
%     'ENT_Pilot_01',...
%     'ENT_Pilot_04',...
%     'Ent_Pilot_05',...
%     'Ent_Pilot_06',...
%     'Ent_Pilot_07',...
%     'Ent_Pilot_08',...
%      'ENT_09',...    
%      'ENT_10',...
%      'ENT_Pilot_11',...
%      'ENT_12',...
%      'ENT_13',...  %only 6 usable trials!
%      'ENT_14',...  
%      'ENT_15',...
%     'ENT_16',... %no behavioral data?
%     'ENT_17',...
%     'ENT_18',...
%     'ENT_19',...
%     'ENT_20',...
%      'ENT_21',...
%      'ENT_22',...
%      'ENT_23',...
%      'ENT_24',...
%      'ENT_25',...
%      'ENT_26',...
%      'ENT_27',...
      'ENT_28',...
    };

% The sessions that each subject ran; the strings in this cell are the
% directories in dirs.dataDir (set below) containing the ns_egis/ns_raw
% directory and, if applicable, other directories (ns_evt, ns_bci, etc).
% They are not necessarily the session directory names where the FieldTrip
% data is saved for each subject because of the option to combine sessions.
% See 'help create_ft_struct' for more information.
exper.sessions = {{'ses1'}};
%exper.sessions = {{'session_1'},{'session_2'}};

%% set up file and directory handling parameters

% directory where the data to read is located
dirs.subDir = '';
% dirs.dataDir = fullfile(exper.name,'EEG','Sessions','face_house_ratings','eppp',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000),dirs.subDir);
dirs.behDir = fullfile(exper.name,'Behavioral','Sessions','pilot',dirs.subDir);
% dirs.dataDir = fullfile(exper.name,'EEG','Sessions','ftpp',sprintf('%d_%d',exper.prepost(1)*1000,exper.prepost(2)*1000),dirs.subDir);
dirs.dataDir = fullfile(exper.name,'EEG','Sessions','pilot',dirs.subDir);
% Possible locations of the data files (dataroot)
%dirs.serverDir = fullfile(filesep,'Volumes','curranlab','Data');
dirs.serverLocalDir = fullfile(filesep,'Volumes','RAID','curranlab','Data');
dirs.dreamDir = fullfile(filesep,'data','projects','curranlab');
dirs.localDir = fullfile(getenv('HOME'),'Documents','Documents','boulder','Entrain_EEG','Analysis','data');

% pick the right dirs.dataroot
if isfield(dirs,'serverDir') && exist(dirs.serverDir,'dir')
  dirs.dataroot = dirs.serverDir;
  %runLocally = 1;
elseif isfield(dirs,'serverLocalDir') && exist(dirs.serverLocalDir,'dir')
  dirs.dataroot = dirs.serverLocalDir;
  %runLocally = 1;
elseif isfield(dirs,'dreamDir') && exist(dirs.dreamDir,'dir')
  dirs.dataroot = dirs.dreamDir;
  %runLocally = 0;
elseif isfield(dirs,'localDir') && exist(dirs.localDir,'dir')
  dirs.dataroot = dirs.localDir;
  %runLocally = 1;
else
  error('Data directory not found.');
end

% Use the FT chan locs file
files.elecfile = 'GSN-HydroCel-129.sfp';
files.locsFormat = 'besa_sfp';
exper.refChan = 'Cz';
ana.elec = ft_read_sens(files.elecfile,'fileformat',files.locsFormat);

%% Convert the data to FieldTrip structs

% raw data
ana.segFxn = 'seg2ft';

%avoid javaclasspath bug !!!path specific!!!  See http://bit.ly/1o7QlRk 
mydir = pwd;
cd('/Users/nketz/Documents/MATLAB/FTgit/external/egi_mff/')
mff_setup
cd(mydir);

ana.offsetMS = 8; %due to anti-aliasing
%ana.offsetMSTCP = 0; %offset due to tcp delay (average based on timing tests)

ana.continuous = 'yes';
% ana.trialFxn = 'space_trialfun';
ana.trialFxn = 'ent_trialfun_mff';
ana.allowTrialOverlap = true;
ana.renumberSamplesContiguous = true;
% files used when adding metadata to segmented trials
ana.useMetadata = true;
ana.metadata.types = {'nsEvt'};
ana.useExpInfo = 1;
% ana.evtToleranceMS = 8; % 2 samples @ 250 Hz
ana.usePhotodiodeDIN = 1;
ana.photodiodeDIN_toleranceMS = 20;
ana.photodiodeDIN_str = 'DIN ';
if ana.useExpInfo
  % possible sessions and phases
  ana.sessionNames = {'ses1'};
  % needed for seg2ft
  ana.phaseNames = {{''}};
  
  % types of event info to store in trialinfo field; must correspond to
  % values listed in exper.eventValues
%  tmp = {'subn', 'blkn', 'phsn', 'trln', 'ordr', 'stmt', 'pfrq', 'itid'};
  tmp = {'subn', 'blkn', 'phsn', 'trln', 'ordr', 'stmt', 'pfrq', 'itid', 'corr'};
  ana.trl_order.flckr0 = tmp;
  ana.trl_order.flckr6 = tmp;
  ana.trl_order.flckr10 = tmp;
  ana.trl_order.flckr20 = tmp;
end


% process the data after segmentation?
ana.ftFxn = 'ft_freqanalysis';
% ftype is a string used in naming the saved files (data_FTYPE_EVENT.mat)
ana.ftype = 'pow';
ana.overwrite.raw = 1;
ana.overwrite.proc = 1;

% wavelet
cfg_proc = [];
cfg_proc.method = 'wavelet';
cfg_proc.width = 4;
%cfg_proc.toi = -0.8:0.04:3.0;
cfg_proc.toi = -1.0:0.04:2.0;
% evenly spaced frequencies, but not as many as foilim makes
freqstep = (exper.sampleRate./(diff(exper.prepost{1}')*exper.sampleRate)) * 2;
%cfg_proc.foi = 3:freqstep:50;
cfg_proc.foi = 3:freqstep:50;
cfg_proc.output = 'pow';



% any preprocessing? (run after processing artifacts)
cfg_pp = [];
% average rereference
cfg_pp.reref = 'yes';
cfg_pp.refchannel = 'all';
cfg_pp.implicitref = exper.refChan;
% do a baseline correction
cfg_pp.demean = 'yes';
cfg_pp.baselinewindow = [-0.2 0];


% preprocess continuous data in these ways
ana.cfg_cont.lpfilter = 'yes';
ana.cfg_cont.lpfreq = 100;
ana.cfg_cont.hpfilter = 'yes';
ana.cfg_cont.hpfreq = 0.1;
ana.cfg_cont.hpfilttype = 'but';
ana.cfg_cont.hpfiltord = 4;
ana.cfg_cont.bsfilter = 'yes';
ana.cfg_cont.bsfreq = [59 61];

ana.artifact.continuousRepair = true;
ana.artifact.continuousReject = false;
ana.artifact.continuousICA = false;

% artifact settings
ana.artifact.reject = 'complete';
ana.artifact.preArtBaseline = 'yes'; % yes=entire trial

% set up for ftManual/ftAuto
%ana.artifact.type = {'none'};
ana.artifact.type = {'ftAuto'};
% ana.artifact.resumeManArtFT = false;
ana.artifact.resumeManArtContinuous = false;
ana.artifact.resumeICACompContinuous = false;
% negative trlpadding: don't check that time (on both sides) for artifacts.
% IMPORTANT: Not used for threshold artifacts. only use if segmenting a lot
% of extra time around trial epochs. Otherwise set to zero.
ana.artifact.trlpadding = 0;
ana.artifact.artpadding = 0.1;
ana.artifact.fltpadding = 0;

% set up for ftAuto after continuous ICA rejection
ana.artifact.checkAllChan = true;
ana.artifact.thresh = true;
ana.artifact.threshmin = -100;
ana.artifact.threshmax = 100;
ana.artifact.threshrange = 200;
ana.artifact.basic_art = true;
ana.artifact.basic_art_z = 40;
ana.artifact.jump_art = true;
ana.artifact.jump_art_z = 40;

% eog_art is only used with ftAuto
ana.artifact.eog_art = true;
ana.artifact.eog_art_z = 4.5;

% single precision to save space
cfg_pp.precision = 'single';

cfg_proc.keeptrials = 'yes';

% set the save directories
[dirs,files] = mm_ft_setSaveDirs(exper,ana,cfg_proc,dirs,files,'pow');

% create the raw and processed structs for each sub, ses, & event value
[exper] = create_ft_struct(ana,cfg_pp,exper,dirs,files);

process_ft_data(ana,cfg_proc,exper,dirs);


%%
% save the analysis details

backup_orig_AD = true;
% whether to sort by subject number
sortBySubj = false;
% whether to overwite existing subjects in the struct
replaceOrig = true;

% concatenate the additional ones onto the existing ones
saveFile = fullfile(dirs.saveDirProc,'analysisDetails.mat');
if ~exist(saveFile,'file')
  fprintf('Saving analysis details: %s...',saveFile);
  save(saveFile,'exper','ana','dirs','files','cfg_proc','cfg_pp');
  fprintf('Done.\n');
else
  additional_AD_file = fullfile(dirs.saveDirProc,sprintf('analysisDetails%s.mat',sprintf(repmat('_%s',1,length(exper.subjects)),exper.subjects{:})));
  fprintf('Temporarily saving new analysis details: %s...',additional_AD_file);
  save(additional_AD_file,'exper','ana','dirs','files','cfg_proc','cfg_pp');
  
  [exper,ana,dirs,files,cfg_proc,cfg_pp] = mm_mergeAnalysisDetails(saveFile,additional_AD_file,backup_orig_AD,sortBySubj,replaceOrig);
end

%%
%save adfile simple
saveFile = fullfile(dirs.saveDirProc,'analysisDetails.mat');
if ~exist(saveFile,'file')
  fprintf('Saving analysis details: %s...',saveFile);
  save(saveFile,'exper','ana','dirs','files','cfg_proc','cfg_pp');
  fprintf('Done.\n');
else
    %making dated backup
    busaveFile = fullfile(dirs.saveDirProc,['analysisDetails_' datestr(now,30) '.mat']);
    fprintf('Saving backup analysis details: %s...',busaveFile);
    status = system(['mv ' saveFile ' ' busaveFile]);
    fprintf('done\n');
    if status==0
        fprintf('Saving analysis details: %s...',saveFile);
        save(saveFile,'exper','ana','dirs','files','cfg_proc','cfg_pp');
        fprintf('done\n');
    else
        fprintf('\n\nerror making backup, aborted saving!!!\n\n');
    end
end


