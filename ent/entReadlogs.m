function outdata = entReadlogs(dirs,cfg)
%Ent analysis script that reads in EntAssoc_*.txt files in 'logs/'
%
%   input: 
%       dirs: struct containing relevent directories for this experiment
%       cfg
%        filefiltstr: regexp filter on file listing, default is 'EntAssoc_*_Sub[0-9]{1,2}\.txt'
%        badsubs: subject numbers to exclude from analysis
%
%   output:
%       outdata: struc with headernames as matricies with Subs on the
%       columns

%find logs
logdir = [dirs.dataroot filesep dirs.behDir];
files = dir([logdir filesep '*.txt']);
files = {files.name};
%filter
if exist('cfg','var')
    if isfield(cfg, 'filefiltstr')
        filesmatch = regexp(files,cfg.filefiltstr);
    else
        filesmatch = regexp(files,'EntAssoc_[0-9]+_Sub[0-9]{1,2}\.txt');
    end       
else
    filesmatch = regexp(files,'EntAssoc_[0-9]+_Sub[0-9]{1,2}\.txt');
    cfg = [];
end
filesmatch = ~cellfun(@isempty,filesmatch);
if sum(filesmatch)==0
    error('no log files found');
end
files = files(filesmatch);

if isfield(cfg,'badsubs')
    %remove bad subs from files
    badsubs = zeros(1,length(files));
    for ilog = 1:length(files)
        mystart = regexp(files{ilog},'Sub');
        myend = regexp(files{ilog},'\.txt');
        subno = str2num(files{ilog}(mystart+3:myend-1));
        if sum(subno==cfg.badsubs)
            badsubs(ilog) = 1;
        end
    end
    files = files(~badsubs);
end

%get var names and types from header
temp = textread([logdir filesep files{1}],'%s');
vars = strread(temp{1},'%s','delimiter',';');
fstring = [];
vartype = {};
for ivar = 1:length(vars)
    switch vars{ivar}(end)
        case '$'
            fstring = [fstring '%s '];
            %outdata.(vars{ivar}(1:end-1)) = {};
            vartype{ivar} = '$';
        case '|'
            fstring = [fstring '%d '];
            %outdata.(vars{ivar}(1:end-1)) = [];
            vartype{ivar} = '|';
        case '#'
            fstring = [fstring '%f '];
            %outdata.(vars{ivar}(1:end-1)) = [];
            vartype{ivar} = '#';
        otherwise 
            fstring = [fstring '%s '];
            %outdata.(vars{ivar}(1:end-1)) = {};
            vartype{ivar} = '$';
    end
    vars{ivar} = vars{ivar}(1:end-1);    
end

%get data from logs
%mydata = cell(1,filesmatch);
cdata = cell(1,length(vars)+1);
for ilog = 1:length(files)    
    %read in log file
    logfilename = [logdir filesep files{ilog}];
    fid = fopen(logfilename,'r');
    mydata = textscan(fid, fstring, 'Headerlines', 1, 'Delimiter', ';', 'TreatAsEmpty', {'na'});        
    fclose(fid);
    
    for ivar = 1:length(vars)
        cdata{ivar} = cat(2,cdata{ivar}, mydata{ivar});
    end
    cdata{ivar+1} = cat(2,cdata{ivar+1},files(ilog));
end

for ivar = 1:length(vars)
    vardata.(vars{ivar}) = cdata{ivar};
end
outdata.data = vardata;
outdata.logname = cdata{ivar+1};
outdata.vartype = vartype;
outdata.cfg = cfg;
