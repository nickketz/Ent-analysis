function out = getBlkData(logdata)
%
% filter indata by block returning a cell array for each block
%
% in:
%   logdata = output structure from entReadlogs
%
% out:
%   out = output cell array with the same structure as indata, but each
%   cell is a single block


nblks = max(max(logdata.data.Block));
indblk = logdata.data.Block(:,1) == 1;
for iblk = 2:nblks
    indblk = cat(2,indblk, logdata.data.Block(:,1) == iblk);
end

fnames = fieldnames(logdata.data);
out = {};
for iblk = 1:nblks
    tmpdata = logdata.data;
    blkind = tmpdata.Block(:,1)==iblk;
    for ifname = 1:length(fnames)
        tmpdata.(fnames{ifname})= tmpdata.(fnames{ifname})(blkind,:);
    end
    out{iblk} = tmpdata;
end
        
        
    