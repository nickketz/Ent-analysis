function out = rmBlkData(logdata,badblocks)
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

tmpdata = logdata.data;
bblkind = zeros(size(tmpdata.Block, 1),1);
for iblk = badblocks
    tmpind = (tmpdata.Block(:,1) == iblk);
    bblkind = bblkind + tmpind;
end
bblkind = bblkind==0;

for ifname = 1:length(fnames)
    tmpdata.(fnames{ifname}) = tmpdata.(fnames{ifname})(bblkind,:);
end
out = tmpdata;

        
        
    