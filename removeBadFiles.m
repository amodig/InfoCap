function fileList = removeBadFiles(directory)
% REMOVEBADFILES remove bad files and subdirs
% (thx to Jonas from stackoverflow)

% select all input files (should be sorted with leading zeros)
fileList = dir(directory);
% remove all folders
isBadFile = cat(1,fileList.isdir); % all directories are bad
% loop to identify hidden files
for iFile = find(~isBadFile)'; % loop only non-dirs
    % on UNIX, hidden files start with a dot
    isBadFile(iFile) = strcmp(fileList(iFile).name(1),'.');
    if ~isBadFile(iFile) && ispc
        % check for hidden Windows files - only works on Windows
        [~,stats] = fileattrib(fullfile(directory,fileList(iFile).name));
        if stats.hidden
            isBadFile(iFile) = true;
        end
    end
end
% remove bad files
fileList(isBadFile) = [];
end