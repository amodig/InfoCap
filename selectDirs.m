function list = selectDirs(directory)
% SELECTDIRS makes a list of subdirs in directory
%
% The MIT License (MIT)
% Copyright (c) 2013 Arttu Modig
% version 1.0

% select all files and dirs
list = dir(directory);
% start marking no-dirs
isNotDir = ~cat(1,list.isdir);
% loop through subdirs
for iDir = find(~isNotDir)'; % loop only dirs
  % pseudo-folders start with a dot
  isNotDir(iDir) = strcmp(list(iDir).name(1),'.');
end
% remove non-dirs
list(isNotDir) = [];
end