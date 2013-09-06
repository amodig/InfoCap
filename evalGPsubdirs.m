function evalGPsubdirs(inputPath, outputFolderName, varargin)
% EVALGPSUBDIRS evaluates GP through subfolders in inputPath
%
% The MIT License (MIT)
% Copyright (c) 2013 Arttu Modig
% version 1.0

dirs = selectDirs(inputPath);

for i = 1:length(dirs),
  dirname = dirs(i).name;
  dirpath = [inputPath filesep dirname];
  outputFolderPath = [dirpath filesep outputFolderName];
  % use evalGPdir
  fprintf(1,'Evaluating subfolder: %s\n', dirpath);
  evalGPdir(dirpath, outputFolderPath);
end