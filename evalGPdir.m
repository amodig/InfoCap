function evalGPdir(inputDirPath, outputFolderPath, varargin)
% EVALGPDIR computes GP models for data files in inputDirPath 
%
% The MIT License (MIT)
% Copyright (c) 2013 Arttu Modig
% version 1.1
    
files = removeBadFiles(inputDirPath);

for i = 1:length(files),
    filename = files(i).name;
    filepath = [inputDirPath filesep filename];
    % use Var-GP-LVM
    fprintf(1,'Computing GP-LVM model for: %s\n', filepath);
    evalVarGP(filepath, outputFolderPath);
end