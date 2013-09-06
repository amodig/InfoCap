function outputFolder = align_subdir(inputPath, outputFolder, varargin)
% ALIGN_subdir Align subdirs of sequences with CTW.
% ALIGN_subdir(DIRECTORY, OUTPUTFOLDER, FEATURES) writes the
% alignment vectors of a subdir. Each subdir must consist of
% repetitions of a certain motion/gesture.
%
% Input:
% directory of sequences,
% output folder,
% features to use (optional)
%
% Output:
% outputfolder location
%
% version 1.01
% The MIT License (MIT)
% Copyright (c) 2013 Arttu Modig

fprintf(1,'align_subdir: Checking arguments...\n');
if isempty(varargin),
    features = [];
elseif (length(varargin)==1),
    features = varargin{1};
else
    error('CTW:argChk','Wrong number of input arguments.');
end

% add path to CTW functions
if ~isdeployed, % if not, the code is not compiled
    addpath(genpath('./ctw'));
    fprintf(1,'Deployed CTW toolbox\n');
end

% select only subdirs
d = dir(inputPath);
isub = [d(:).isdir]; % returns logical vector
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..','__MACOSX','alignment','.svn'})) = []; % remove dumb dirs

% go through subdir tree, with align_dir
for i = 1:length(nameFolds),
  fprintf(1,'Aligning dir: %s\n',nameFolds{i});
  align_dir(sprintf('%s/%s',inputPath,nameFolds{i}),sprintf('%s/%s/%s',inputPath,nameFolds{i},outputFolder));
end

end

