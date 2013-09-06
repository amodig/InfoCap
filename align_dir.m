function outputFolder = align_dir(inputPath, outputFolder, varargin)
% ALIGN_dir Align a subdir of sequences with CTW.
% ALIGN_dir(DIRECTORY, OUTPUTFOLDER, FEATURES) writes the
% alignment vectors of a subdir. Subdir must consist of
% repetitions of a single sequence. I.e. In a folder of {seq1, seq2,
% seq3,...} seq1, seq2, seq3, and so on, are repetitions.
%
% Input:
% directory of sequences,
% output folder,
% features to use (optional)
%
% Output:
% outputfolder location
%
% version 2.01
% The MIT License (MIT)
% Copyright (c) 2013 Arttu Modig

fprintf(1,'align_dir: Checking arguments...\n');
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

if ~exist('removeBadFiles.m','file')
    addpath('/var/www/Infocap-Service')
end
% remove all subfolders and hidden files
fileList = removeBadFiles(inputPath);
% check files
if (length(fileList) < 1),
    error('align_pairdir:dataChk','Did not find any data files!');
else if (length(fileList) == 1),
    error('align_pairdir:dataChk','Found only 1 file, CTW needs 2!')
end

% read tab-delimited data (should be sorted with leading-zero-values)
data = cell(1,length(fileList));
for i=1:length(fileList),
    filepath = fullfile(inputPath,fileList(i).name);
    data{i} = dlmread(filepath,'\t');
    fprintf(1,'Read file: %s\n', filepath);
end

% create outputfolder (if not exist)
fprintf(1,'Output folder: %s\n', outputFolder);
try
    if ~exist(outputFolder, 'dir')
      mkdir(outputFolder);
      fprintf(1,'Created successfully!\n');
    end;
catch err,
    rethrow(err);
end

% algorithm parameter
parDtw = st('dp', 'c');
parCtw = st('th', 0, 'debg', 'n');
parCca = st('d', .8, 'lams', .6);

fprintf(1,'Start CTW alignment...\n');
tic1 = tic;
for i = 1:length(data),
    for j = 1:length(data),
      if i~=j,
        if isempty(features)
            As = data{i}';
            Bs = data{j}';
        else
            As = data{i}(:,features)';
            Bs = data{j}(:,features)';
        end
        Xs = {As, Bs};
        tic2 = tic;

        % init from dtw
        aliDtw = dtw(Xs, [], parDtw);
        % ctw
        aliCtw = ctw(Xs, aliDtw, [], parCtw, parCca, parDtw);

        evalTime = toc(tic2);
        fprintf(1,'Pair (%d&%d) CTW evaluation time: %f seconds.\n',i,j,evalTime);

        ali1 = aliCtw.P(:,1);
        ali2 = aliCtw.P(:,2);

        % write alignment indeces
        % NOTE: ignores original filenames
        try
            dlmwrite(sprintf('%s/%d_ali_%d.txt', outputFolder, i,j), ali1,'\t');
            dlmwrite(sprintf('%s/%d_ali_%d.txt', outputFolder, j,i), ali2,'\t');
            fprintf(1,'Alignment index vectors written succesfully to output folder %s\n', outputFolder);
        catch err,
            % this could be somehow more helpfull...
            fprintf(1,'Could not write output.\n');
            rethrow(err);
        end
      end
    end
end
totalTime = toc(tic1);
fprintf(1,'Total CTW evaluation time: %f seconds.\n',totalTime);

end

