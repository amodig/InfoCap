function align_pair(file1, file2, outputFolder, varargin)
% ALIGN_PAIR Align a pair of sequences with CTW.
% ALIGN_PAIR(FILE1, FILE2, OUTPUTFOLDER) writes and returns the
% alignment vectors of two input sequences.
%
% Input:
% sequence files,
% output folder,
% features to use (optional)
%
% Output:
% alignment vectors
%
% version 1.12
% The MIT License (MIT)
% Copyright (c) 2013 Arttu Modig

fprintf(1,'MATLAB user: ');
[status,cmdout] = system('whoami','-echo');

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
end

% read tab-delimited data pair
fprintf(1,'Reading data file: %s\n',file1);
seq1 = dlmread(file1,'\t');
fprintf(1,'Reading data file: %s\n',file2);
seq2 = dlmread(file2,'\t');

% create outputfolder (if not exist)
try
    if ~exist(outputFolder, 'dir')
      mkdir(outputFolder);
    end;
catch err,
    rethrow(err);
end

% algorithm parameter
parDtw = st('dp', 'c');
parCtw = st('th', 0, 'debg', 'n');
parCca = st('d', .8, 'lams', .6);

tic1 = tic;

if isempty(features)
    As = seq1';
    Bs = seq2';
else
    As = seq1(:,features)';
    Bs = seq2(:,features)';
end
Xs = {As, Bs};

% init from dtw
aliDtw = dtw(Xs, [], parDtw);
% ctw
aliCtw = ctw(Xs, aliDtw, [], parCtw, parCca, parDtw);

evalTime = toc(tic1);
fprintf(1,'CTW evaluation time: %f seconds.\n',evalTime);

ali1 = aliCtw.P(:,1);
ali2 = aliCtw.P(:,2);

% write alignment indeces
% NOTE: ignores original filenames
try
    dlmwrite(sprintf('%s/%d_ali_%d.txt', outputFolder, 1,2), ali1,'\t');
    dlmwrite(sprintf('%s/%d_ali_%d.txt', outputFolder, 2,1), ali2,'\t');
    fprintf(1,'Alignment index vectors written succesfully to output folder %s\n', outputFolder);
catch err,
    % could be somehow more helpful...
    fprintf(1,'Could not write output.\n');
    rethrow(err);
end

end

