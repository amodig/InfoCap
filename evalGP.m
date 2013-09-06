function model = evalGP(inputFilePath, outputPath, varargin)
% EVALGP perform GP-LVM
% EVALGP
%
% input: input file path, output folder, amount of iterations
% output: GP-LVM model of the input data
%
% version 1.0
% The MIT License (MIT)
% Copyright (c) 2013 Arttu Modig

% Parameters:
% do normalize, no diff, do learn

latentDim = 6;
scale = 1;

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

if ~isdeployed, % if not, the code is not compiled
    addpath(genpath('./gp/'))
end

% read input (residual) data
Y = dlmread(inputfilePath,'\t');

if scale,
    origBias = mean(Y, 1);
    origScale = 1./sqrt(var(Y, 1));
    Y = Y - repmat(origBias, size(Y, 1), 1);
    Y = Y.*repmat(origScale, size(Y, 1), 1);
end

% Set up model (from hgplvm)
options = fgplvmOptions('ftc');

d = size(Y, 2);
model = fgplvmCreate(latentDim, d, Y, options);

% Add dynamics (from hgplvm)
optionsDyn = gpOptions('ftc');

seqInd = 1:size(Y,1);
t = 1:max(seqInd);
t = t(seqInd)';

kern = kernCreate(t, {'rbf', 'white'});
kern.comp{2}.variance = 1e-5;
kern.whiteVariance = kern.comp{2}.variance;
kern.comp{1}.inverseWidth = 5e-3;
kern.comp{1}.variance = 1;
optionsDyn.kern = kern;
diff = 0;
learn = 1;
model = fgplvmAddDynamics(model, 'gpTime', optionsDyn, t, diff, learn);

% Optimise the model.
iters = 1000;
display = 1;

model = fgplvmOptimise(model, display, iters);

% Save the results.
[~, name, ~] = fileparts(inputFilePath);
capName = name;
capName(1) = upper(capName(1));
if scale,
    save(['gpmodel' capName '.mat'], 'model', 'origScale', 'origBias');
else
    save(['gpmodel' capName '.mat'], 'model');
end
fprintf(1,'The model has been saved to %s.',['gpmodel' capName '.mat']);

dlmwrite(sprintf('%s/%s.txt', outputPath, name), model.X,'\t');
fprintf(1,'The laten-GP-data has been written to %s.', ...
    sprintf('%s/%s.txt', outputPath, name));

end