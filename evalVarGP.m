function model = evalVarGP(inputFilePath, outputFolderPath, varargin)
% EVALVARGP perform variational bayesian VAR-GP-LVM
%
% input: input file path, output folder
% output: VAR-GP-LVM model of the input data
%
% version 1.31
% The MIT License (MIT)
% Copyright (c) 2013 Arttu Modig

scale = 1;

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

if ~isdeployed, % if not, the code is not compiled
    addpath(genpath('./gp/'))
end

% read input (residual) data
fprintf(1,'Reading data file:  %s\n', inputFilePath);
Y = dlmread(inputFilePath,'\t');
% delete zero features
Y = Y(:,any(Y)); % (any ignores entries that are NaN)
% read dimensions
[~,ncol] = size(Y);

% create outputfolder (if not exist)
fprintf(1,'Output folder: %s\n', outputFolderPath);
try
    if ~exist(outputFolderPath, 'dir')
      mkdir(outputFolderPath);
      fprintf(1,'Created successfully!\n');
    end;
catch err,
    rethrow(err);
end

if scale,
    origBias = mean(Y, 1);
    origScale = 1./sqrt(var(Y, 1));
    Y = Y - repmat(origBias, size(Y, 1), 1);
    Y = Y.*repmat(origScale, size(Y, 1), 1);
end

% Set up model (from hgplvm)
%options = fgplvmOptions('ftc');
% Set up model (from demOilVargplvm)
options = vargplvmOptions('dtcvar');
options.kern = {'rbfard2', 'bias', 'white'};
options.numActive = 100; 
options.optimiser = 'scg';
dynUsed = 1;
% select the preliminary number of latent dimensions
if (ncol < 10),
  latentDim = ncol;
else
  latentDim = 10;
end

model = vargplvmCreate(latentDim, ncol, Y, options);

model = vargplvmParamInit(model, model.m, model.X); 
model.vardist.covars = 0.5*ones(size(model.vardist.covars)) + 0.001*randn(size(model.vardist.covars));

dynamicKern = {'rbf', 'bias', 'white'};
seqInd = 1:size(Y,1);
t = 1:max(seqInd);
t = t(seqInd)';

%-------- Add dynamics to the model -----
if dynUsed
    optionsDyn.type = 'vargpTime';
    optionsDyn.t=t;
    optionsDyn.inverseWidth=30;
    %   optionsDyn.vardistCovars = vardistCovarsMult;
    %optionsDyn.seq = seq;
    optionsDyn.learnVariance = 0;
    optionsDyn.initX = 'ppca';
    optionsDyn.regularizeMeans = 0;
    
    % Dynamic kernel:
    kern = kernCreate(t, dynamicKern); % Default: {'rbf','white','bias'}
    
    if strcmp(kern.comp{2}.type, 'white')
        kern.comp{2}.variance = 1e-2; % Usual values: 1e-1, 1e-3
    end
    
    if strcmp(kern.comp{2}.type, 'whitefixed')
        if ~exist('whiteVar')
            whiteVar = 1e-6;
        end
        kern.comp{2}.variance = whiteVar;
        fprintf(1,'# fixedwhite variance: %d\n',whiteVar);
    end
    
    if strcmp(kern.comp{1}.type, 'rbfperiodic')
        if exist('periodicPeriod')
            kern.comp{1}.period = periodicPeriod;
        end
        fprintf(1,'# periodic period: %d\n',kern.comp{1}.period);
    end
    
    % The following is related to the expected number of
    % zero-crossings.(larger inv.width numerator, rougher func)
    if ~strcmp(kern.comp{1}.type,'ou')
        kern.comp{1}.inverseWidth = optionsDyn.inverseWidth./(((max(t)-min(t))).^2);
        kern.comp{1}.variance = 1;
    end
    optionsDyn.kern = kern;
    
    if exist('vardistCovarsMult')
        optionsDyn.vardistCovars = vardistCovarsMult;
    end
    
    % Fill in with default values whatever is not already set
    optionsDyn = vargplvmOptionsDyn(optionsDyn);
    model = vargplvmAddDynamics(model, 'vargpTime', optionsDyn, optionsDyn.t, 0, 0,optionsDyn.seq);
    
    fprintf(1,'# Further calibration of the initial parameters...\n');
    model = vargplvmInitDynamics(model,optionsDyn);
end
modelInit = model;

fixedBetaIters = 50;
% do not learn beta for few iterations for intitilization
if fixedBetaIters > 0
    model.learnBeta = 0;
    display = 1;
    fprintf(1,'# Intitiliazing the model (fixed beta) %d iterations...\n',fixedBetaIters);
    model = vargplvmOptimise(model, display, fixedBetaIters);
    model.fixedBetaIters = fixedBetaIters;
    %disp('# Saving model after optimising beta...')
    %modelWriteResult(model, dataSetName, experimentNo);
end

% Optimise the model.
iters = 1000;
display = 1;
model.learnBeta = 1;
model.iters = 0;

fprintf(1,'# Optimising the model for %d iterations...\n',iters);
model = vargplvmOptimise(model, display, iters);
model.iters = model.iters + iters;
model.date = date;

% Save the results.
[~, name, ~] = fileparts(inputFilePath);
capName = name;
capName(1) = upper(capName(1));
if scale,
    save(['vargpmodel' capName '.mat'], 'model', 'origScale', 'origBias');
else
    save(['vargpmodel' capName '.mat'], 'model');
end
fprintf(1,'The model has been saved to %s.',['vargpmodel' capName '.mat']);

% Reduce dimensions (leaves only non-zero components)
redX = model.X(:, var(model.X) > 1e-6);

% Write reduced data (use original filename)
try
  dlmwrite(sprintf('%s/%s.txt', outputFolderPath, name), redX,'\t');
  fprintf(1,'The latent-vargplvm-data has been written to %s.', ...
    sprintf('%s/%s.txt', outputFolderPath, name));
catch err,
  % could be somehow more helpful...
  fprintf(1,'Could not write output.\n');
  rethrow(err);
end

end
