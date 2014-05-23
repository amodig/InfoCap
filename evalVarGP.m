function model = evalVarGP(inputFilePath, outputFolderPath, latentDim, iters)
% EVALVARGP perform variational bayesian VAR-GP-LVM
%
% input: input file path, output folder, latent dims (force to, optional), dimensions (to use, optional) 
% output: VAR-GP-LVM model of the input data
%
% version 1.4
% The MIT License (MIT)
% Copyright (c) 2013 Arttu Modig

if ~isdeployed, % if not, the code is not compiled
    addpath(genpath('./gp/'))
end

% Fix seeds
randn('seed', 1e5);
rand('seed', 1e5);

% Define constants
scale = 1;
%iters = 1000;
indPoints = 100; % inducing points
%latentDim = 9; % the #latentdims with ARD
dynUsed = 1; % use dynamics
fixedBetaIters = 50;
fixInd = 0;
baseKern = {'rbfard2', 'white'};
dynamicKern = {'rbf', 'bias', 'white'};
reconstrIters = 2000;
learnVariance = 0;
initX ='ppca';
doReconstr = 0;

%% Load data

% read input (residual) data
fprintf(1,'# Reading data file:  %s\n', inputFilePath);
Y = dlmread(inputFilePath,'\t');
% delete zero features
Y = Y(:,any(Y)); % (any ignores entries that are NaN)
% read dimensions
[~,d] = size(Y);

% create outputfolder (if not exist)
fprintf(1,'# Output folder: %s\n', outputFolderPath);
try
    if ~exist(outputFolderPath, 'dir') % searches also subdirectories!
      mkdir(outputFolderPath);
      fprintf(1,'# Created successfully!\n');
    end;
catch err,
    rethrow(err);
end

%% Set up model

if scale,
    origBias = mean(Y, 1);
    origScale = 1./sqrt(var(Y, 1));
    Y = Y - repmat(origBias, size(Y, 1), 1);
    Y = Y.*repmat(origScale, size(Y, 1), 1);
end

% Set up model (from demOilVargplvm)
%options = vargplvmOptions('dtcvar');
%options.kern = {'rbfard2', 'bias', 'white'};
%options.numActive = 100; 
%options.optimiser = 'scg';
%dynUsed = 1;

% Set up model (from demCmu35gplvmVargplvm3.m)
options = vargplvmOptions('dtcvar');
options.kern = baseKern; % default: {'rbfard2', 'white'}
options.numActive = indPoints; % default: 100
options.optimiser = 'scg';

% select the preliminary number of latent dimensions (default: 9)
if (d < latentDim),
  latentDim = d;
end

fprintf(1,'# Creating the model...\n');
if fixInd
    options.fixInducing = 1;
    options.fixIndices=1:indPoints;
end
model = vargplvmCreate(latentDim, d, Y, options);
model = vargplvmParamInit(model, model.m, model.X); 
model.beta=1/(0.01*var(model.m(:)));
model.vardist.covars = 0.5*ones(size(model.vardist.covars)) + 0.001*randn(size(model.vardist.covars));
model.reconstrIters = reconstrIters;

% set timestamps
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
    optionsDyn.learnVariance = learnVariance; % default: 0
    optionsDyn.initX = initX; % default: 'ppca'
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

%% Optimise the model.
display = 1;
model.learnBeta = 1;
model.iters = 0;

fprintf(1,'# Optimising the model for %d iterations...\n',iters);
model = vargplvmOptimise(model, display, iters);
model.iters = model.iters + iters;
model.date = date;

%% Save the model.
fprintf(1,'# Saving model after doing %d iterations... \n',iters);
[~, name, ~] = fileparts(inputFilePath);
capName = name;
capName(1) = upper(capName(1));
if scale,
    save([outputFolderPath filesep 'model-vargplvm-' capName '.mat'], 'model', 'origScale', 'origBias');
else
    save([outputFolderPath filesep 'model-vargplvm-' capName '.mat'], 'model');
end
fprintf(1,'# The model has been saved to %s.\n',['vargpModel' capName '.mat']);

end % end function
