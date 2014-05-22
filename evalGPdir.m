function evalGPdir(inputDirPath, outputFolderPath, varargin)
% EVALGPDIR computes GP models for data files in inputDirPath 
%
% The MIT License (MIT)
% Copyright (c) 2013 Arttu Modig
% version 1.2
    
%% Compute models

files = removeBadFiles(inputDirPath);
models = cell(length(files),1);
for i = 1:length(files),
    filename = files(i).name;
    filepath = [inputDirPath filesep filename];
    % use Bayesian GP-LVM (vargplvm)
    fprintf(1,'Computing vargplvm model for file: %s\n', filepath);
    model = evalVarGP(filepath, outputFolderPath, 15, 1000);
    model.filename = filename;
    models{i} = model;
end

%% Saving all models

fprintf(1,'# Saving GP-LVM models.\n');
try
  save([outputFolderPath filesep 'models-vargplvm.mat'], 'models');
  fprintf(1,'# All vargplvm models have been saved to %s.\n',[outputFolderPath filesep 'models-vargplvm.mat']);
catch err,
  fprintf(1,'# ERROR: Could not save models.\n');
  rethrow(err);
end

%% Reduce models and save the latent components

for i = 1:length(models),
  for j = 1:length(models),
    if i~=j,
      [redmodel1, P1, dims1] = reduce_model(models{i});
      [redmodel2, ~, ~] = reduce_model(models{j}, P1, dims1); % take P1 latent dim. from dims1
      name1 = redmodel1.filename;
      name2 = redmodel2.filename;
      
      %% TODO: kirjoita paridata (i,j) niin, ett? se tallentaa kaikki redusoidut latentit komponentit,
      % eli kaksi latenttia datasetti? per vertailu (i,j)
      
      % Write reduced data (use original filename)
      pathstr1 = sprintf('%s/%s.txt', outputFolderPath, name1);
      pathstr2 = sprintf('%s/%s.txt', outputFolderPath, name2);

      try
        dlmwrite(sprintf('%s/%s.txt', outputFolderPath, name1), redmodel1.X, '\t');
        fprintf(1,'# The latent-vargplvm-data has been written to: %s\n', pathstr1);
      catch err,
        fprintf(1,'# ERROR: Could not write output to: %s\n', pathstr1);
        rethrow(err);
      end

      try
        dlmwrite(sprintf('%s/%s.txt', outputFolderPath, name2), redmodel2.X, '\t');
        fprintf(1,'# The latent-vargplvm-data has been written to: %s\n', pathstr2);
      catch err,
        fprintf(1,'# ERROR: Could not write output to: %s\n', pathstr2);
        rethrow(err);
      end
    end
  end
end

end
