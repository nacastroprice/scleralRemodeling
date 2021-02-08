%% Script to rename all files in a directory
% Inputs to change have comments
imageSchema = '^2043_2(\d* )(\d*).*'; % schema of the old names
% directory of the files
projectdir = 'R:\funded_projects\Grytz-R01EY026588-Scleral_remodeling_in_myopia\organ_culture\2043\2\images\';
   dinfo = dir( fullfile(projectdir, '*.tif') ); % assumes files are images
   oldnames = {dinfo.name};
   unwanted = cellfun(@isempty, regexp(oldnames, imageSchema) );
   oldnames(unwanted) = [];
   newnames = regexprep(oldnames, '^(\d*)_2', '20432_');% change last two inputs 
                                                        % (1) Schema of section you want to change (if you want to change it completely it equals initial schema)
                                                        % (2) what you want to change it to
   for K = 1 : length(oldnames)
      movefile( fullfile(projectdir, oldnames{K}), fullfile(projectdir, newnames{K}) );
      disp(length(oldnames)-K)
   end