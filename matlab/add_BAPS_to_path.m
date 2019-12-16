baps_path = cd;
addpath(genpath(baps_path));

% % addpath(genpath('/home/ai2/murphyk/matlab/BNT'))
% % fails to add directories which only contain directories but no regular files
% % e.g., BNT/inference
% % This bug has been fixed in matlab 6.5
% 
% global BNT_HOME
% % BNT_HOME = 'C:\kpmurphy\matlab\BNT';
% % BNT_HOME = '/home/ai2/murphyk/matlab/BNT';
% BNT_HOME = 'D:\Matlab\BNT';
% 
% 
% files = {'CPDs', 'general', 'misc', 'graph', 'graph/C', 'Graphics', 'stats1', 'stats2', ...
% 	 'netlab2', 'HMM', 'Kalman', 'Entropic', 'Entropic/Brand', ...
% 	 'inference', 'inference/static', 'inference/dynamic', 'inference/online', ...
% 	 'learning', 'potentials', 'potentials/Tables', ...
% 	 'examples/dynamic', 'examples/dynamic/HHMM', 'examples/dynamic/HHMM/Square', ...
% 	 'examples/dynamic/HHMM/Map', ...
% 	 'examples/dynamic/HHMM/Motif', 'examples/dynamic/SLAM', 'examples/limids', ...
% 	 'examples/static', 'examples/static/Misc',  'examples/static/Models', ...
% 	 'examples/static/Belprop', ...
% 	 'examples/static/Zoubin', 'examples/static/HME', 'examples/static/SCG', ...
% 	 'examples/static/dtree', ...
% 	 'examples/static/StructLearn', 'examples/static/fgraph'};
% 
%  
% %eval(sprintf('addpath %s', BNT_HOME));
% eval(sprintf('addpath ''%s'' ', BNT_HOME));
% % to cope with filenames with spaces, we must use quotes
% % Can use fullfile to get either / for unix or \ for windows
% % use 'isunix' to determine the operating system
% 
% for i=1:length(files)
%   f = files{i};
%   eval(sprintf('addpath ''%s''/%s', BNT_HOME, f));
% end
% 
% 
% 
