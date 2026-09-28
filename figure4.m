%% community detection code
% Initializations for community detection
%   gamma,resolution parameter
%   gamma>1,        detects smaller modules
%   0<=gamma<1,     detects larger modules
%   gamma=1,        classic modularity (default)
gamm_range = 1.3:0.01:2.5;
all_regions = postop_subj_mats.pat05_03.fMRI_struct.region_names;
out_folder = "Y:\NBMProject_organized\recreation_hernanspaper\data\communitydetection_outputs\subjSpecific\age_corrected_avg_zscore_ic";

% load community mat data
[mat_paths, mat_names] = get_files_in_dir("Y:\NBMProject_organized\recreation_hernanspaper\data\communitydetection_outputs\subjSpecific\age_corrected_avg_zscore_ic\", "*.mat");
for ii = 1:height(mat_paths)
    load(mat_paths(ii))
end

% still using code from Multi_Scale_Community_133.m

% put into list of matrices (where you can click and open each one)
WW_Ctrl_ACC = squeeze(num2cell(CtrlsZ_ACC, [1,2]));
% use new version of multiscale community
MR_Ctrl_ACC = multiscale_community_new(WW_Ctrl_ACC, gamm_range); % shows you the gamma it's at each iteration
% get partition distance
[~, mi_Ctrl_ACC] = partition_distance(MR_Ctrl_ACC');

% plot results
% stable partitions
figure
imagesc(mi_Ctrl_ACC==1);
colormap('jet'); colorbar;
ax1=gca;
ax1.XTick=1:5:length(gamm_range);
ax1.XTickLabel=(gamm_range(1:5:end));
ax1.YTick=1:5:length(gamm_range);
ax1.YTickLabel=(gamm_range(1:5:end));
xlabel('gamma');
ylabel('gamma');
title(strcat('Controls', ' normalized mutual information'))

% filename = 'stablepartitions_Ctrl_ACC';
% savefig(gcf, fullfile(out_folder, [filename '.fig']));
% saveas(gcf, fullfile(out_folder, [filename, '.pdf']));

% plot smoothed stability
figure
plot(gamm_range, sum(mi_Ctrl_ACC==1));
xlabel('gamma')
ylabel('smoothed sum of partition distance')
title('Controls smoothed stability')

ROIs_header = readcell("Y:\NBMProject_organized\recreation_hernanspaper\notes\atlases\subjSpecific_ROIs_111.xlsx");
region_names = string(ROIs_header);

gamma_idx = 71;
idxs = unique(MR_Ctrl_ACC(gamma_idx,:));

com_regions_ACC = struct();
for cc = 1:length(idxs)
    com_regions_ACC.(sprintf('Community_%d', cc)) = ...
        region_names(MR_Ctrl_ACC(gamma_idx,:) == cc);
end

disp(com_regions_ACC)

%% FIGURE 5A: Community Detection

% "Y:\NBMProject_organized\recreation_hernanspaper\data\communitydetection_outputs\subjSpecific\age_corrected_avg_zscore_ic\com_regions_ACC_Postop38.mat"
% "Y:\NBMProject_organized\recreation_hernanspaper\data\communitydetection_outputs\subjSpecific\age_corrected_avg_zscore_ic\com_regions_ACC_Preop38.mat"

%% FIGURE 5B: Glass Brain Creation
addpath Y:\NBMProject\glassbrain_communities\111925_glassbrainvisual
% from: create_glass_brain_visual_edited.m
% If you need to create .node file, go to Y:\NBMProject\code\get_centroids_of_segmentation.m
addpath('Z:\shared_toolboxes\BrainNetViewer_20191031');
base_folder = "Y:\NBMProject\glassbrain_communities\111925_glassbrainvisual\";

%% upload .xlsx and save as .txt
nodes = readtable(strcat(base_folder, "control_subcortical_communities_redNBM_fixedmidline_NO_BG.xlsx"));
writetable(nodes, strcat(base_folder, "node_file_control_subcortical_communities_redNBM_fixedmidline_NO_BG_TABDELIMITED.txt"), 'Delimiter', '\t')

%% create clean .node file

inputFilePath = strcat(base_folder, 'node_file_postop_community_5_TABDELIMITED_NO_BG.txt');
outputFilePath = strcat(base_folder, 'node_file_postop_community_5_TABDELIMITED_NO_BG.node');

inputFile = fopen(inputFilePath, 'r');
outputFile = fopen(outputFilePath, 'w');

while ~feof(inputFile)
    line = fgetl(inputFile);
    parts = strsplit(line);
    formattedLine = '';
    for i = 1:length(parts)
        num = str2double(parts{i});
        if ~isnan(num)
            formattedPart = sprintf('%.6f', num);
        else
            formattedPart = parts{i};
        end
        formattedLine = [formattedLine, formattedPart, ' '];
    end
    fprintf(outputFile, '%s\n', strtrim(formattedLine));
end

fclose(inputFile);
fclose(outputFile);

disp(['Formatted .node file saved as ', outputFilePath]);

%% edit options file
cd(base_folder)
opts = load('control_options_FIXED.mat');   % loads S and EC
opts = load('preop_options_final_v2.mat');   % loads S and EC

%% Load new options file to run below
optionsFilePath = strcat(base_folder,'control_options_final_v2.mat');

%% CREATE GLASS BRAIN

% BEFORE RUNNING, YOU HAVE TO DELETE FIRST ROW OF VARS FROM NODE FILE

% file paths on ernie
% nodesNodePath = 'Y:\NBMProject\glassbrain_communities\preop_postop_communityaffiliations_020525\node_file_postop_community_5_TABDELIMITED.node';
% surfFilePath = "Z:\shared_toolboxes\BrainNetViewer_20191031\Data\SurfTemplate\BrainMesh_Ch2withCerebellum.nv";
% optionsFilePath = 'Y:\NBMProject\glassbrain_communities\preop_options.mat';
% figSavePath = 'Y:\NBMProject\glassbrain_communities\preop_postop_communityaffiliations_020525\postop_community_5_111925.png';
% outputFilePath = "Y:\NBMProject\glassbrain_communities\111925_glassbrainvisual\node_file_control_subcortical_communities_redNBM_fixedmidline_NO_BG_TABDELIMITED.node";
outputFilePath = "Y:\NBMProject\glassbrain_communities\111925_glassbrainvisual\node_file_PREOPPOSTOP_subcortical_communities_redNBM_fixedmidline_NO_BG_TABDELIMITED.node";

preop_postop = "control";
community = "subcortical_communities_allbrainstem_fixedmidline_no_BG";
nodesNodePath = outputFilePath;
surfFilePath = "Z:\shared_toolboxes\BrainNetViewer_20191031\Data\SurfTemplate\BrainMesh_Ch2withCerebellum.nv";
% optionsFilePath = strcat(base_folder, 'preop_options.mat');
figSavePath_png = strcat(base_folder, preop_postop, "_", community, "121725.png");
figSavePath_pdf = strcat(base_folder, preop_postop, "_", community, "121725.pdf");

% cd to where options file is located
% cd(fileparts(optionsFilePath));
load(optionsFilePath);

% now cd to brain net viewer path
cd('Z:\shared_toolboxes\BrainNetViewer_20191031');

% create figure
figureHandle = BrainNet_MapCfg(surfFilePath, nodesNodePath, optionsFilePath);

% save as png
exportgraphics(figureHandle, figSavePath_png, 'Resolution', 300);

% save as pdf
% exportgraphics(figureHandle, figSavePath_pdf, 'ContentType', 'vector'); 