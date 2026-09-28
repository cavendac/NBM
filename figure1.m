%% nbm paper code part two

%% load FPAC regions
all_regions = readmatrix("Y:\NBMProject_organized\recreation_hernanspaper\notes\atlases\subjSpecific_ROIs_111.xlsx", "OutputType","string");
load("Z:\parcellations\Networks\patient_space\new_fmri_fpac_dk_parc.mat") % this loads fs_lobes
target_regions = vertcat(fs_lobes{1,2}, fs_lobes{2, 2});   % should be a cell array of strings
fpac_indices = nan(size(target_regions));
for i = 1:length(target_regions)
    idx = find(strcmp(all_regions, target_regions{i}));
    if ~isempty(idx)
        fpac_global(i) = idx; % store the index
    end
end
disp(fpac_global)

%% region names
ipsi_nbm = 83;
contra_nbm = 84;
fpac_global = [16;17;45;21;23;25;27;31;33;32;36;38;39;40;41;42;44;48;50;51;79;55;57;59;61;65;66;67;70;72;73;74;75;76;78;82];
WB_regions = 1:111;

%% load new mats in same format as HFG

folder_to_load = "Y:\NBMProject_organized\recreation_hernanspaper\data\subjSpecific_hernanformats_031726";
load("Y:\NBMProject_organized\recreation_hernanspaper\data\subjSpecific_structs_03172026\allcon01_subj_mats_age.mat");
load("Y:\NBMProject_organized\recreation_hernanspaper\data\subjSpecific_structs_03172026\alllatest_postop_subj_mats_age.mat");
load("Y:\NBMProject_organized\recreation_hernanspaper\data\subjSpecific_structs_03172026\alllatest_preop_subj_mats_age.mat");
load("Y:\NBMProject_organized\recreation_hernanspaper\data\subjSpecific_structs_03172026\allpat01_subj_mats_age.mat");

AllPreOpPatsZ_ACC = build_matrix(allpat01_subj_mats_age, all_preop_pats_84, 'age_corrected_avg_zscore_ic', "", true, ...
    folder_to_load, ...
    "AllPreOpPatsZ_ACC");

CtrlZ_all_ACC = build_matrix(allcon01_subj_mats_age, connames, 'age_corrected_avg_zscore', "", true, ...
    folder_to_load, ...
    "CtrlZ_all_ACC");

PreOpZ_ACC = build_matrix(alllatest_preop_subj_mats_age, all_latestpreop, 'age_corrected_avg_zscore_ic', "", true, ...
    folder_to_load, ...
    "PreOpZ_ACC");

PostOpZ_ACC = build_matrix(alllatest_postop_subj_mats_age, all_latestpostop, 'age_corrected_avg_zscore_ic', "", true, ...
    folder_to_load, ...
    "PostOpZ_ACC");

% these szfree and nonszfree vars come from a previous script
PreOpZ_SzFree_ACC = build_matrix(alllatest_preop_subj_mats_age, szfree_preop, 'age_corrected_avg_zscore_ic', "", true, ...
    folder_to_load, ...
    "PreOpZ_SzFree_ACC");

PostOp_SzFree_Z_ACC = build_matrix(alllatest_postop_subj_mats_age, szfree_postop, 'age_corrected_avg_zscore_ic', "", true, ...
    folder_to_load, ...
    "PostOp_SzFree_Z_ACC");

PreOpZ_NonSzFree_ACC = build_matrix(alllatest_preop_subj_mats_age, nonszfree_preop, 'age_corrected_avg_zscore_ic', "", true, ...
    folder_to_load, ...
    "PreOpZ_NonSzFree_ACC");

PostOp_NonSzFree_Z_ACC = build_matrix(alllatest_postop_subj_mats_age, nonszfree_postop, 'age_corrected_avg_zscore_ic', "", true, ...
    folder_to_load, ...
    "PostOp_NonSzFree_Z_ACC");

%% check distribution of FC values

mats = {"AllPreOpPatsZ_ACC.mat";"CtrlZ_all_ACC.mat";"PostOp_NonSzFree_Z_ACC";"PostOpZ_ACC.mat";
    "PostOp_SzFree_Z_ACC.mat";"PreOpZ_ACC.mat";"PreOpZ_NonSzFree_ACC.mat";"PreOpZ_SzFree_ACC.mat"};

num_mats = numel(mats);

% to visually inspect
for ii = 1:num_mats
    varname = mats{ii};
    if contains(varname, ".mat")
        varname = strrep(varname, ".mat", "");
    end
    vals = eval(varname);  % size: [111 111 n_subj]
    n_subj = size(vals, 3);

    % figure out subplot grid
    n_cols = ceil(sqrt(n_subj));
    n_rows = ceil(n_subj / n_cols);

    figure('Name', varname, 'NumberTitle', 'off');
    sgtitle(varname, 'Interpreter', 'none');

    for ss = 1:n_subj
        subplot(n_rows, n_cols, ss);
        imagesc(vals(:,:,ss));
        colorbar;
        title(sprintf('Subj %d', ss));
        axis square;
    end
end

num_mats = numel(mats);
f = figure;
for ii = 1:num_mats
    varname = mats{ii}; 
    if contains(varname, "mat")
        varname = strrep(varname, ".mat", "");
    end
    % string with variable name
    vals = eval(varname);            % get the matrix contents
    vals = vals(:);                  % flatten to vector
    subplot(2, 6, ii);
    histogram(vals, 100, 'Normalization', 'pdf');
    title(sprintf('Matrix %s', varname));
    xlabel('FC Values');
    ylabel('Density');
end
sgtitle('Distribution of FC Matrix Values - Hernan Data Recreation - age_corrected_avg_zscore_ic');
folder = "Y:\NBMProject_organized\recreation_hernanspaper";
cd(folder)

% exportgraphics(f, 'hernandatarecreation_distributionofFCmatrixvalues_fixedzscore_agecorrected.pdf', 'ContentType','vector','BackgroundColor','none');
% savefig(f, 'hernandatarecreation_distributionofFCmatrixvalues_fixedzscore_agecorrected.fig');

%% FIGURE 1: all preop vs all controls
addpath Y:\NBMProject_organized\recreation_hernanspaper\code\final

roi_list = {fpac_global, WB_regions};
roi_names = {"FPAC","WB"};
ipsi_nbm = 83;
contra_nbm = 84;

% get rid of participants not supposed to be in list
% all_preop_pats_84(45) = [];
% all_preop_pats_84(58) = [];
% all_preop_pats_84(64) = [];
% all_preop_pats_84(71) = [];

figure('Name','NBM Comparisons','Units','normalized','Position',[0.1 0.1 0.8 0.6]);
for s = 1:2 % loop over ipsi + contra NBM
    if s == 1
        nbm_side = "ipsi";
        nbm_idx = ipsi_nbm;
    else
        nbm_side = "contra";
        nbm_idx = contra_nbm;
    end

    for r = 1:numel(roi_list) % loop over ROIs
        roi_idx = roi_list{r};
        roi_name = roi_names{r};

        X_control = nan(height(all_controls),1);
        for ii = 1:height(all_controls)
            current_mat = CtrlZ_all_ACC(:,:,ii);
            avg_zscore = current_mat(nbm_idx, roi_idx);   % NBM vs ROI
            X_control(ii) = mean(avg_zscore, 'omitnan');
        end

        X_preop = nan(height(all_preop_pats_84),1);
        for ii = 1:height(all_preop_pats_84)
            current_mat = AllPreOpPatsZ_ACC(:,:,ii);
            avg_zscore = current_mat(nbm_idx, roi_idx);
            X_preop(ii) = mean(avg_zscore, 'omitnan');
        end

        subplot(2,2,(s-1)*2 + r); % index where to go: row = nbm_side, col = roi

        grp1_name = "Preop";
        grp1_FC = X_preop;
        grp2_name = "Controls";
        grp2_FC = X_control;
        paired = 0;
        save = 0; % don't save inside the loop

        node1 = sprintf("%s NBM", nbm_side);  % "ipsi NBM" or "contra NBM"
        node2 = roi_name;
        typeFC = "avg-zscore";

        [new_violinplot] = compareFCacross2grps(...
            grp1_name, grp1_FC, grp2_name, grp2_FC, ...
            paired, node1, node2, save);

        title(sprintf("%s to %s", node1, roi_name));
    end
end

sgtitle('NBM Functional Connectivity Comparisons', ...
        'FontWeight','bold','FontSize',16);

annotation('textbox',[0 0.91 1 0.05], ...
           'String','All Preop (n=84) vs All Controls (n=106); age-corrected FC', ...
           'EdgeColor','none', ...
           'HorizontalAlignment','center', ...
           'FontSize',14, 'FontAngle','italic');

set(gcf,'Units','inches','Position',[1 1 11 6]);  % width=11, height=6 → landscape
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 6]);

% SAVE above
% save_dir = 'C:\Users\cavendac\OneDrive - Vanderbilt\';
% file_title = "Fig1_NBM_FC_allpreop_allctrls_agecorrected_121625";
% print(gcf, sprintf('%s%s.pdf', save_dir, file_title), '-dpdf', '-bestfit');

p = [0.000714, 0.000318, 0.009, 0.0297];
alpha = 0.05;
[p_sorted, idx] = sort(p);
m = numel(p);

p_fdr_sorted = p_sorted .* m ./ (1:m);
p_fdr_sorted = min(1, cummin(flip(p_fdr_sorted)));
p_fdr_sorted = flip(p_fdr_sorted);

p_fdr = zeros(size(p));
p_fdr(idx) = p_fdr_sorted;

sig_fdr = p_fdr < alpha;

T = table(p(:), p_fdr(:), sig_fdr(:), 'VariableNames', {'p_uncorrected','p_fdr','sig_fdr'});