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

%% load new mats in same format as hernan's

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

%% FIGURE 2: 35 preop, 35 postop, all ctrls
roi_list = {fpac_global, WB_regions};
roi_names = {"FPAC","WB"};

figure('Name','NBM Comparisons','Units','normalized','Position',[0.1 0.1 0.8 0.6]);
for s = 1:2 % loop over ipsi and contra NBM
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
        X_preop = nan(height(preop_pats),1);
        for ii = 1:height(preop_pats)
            current_mat = PreOpZ_ACC(:,:,ii);
            avg_zscore = current_mat(nbm_idx, roi_idx);
            X_preop(ii) = mean(avg_zscore, 'omitnan');
        end
        X_postop = nan(height(postop_pats),1);
        for ii = 1:height(postop_pats)
            current_mat = PostOpZ_ACC(:,:,ii);
            avg_zscore = current_mat(nbm_idx, roi_idx);
            X_postop(ii) = mean(avg_zscore, 'omitnan');
        end
        subplot(2,2,(s-1)*2 + r);
        grp1_name = "Preop";
        grp1_FC = X_preop;
        grp2_name = "Postop";
        grp2_FC = X_postop;
        grp3_name = "Controls";
        grp3_FC = X_control;
        paired = 0;
        save = 0; % don't save inside the loop
        node1 = sprintf("%s NBM", nbm_side);  % "ipsi NBM" or "contra NBM"
        node2 = roi_name;
        typeFC = "avg-zscore";
        [new_violinplot] = compareFCacross3grps(...
            grp1_name, grp1_FC, grp2_name, grp2_FC, grp3_name, grp3_FC, ...
            node1, node2, save);
        title(sprintf("%s to %s", node1, roi_name));
    end
end

sgtitle('NBM Functional Connectivity Comparisons', ...
        'FontWeight','bold','FontSize',16);

annotation('textbox',[0 0.91 1 0.05], ...
           'String','Preop (n=38) vs Postop (n=38) vs All Controls (n=106)', ...
           'EdgeColor','none', ...
           'HorizontalAlignment','center', ...
           'FontSize',14, 'FontAngle','italic');

set(gcf,'Units','inches','Position',[1 1 11 6]);  % width=11, height=6 → landscape
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 6]);


%% manual fdr correction
%Ipsi NBM to FPAC:
% Preop-postop: p=0.0158
% Preop-controls: p=0.00713
% 
% Contra NBM to FPAC:
% Preop-postop: p=0.000976
% Postop-controls: p=0.0462
% 
% Ipsi NBM to Thal:
% Preop-postop: p=0.00679
% Postop-controls: p=0.000929
% 
% Ipsi NBM to WB:
% Preop-postop: p=0.00176
% Preop-controls: p=0.00888
% 
% Contra NBM-WB:
% Preop-postop: 0.00478

p = [0.0158, 0.00713, 0.000976, 0.0462, 0.00679, 0.000929, 0.00176, 0.00888, 0.00478];
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

%% SAVE above
% file_title = "Fig2_NBM_FC_matched_preoppostop_allctrls_agecorrected_121625";
% print(gcf, sprintf('%s%s.pdf', save_dir, file_title), '-dpdf', '-bestfit');

%% FIGURE 3: szfree preop/postop/change
bfore_idx = 83;

%% preop plot
for ii = 1:length(sz_free_preop)
    current_mat = PreOpZ_SzFree_ACC(:, :, ii);
    avg_zscore = current_mat(bfore_idx, fpac_global);
    mean_avg_zscore = mean(avg_zscore, 'omitnan');
    disp(mean_avg_zscore);
    X_szfreepreop(ii) = mean_avg_zscore;
end

for ii = 1:length(non_sz_free_preop)
    current_mat = PreOpZ_NonSzFree_ACC(:, :, ii);
    avg_zscore = current_mat(bfore_idx, fpac_global);
    mean_avg_zscore = mean(avg_zscore, 'omitnan');
    disp(mean_avg_zscore);
    X_nonszfree_preop(ii) = mean_avg_zscore;
end

grp1_name = "szfree";
grp1_FC = X_szfreepreop;
grp2_name = "nonszfree";
grp2_FC = X_nonszfree_preop;
paired = 0;
save = 0;
node1 = "ipsiNBM";
node2 = "FPAC";
typeFC = "avg_zscore_ic";
[new_violinplot] = compareFCacross2grps(grp1_name,grp1_FC,grp2_name,grp2_FC,paired,node1,node2,0);

%% postop plot
for ii = 1:length(sz_free_postop)
    current_mat = PostOp_SzFree_Z_ACC(:, :, ii);
    avg_zscore = current_mat(bfore_idx, fpac_global);
    mean_avg_zscore = mean(avg_zscore, 'omitnan');
    disp(mean_avg_zscore);
    X_szfreepostop(ii) = mean_avg_zscore;
end

for ii = 1:length(non_sz_free_postop)
    current_mat = PostOp_NonSzFree_Z_ACC(:, :, ii);
    avg_zscore = current_mat(bfore_idx, fpac_global);
    mean_avg_zscore = mean(avg_zscore, 'omitnan');
    disp(mean_avg_zscore);
    X_nonszfree_postop(ii) = mean_avg_zscore;
end

grp1_name = "szfree";
grp1_FC = X_szfreepostop;
grp2_name = "nonszfree";
grp2_FC = X_nonszfree_postop;
paired = 0;
save = 0;
node1 = "ipsiNBM";
node2 = "FPAC";
typeFC = "avg_zscore_ic";
[new_violinplot] = compareFCacross2grps(grp1_name,grp1_FC,grp2_name,grp2_FC,paired,node1,node2,0);

%% postop-preop change plot ipsi NBM-FPAC
for ii = 1:length(sz_free_postop)
    current_mat_post = PostOp_SzFree_Z_ACC(:, :, ii);
    avg_zscore_post = current_mat_post(bfore_idx, fpac_global);
    mean_avg_zscore_post = mean(avg_zscore_post, 'omitnan');
    disp(mean_avg_zscore_post);
    current_mat_pre = PreOpZ_SzFree_ACC(:, :, ii);
    avg_zscore_pre = current_mat_pre(bfore_idx, fpac_global);
    mean_avg_zscore_pre = mean(avg_zscore_pre, 'omitnan');
    disp(mean_avg_zscore_pre);
    X_szfreepostop_preop(ii) = mean_avg_zscore_post - mean_avg_zscore_pre;
end

for ii = 1:length(non_sz_free_postop)
    current_mat = PostOp_NonSzFree_Z_ACC(:, :, ii);
    avg_zscore_post = current_mat(bfore_idx, fpac_global);
    mean_avg_zscore_post = mean(avg_zscore_post, 'omitnan');
    disp(mean_avg_zscore_post);
    current_mat_pre = PreOpZ_NonSzFree_ACC(:, :, ii);
    avg_zscore_pre = current_mat_pre(bfore_idx, fpac_global);
    mean_avg_zscore_pre = mean(avg_zscore_pre, 'omitnan');
    disp(mean_avg_zscore_pre);
    X_nonszfree_postop_preop(ii) = mean_avg_zscore_post - mean_avg_zscore_pre;
end

grp1_name = "change-postop-preop-szfree";
grp1_FC = X_szfreepostop_preop;
grp2_name = "change-postop-preop-nonszfree";
grp2_FC = X_nonszfree_postop_preop;
paired = 0;
save = 0;
node1 = "ipsiNBM";
node2 = "FPAC";
typeFC = "avg_zscore_ic";
[new_violinplot] = compareFCacross2grps(grp1_name,grp1_FC,grp2_name,grp2_FC,paired,node1,node2,0);

%% postop-preop change plot ipsi NBM-WHOLE BRAIN
for ii = 1:length(sz_free_postop)
    current_mat_post = PostOp_SzFree_Z_ACC(:, :, ii);
    avg_zscore_post = current_mat_post(bfore_idx, 1:111);
    mean_avg_zscore_post = mean(avg_zscore_post, 'omitnan');
    disp(mean_avg_zscore_post);
    current_mat_pre = PreOpZ_SzFree_ACC(:, :, ii);
    avg_zscore_pre = current_mat_pre(bfore_idx, 1:111);
    mean_avg_zscore_pre = mean(avg_zscore_pre, 'omitnan');
    disp(mean_avg_zscore_pre);
    X_szfreepostop_preop(ii) = mean_avg_zscore_post - mean_avg_zscore_pre;
end

for ii = 1:length(non_sz_free_postop)
    current_mat = PostOp_NonSzFree_Z_ACC(:, :, ii);
    avg_zscore_post = current_mat(bfore_idx, 1:111);
    mean_avg_zscore_post = mean(avg_zscore_post, 'omitnan');
    disp(mean_avg_zscore_post);
    current_mat_pre = PreOpZ_NonSzFree_ACC(:, :, ii);
    avg_zscore_pre = current_mat_pre(bfore_idx, 1:111);
    mean_avg_zscore_pre = mean(avg_zscore_pre, 'omitnan');
    disp(mean_avg_zscore_pre);
    X_nonszfree_postop_preop(ii) = mean_avg_zscore_post - mean_avg_zscore_pre;
end

grp1_name = "change-postop-preop-szfree";
grp1_FC = X_szfreepostop_preop;
grp2_name = "change-postop-preop-nonszfree";
grp2_FC = X_nonszfree_postop_preop;
paired = 0;
save = 0;
node1 = "ipsiNBM";
node2 = "WB";
typeFC = "avg_zscore_ic";
[new_violinplot] = compareFCacross2grps(grp1_name,grp1_FC,grp2_name,grp2_FC,paired,node1,node2,0);

set(gcf,'Units','inches','Position',[1 1 11 6]);  % width=11, height=6 → landscape
set(gcf,'PaperOrientation','landscape');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 11 6]);

%% SAVE above
file_title = "FigX_preoppostopchange_ipsiNBM-WB_121625";
print(gcf, sprintf('%s%s.pdf', save_dir, file_title), '-dpdf', '-bestfit');

%% Linear model: FC ~ Outcome + Volume
% needs:
%   FCstruct.preop.preop_pat_list
%   FCstruct.preop.ipsiNBM_to_fpac_global
%   master_table_35.SubjectID_data
%   master_table_35.EngelOutcome1yrNum
%   data_ipsi_nbm_vol_mm3_pre   (volume vector aligned to master_table_35 row order)
%   out_dir (string)

% match master_table rows to order of FC
fc_ids = string(FCstruct.preop.preop_pat_list);
fc_ids_post = string(FCstruct.postop.postop_pat_list); 
master_ids = string(master_table_35.SubjectID_data);
master_ids_01 = upper(strtrim(master_ids + "_01"));
master_ids_03 = upper(strtrim(master_ids + "_03"));

[isMatch, fc_idx] = ismember(master_ids_01, upper(strtrim(fc_ids)));
[isMatch_post, fc_idx_post] = ismember(master_ids_03, upper(strtrim(fc_ids_post)));

master_idx = find(isMatch);
fc_idx = fc_idx(isMatch);

master_idx_post = find(isMatch_post);
fc_idx_post = fc_idx_post(isMatch_post);

% because do not have volume for Epat13_01 which is idx 2 and sz_free
X_szfreepreop_test = X_szfreepreop(:, [1 3:18]); 
X_szfree_postop_test = X_szfreepostop(:, [1 3:18]);
X_szfreepostop_preop_test = X_szfreepostop_preop(:, [1 3:18]);

y_all_pre = [X_szfreepreop_test X_nonszfree_preop]';
y_all_post = [X_szfreepostop X_nonszfree_postop]';
y_all_change = [X_szfreepostop_preop_test X_nonszfree_postop_preop]';

y_all_final = y_all_change;

y_label = "ipsiNBM_to_fpac_global_post-pre";

% covariate: volumes
    % to get volume, use variable (Tboth) from NBM_CODE_FINAL
    upper_ids = vertcat(sz_free_ids, non_sz_free_ids);
    [tf, loc] = ismember(upper_ids, upper(Tboth.Patient));
    upper_ids = upper_ids(tf);   % keep only matched IDs
    loc = loc(tf);               % keep only valid indices
    ipsi_nbm_volume = Tboth.Ipsi_NBM_vol_mm3_preop(loc);
    % ids = find(ismember(upper(Tboth.Patient), upper_ids));
    % ipsi_nbm_volume = Tboth.Ipsi_NBM_vol_mm3_preop(ids);


% Predictor: outcome (szfree vs nonszfree)
outcome = categorical([ ...
    repmat("Sz-Free", 18, 1); ...
    repmat("Non-Sz-Free", 16, 1) ...
]);

FC_vec = y_all_final;
Vol_vec = ipsi_nbm_volume;
Out_vec = categorical(outcome);

FC_vec  = FC_vec(:);
Vol_vec = Vol_vec(:);
Out_vec = Out_vec(:);

% build table
tbl = table(FC_vec, Out_vec, Vol_vec, 'VariableNames', {'FC','Outcome','Volume'});
tbl.Outcome = reordercats(tbl.Outcome, {'Sz-Free','Non-Sz-Free'});
mdl = fitlm(tbl, 'FC ~ Outcome + Volume');

disp(mdl)
fprintf("\nCoefficients (Non-Sz-Free vs Sz-Free), adjusted for Volume:\n");
disp(mdl.Coefficients)

% Plot: FC by outcome (raw points) + model-adjusted group means (holding volume at its mean)
fig = figure('Color','w'); hold on
g = tbl.Outcome;
x = double(g);  % Sz-Free=1, Non-Sz-Free=2
jitter = (rand(size(x)) - 0.5) * 0.15;
scatter(x(g=="Sz-Free")     + jitter(g=="Sz-Free"),     tbl.FC(g=="Sz-Free"),     60, 'filled');
scatter(x(g=="Non-Sz-Free") + jitter(g=="Non-Sz-Free"), tbl.FC(g=="Non-Sz-Free"), 60, 'filled');
vol0 = mean(tbl.Volume, 'omitnan');
pred_tbl = table( ...
    categorical(["Sz-Free"; "Non-Sz-Free"], categories(tbl.Outcome)), ...
    [vol0; vol0], ...
    'VariableNames', {'Outcome','Volume'});
pred_fc = predict(mdl, pred_tbl);
plot([1 2], pred_fc, '-o', 'LineWidth', 2, 'MarkerSize', 8);
xlim([0.5 2.5]);
xticks([1 2]);
xticklabels({'Sz-Free','Non-Sz-Free'});
ylabel(strrep(y_label,'_',' '));
title(sprintf('%s ~ Outcome + Volume', strrep(y_label,'_',' ')));
legend({'Sz-Free (raw)','Non-Sz-Free (raw)','Model-adjusted mean (Vol=mean)'}, ...
       'Location','best');
set(gca, 'FontSize', 14, 'LineWidth', 2, 'FontWeight', 'bold');
box off

% safe_name = regexprep(string(y_label), '[^\w\-]', '_');
% out_pdf = fullfile(out_dir, safe_name + "_LM_outcome_vol.pdf");
% exportgraphics(fig, out_pdf, "ContentType", "vector");
% fprintf("Saved: %s\n", out_pdf);

%% analyze NBM volume sz-free vs non-sz-free

% split ipsi_nbm_volume var
sz_ids  = upper(strtrim(string(sz_free_ids)));
non_ids = upper(strtrim(string(non_sz_free_ids)));

is_sz  = ismember(upper_ids, sz_ids);
is_non = ismember(upper_ids, non_ids);

ipsi_nbm_szfree     = ipsi_nbm_volume(is_sz);
ipsi_nbm_nonszfree  = ipsi_nbm_volume(is_non);

[h, p] = ttest2(ipsi_nbm_szfree, ipsi_nbm_nonszfree);
% p = 0.4964

%% Figure 4 set-up: compare preop vs postop structural (all patients, independent samples)
patients_structdata = readtable('Y:\NBMProject_organized\recreation_hernanspaper\data\NBM_IpsiContraMetrics_PrePost_combined.csv');
metrics = {'Ipsi_NBM_vol_mm3','Contra_NBM_vol_mm3'};
controls_structdata = readtable('Y:\NBMProject_organized\recreation_hernanspaper\data\NBM_ControlStructuralNBMMetrics.csv');
control_metrics = {'Left_NBM_vol_mm3', 'Right_NBM_vol_mm3'};

patients_pre  = patients_structdata(patients_structdata.HasPreop == 1 & patients_structdata.SideofSeizures ~= 3, :);
patients_post = patients_structdata(patients_structdata.HasPostop == 1 & patients_structdata.SideofSeizures ~= 3, :);

group_labels = ["preop","postop"];
group_colors = [
    0,0.4470,0.7410;   % blue (pre)
    0.4940,0.1840,0.5560  % purple (post)
];

fprintf('\n Pre-op vs Post-op (independent samples) \n');
figure('Position',[100 100 320*numel(metrics) 400]);
tiledlayout(1,numel(metrics),"Padding","compact","TileSpacing","compact");

for k = 1:numel(metrics)
    nexttile; hold on;
    metric = metrics{k};
    pre_vals  = patients_pre.([metric '_preop']);
    post_vals = patients_post.([metric '_postop']);
    all_data = [pre_vals; post_vals];     % plot points + boxcharts
    all_groups = [repmat("preop",numel(pre_vals),1);
                  repmat("postop",numel(post_vals),1)];
    group_cat = categorical(all_groups, group_labels,'Ordinal',true);
    for g = 1:numel(group_labels)
        idx = group_cat == group_labels(g);
        if any(idx)
            xvals = g + 0.15*randn(sum(idx),1);
            scatter(xvals, all_data(idx), 40, group_colors(g,:), ...
                'filled','MarkerFaceAlpha',0.6);
            boxchart(g*ones(sum(idx),1), all_data(idx), ...
                'BoxFaceColor', group_colors(g,:), 'BoxWidth',0.6, 'MarkerStyle','none');
        end
    end
    pre_clean = pre_vals(~isnan(pre_vals));
    post_clean = post_vals(~isnan(post_vals));
    [~,pval] = ttest2(pre_clean, post_clean);
    cleanMetric = strrep(metric,'_',' ');
    set(gca,'XTick',1:2,'XTickLabel',cellstr(group_labels));
    ylabel(cleanMetric);
    title(cleanMetric,'FontWeight','bold');
    set(gca,'FontWeight','bold','LineWidth',1.5,'TickDir','out'); box off;
    fprintf('%-25s | p = %.4f (independent)\n', metric, pval);
end
sgtitle('NBM Structural Metrics: Pre-op vs Post-op (Independent)','FontWeight','bold');


%% FIGURE 4: compare preop vs postop (paired) + controls (independent ttests)
% patients with both sessions
Tboth = patients_structdata(patients_structdata.HasPreop == 1 & ...
                            patients_structdata.HasPostop == 1 & ...
                            patients_structdata.SideofSeizures ~= 3, :);

controls_valid = controls_structdata(controls_structdata.HasData == 1, :);

% control laterality matching (29% Left, 71% Right)
prop_left = 0.29;

controls_valid = sortrows(controls_valid, 'ControlID');  % or another stable ID

n_ctrl  = height(controls_valid);
n_left  = round(prop_left * n_ctrl);

rng(42);
perm = randperm(n_ctrl);
ctrl_left_idx  = false(n_ctrl,1);
ctrl_left_idx(perm(1:n_left)) = true;
ctrl_right_idx = ~ctrl_left_idx;

fprintf('Controls: %d total | %d Left (%.1f%%) | %d Right (%.1f%%)\n', ...
    n_ctrl, sum(ctrl_left_idx), 100*sum(ctrl_left_idx)/n_ctrl, ...
    sum(ctrl_right_idx), 100*sum(ctrl_right_idx)/n_ctrl);

fprintf('\n Paired Patients + Controls \n');

metrics = {'Ipsi_NBM_vol_mm3','Contra_NBM_vol_mm3'};

group_colors = [
    0,0.4470,0.7410;      % preop
    0.4940,0.1840,0.5560; % postop
    0.165,0.494,0.263     % controls
];

figure('Position',[100 100 320*numel(metrics) 420]);
tiledlayout(1,numel(metrics),"Padding","compact","TileSpacing","compact");

for k = 1:numel(metrics)
    nexttile; hold on;
    metric = metrics{k};

    % extract data
    pre_vals  = Tboth.([metric '_preop']);
    post_vals = Tboth.([metric '_postop']);

    validMask = ~isnan(pre_vals) & ~isnan(post_vals);
    pre_vals  = pre_vals(validMask);
    post_vals = post_vals(validMask);

    ctrl_vals = NaN(n_ctrl,1);

    if contains(metric,'Ipsi')
        % ipsilateral: left controls use LEFT, right controls use RIGHT
        ctrl_vals(ctrl_left_idx)  = controls_valid.Left_NBM_vol_mm3(ctrl_left_idx);
        ctrl_vals(ctrl_right_idx) = controls_valid.Right_NBM_vol_mm3(ctrl_right_idx);
    elseif contains(metric,'Contra')
        % contralateral: left controls use RIGHT, right controls use LEFT
        ctrl_vals(ctrl_left_idx)  = controls_valid.Right_NBM_vol_mm3(ctrl_left_idx);
        ctrl_vals(ctrl_right_idx) = controls_valid.Left_NBM_vol_mm3(ctrl_right_idx);
    else
        ctrl_vals(:) = NaN;
    end

    ctrl_vals = ctrl_vals(~isnan(ctrl_vals));

    [~,p_pre_post]  = ttest(pre_vals, post_vals);              % paired
    [~,p_pre_ctrl]  = ttest2(pre_vals,  ctrl_vals);            % independent
    [~,p_post_ctrl] = ttest2(post_vals, ctrl_vals);            % independent

    % plot paired patient lines
    % for i = 1:numel(pre_vals)
    %     plot([1 2],[pre_vals(i) post_vals(i)],'-','Color',[0.75 0.75 0.75]);
    % end

    % scatter + boxcharts
    jitter = 0.15;
    
    scatter(1 + jitter*randn(numel(pre_vals),1), ...
            pre_vals, 40, group_colors(1,:), ...
            'filled','MarkerFaceAlpha',0.6,'MarkerEdgeColor','none');
    
    scatter(2 + jitter*randn(numel(post_vals),1), ...
            post_vals, 40, group_colors(2,:), ...
            'filled','MarkerFaceAlpha',0.6,'MarkerEdgeColor','none');
    
    scatter(3 + jitter*randn(numel(ctrl_vals),1), ...
            ctrl_vals, 40, group_colors(3,:), ...
            'filled','MarkerFaceAlpha',0.6,'MarkerEdgeColor','none');

    boxchart(ones(numel(pre_vals),1),  pre_vals,  ...
        'BoxFaceColor',group_colors(1,:), 'BoxWidth',0.4,'MarkerStyle','none');
    
    boxchart(2*ones(numel(post_vals),1), post_vals, ...
        'BoxFaceColor',group_colors(2,:), 'BoxWidth',0.4,'MarkerStyle','none');
    
    boxchart(3*ones(numel(ctrl_vals),1), ctrl_vals, ...
        'BoxFaceColor',group_colors(3,:), 'BoxWidth',0.4,'MarkerStyle','none');

    % significance bars
    allY = [pre_vals; post_vals; ctrl_vals];
    yMax = max(allY);
    yRange = range(allY);
    if yRange == 0, yRange = 1; end

    y1 = yMax + 0.08*yRange;
    y2 = y1   + 0.08*yRange;
    y3 = y2   + 0.08*yRange;

    plot([1 2],[y1 y1],'k','LineWidth',1.5);
    text(1.5, y1+0.02*yRange, getStars(p_pre_post),'FontWeight','bold','HorizontalAlignment','center');

    plot([1 3],[y2 y2],'k','LineWidth',1.5);
    text(2, y2+0.02*yRange, getStars(p_pre_ctrl),'FontWeight','bold','HorizontalAlignment','center');

    plot([2 3],[y3 y3],'k','LineWidth',1.5);
    text(2.5, y3+0.02*yRange, getStars(p_post_ctrl),'FontWeight','bold','HorizontalAlignment','center');

    set(gca,'XTick',1:3,'XTickLabel',{'Preop','Postop','Control'});
    ylabel(strrep(metric,'_',' '));
    title(strrep(metric,'_',' '),'FontWeight','bold');
    set(gca,'FontWeight','bold','LineWidth',1.5,'TickDir','out');
    ylim([0 300]);
    box off;
    fprintf('%-25s | Pre–Post p=%.4f | Pre–Ctrl p=%.4f | Post–Ctrl p=%.4f\n', ...
        metric, p_pre_post, p_pre_ctrl, p_post_ctrl);
end
sgtitle('NBM Structural Metrics: Paired Patients with Control Comparison','FontWeight','bold');

% manual fdr correction
% Ipsi_NBM_vol_mm3
% Pre–Post p=0.0843 
% Pre–Ctrl p=0.0015 
% Post–Ctrl p=0.0000

% Contra_NBM_vol_mm3        
% Pre–Post p=0.6455 
% Pre–Ctrl p=0.2503 
% Post–Ctrl p=0.1626

p = [0.0843, 0.0015, 0.000, 0.6455, 0.2503, 0.1626];
alpha = 0.05;
[p_sorted, idx] = sort(p);
m = numel(p);

p_fdr_sorted = p_sorted .* m ./ (1:m);
p_fdr_sorted = min(1, cummin(flip(p_fdr_sorted)));
p_fdr_sorted = flip(p_fdr_sorted);

p_fdr = zeros(size(p));
p_fdr(idx) = p_fdr_sorted;

sig_fdr = p_fdr < alpha;

T = table(p(:), p_fdr(:), sig_fdr(:), ...
    'VariableNames', {'p_uncorrected','p_fdr','sig_fdr'});


%% FIGURE 4 (violin style): Preop vs Postop (paired) + Controls (independent)
% patients with both sessions
Tboth = patients_structdata(patients_structdata.HasPreop == 1 & ...
                            patients_structdata.HasPostop == 1 & ...
                            patients_structdata.SideofSeizures ~= 3, :);

controls_valid = controls_structdata(controls_structdata.HasData == 1, :);

% control laterality matching (29% Left, 71% Right)
prop_left = 0.29;
controls_valid = sortrows(controls_valid, 'ControlID');

n_ctrl  = height(controls_valid);
n_left  = round(prop_left * n_ctrl);

rng(42);
perm = randperm(n_ctrl);
ctrl_left_idx  = false(n_ctrl,1);
ctrl_left_idx(perm(1:n_left)) = true;
ctrl_right_idx = ~ctrl_left_idx;

fprintf('Controls: %d total | %d Left (%.1f%%) | %d Right (%.1f%%)\n', ...
    n_ctrl, sum(ctrl_left_idx), 100*sum(ctrl_left_idx)/n_ctrl, ...
    sum(ctrl_right_idx), 100*sum(ctrl_right_idx)/n_ctrl);

fprintf('\n Paired Patients + Controls (VIOLIN STYLE) \n');

metrics = {'Ipsi_NBM_vol_mm3','Contra_NBM_vol_mm3'};

for k = 1:numel(metrics)
    metric = metrics{k};

    pre_vals  = Tboth.([metric '_preop']);
    post_vals = Tboth.([metric '_postop']);

    if contains(metric, 'Contra_NBM_vol_mm3')
        data_contra_nbm_vol_mm3_pre = pre_vals;
        data_contra_nbm_vol_mm3_post = post_vals;
    elseif contains(metric, 'Ipsi_NBM_vol_mm3')
        data_ipsi_nbm_vol_mm3_pre = pre_vals;
        data_ipsi_nbm_vol_mm3_post = post_vals;
    end

    validMask = ~isnan(pre_vals) & ~isnan(post_vals);
    pre_vals  = pre_vals(validMask);
    post_vals = post_vals(validMask);

    ctrl_vals = NaN(n_ctrl,1);

    if contains(metric,'Ipsi')
        ctrl_vals(ctrl_left_idx)  = controls_valid.Left_NBM_vol_mm3(ctrl_left_idx);
        ctrl_vals(ctrl_right_idx) = controls_valid.Right_NBM_vol_mm3(ctrl_right_idx);

    elseif contains(metric,'Contra')
        ctrl_vals(ctrl_left_idx)  = controls_valid.Right_NBM_vol_mm3(ctrl_left_idx);
        ctrl_vals(ctrl_right_idx) = controls_valid.Left_NBM_vol_mm3(ctrl_right_idx);
    else
        ctrl_vals(:) = NaN;
    end

    ctrl_vals = ctrl_vals(~isnan(ctrl_vals));

    node1 = "NBM";
    node2 = string(metric);

    figure('Position',[100 100 520 520]);
    compareVOLacross3grps("Preop",  pre_vals, ...
                         "Postop", post_vals, ...
                         "Controls", ctrl_vals, ...
                         node1, node2, 0);
    ylim([0 300]);
    title(strrep(metric,'_',' '),'FontWeight','bold');
    [~,p_pre_post]  = ttest(pre_vals, post_vals);
    [~,p_pre_ctrl]  = ttest2(pre_vals,  ctrl_vals);
    [~,p_post_ctrl] = ttest2(post_vals, ctrl_vals);

    fprintf('%-25s | Pre–Post p=%.4f | Pre–Ctrl p=%.4f | Post–Ctrl p=%.4f\n', ...
        metric, p_pre_post, p_pre_ctrl, p_post_ctrl);
end

sgtitle('NBM Structural Metrics: Violin style (Preop/Postop/Controls)','FontWeight','bold');

% ctrl_right_NBM_vol_mm3 = controls_structdata.Right_NBM_vol_mm3;
% ctrl_left_NBM_vol_mm3 = controls_structdata.Left_NBM_vol_mm3;
% save("Y:/NBMProject_organized/recreation_hernanspaper/data/pat_ipsi_NBM_vol_mm3_post.mat", "data_ipsi_nbm_vol_mm3_post");

%% Hernan Data Community Detection Recreation
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
% Mika Code

% put into list of matrices (where you can click and open each one)
WW_Ctrl_H = squeeze(num2cell(CtrlsZ_H, [1,2]));
% new version of multiscale community
MR_Ctrl_H = multiscale_community_new(WW_Ctrl_H, gamm_range); % shows you the gamma it's at each iteration
% get partition distance
[~, mi_Ctrl_H] = partition_distance(MR_Ctrl_H');

% HERNANS DATA: plot results
% stable partitions
figure
imagesc(mi_Ctrl_H==1);
colormap('jet'); colorbar;
ax1=gca;
ax1.XTick=1:5:length(gamm_range);
ax1.XTickLabel=(gamm_range(1:5:end));
ax1.YTick=1:5:length(gamm_range);
ax1.YTickLabel=(gamm_range(1:5:end));
xlabel('gamma'); %0xticklabels(gamm_range);
ylabel('gamma'); %yticklabels(gamm_range);
title(strcat('Hernan Controls', ' normalized mutual information'))

% filename = 'stablepartitions_Ctrl_H';
% savefig(gcf, fullfile(out_folder, [filename '.fig']));
% saveas(gcf, fullfile(out_folder, [filename, '.pdf']));

% plot smoothed stability
figure
plot(gamm_range, sum(mi_Ctrl_H==1));
xlabel('gamma')
ylabel('smoothed sum of partition distance')
title('Controls smoothed stability')

% filename = 'stablepartitions_Ctrl_H';
% savefig(gcf, fullfile(out_folder, [filename '.fig']));
% saveas(gcf, fullfile(out_folder, [filename, '.pdf']));
% filename = 'smoothedstability_Ctrl_H';
% savefig(gcf, fullfile(out_folder, [filename '.fig']));
% saveas(gcf, fullfile(out_folder, [filename, '.pdf']));

% HERNANS DATA: Interrogating the communities
% 71 aka gamma = 2 is optimal and gives us 9 communities
% idxs = unique(MR_Ctrl_H(71,:));

ROIs_header = readcell("Y:\NBMProject_organized\recreation_hernanspaper\notes\atlases\GT133_ROIs_names_IpCo.xlsx");
region_names = string(ROIs_header);

% Get communities at gamma = 2 (index 71) for hernan's controls
gamma_idx = 71; % gamma = 2.0
idxs = unique(MR_Ctrl_H(gamma_idx,:));

com_regions_H = struct();
for cc = 1:length(idxs)
    com_regions_H.(sprintf('Community_%d', cc)) = ...
        region_names(MR_Ctrl_H(gamma_idx,:) == cc);
end

disp(com_regions_H)

%% FIGURE 5A: Community Detection

% "Y:\NBMProject_organized\recreation_hernanspaper\data\communitydetection_outputs\subjSpecific\age_corrected_avg_zscore_ic\com_regions_ACC_Postop38.mat"
% "Y:\NBMProject_organized\recreation_hernanspaper\data\communitydetection_outputs\subjSpecific\age_corrected_avg_zscore_ic\com_regions_ACC_Preop38.mat"

%% FIGURE 5B: Glass Brain Creation
addpath Y:\NBMProject\glassbrain_communities\111925_glassbrainvisual
% from: create_glass_brain_visual_edited.m
% This is to create glass brain with all necessary tables made.
% If you need .node file, go to Y:\NBMProject\code\get_centroids_of_segmentation.m
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

% save
exportgraphics(figureHandle, figSavePath_png, 'Resolution', 300);
% exportgraphics(figureHandle, figSavePath_pdf, 'ContentType', 'vector'); 

