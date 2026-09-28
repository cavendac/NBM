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