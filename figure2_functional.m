% FIGURE 2: 34 preop, 34 postop, all ctrls
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