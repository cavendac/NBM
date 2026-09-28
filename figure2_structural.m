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
