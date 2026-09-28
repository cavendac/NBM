%% TABLE 1: Demographics and Disease Characteristics 

%% first, split into szfree and nonszfree
% define Engel columns by postop time
engel_col_map = containers.Map;
engel_col_map('_02') = 'EngelOutcome6moNum_0_Engel1A_1_Engel1_2_Engel2_3_Engel3_4Engel4';
engel_col_map('_03') = 'EngelOutcome1yrNum_0_Engel1A_1_Engel1_2_Engel2_3_Engel3_4Engel4';
engel_col_map('_04') = 'EngelOutcome2yrNum_0_Engel1A_1_Engel1_2_Engel2_3_Engel3_4Engel4';
engel_col_map('_05') = 'EngelOutcome3yrNum_0_Engel1A_1_Engel1_2_Engel2_3_Engel3_4Engel4';

% new column
combined_table_updated.LatestEngelOutcome = nan(height(combined_table_updated),1);

% loop through each row to assign outcome
for i = 1:height(combined_table_updated)
    scan_name = combined_table_updated.LatestScan_ScanFolderName{i};
    if isempty(scan_name) || scan_name == "None" % skip if empty
        continue
    end
    suffix = scan_name(end-2:end);
    if ~isKey(engel_col_map, suffix)
        warning("No Engel column for suffix '%s' in row %d", suffix, i);
        continue
    end
    engel_col = engel_col_map(suffix);
    patient_id = combined_table_updated.NormalizedID{i};
    rn = find(string(nih_clinical_info_updated.Patient) == patient_id, 1, 'first'); % find matching row
    if isempty(rn)
        warning("Patient %s not found in nih_clinical_info_updated", patient_id);
        continue
    end
    combined_table_updated.LatestEngelOutcome(i) = nih_clinical_info_updated.(engel_col)(rn); % extract engel outcome and assign
end

% fix table mistake
idx_3 = combined_table_updated.SubjectID == "pat18";
combined_table_updated.LatestEngelOutcome(idx_3) = 4;

is_szfree = combined_table_updated.LatestEngelOutcome == 0;
is_nonszfree = combined_table_updated.LatestEngelOutcome ~= 0;

fprintf("\nFound %d seizure-free patients.", sum(is_szfree))
fprintf("\nFound %d non-seizure-free patients.\n", sum(is_nonszfree))

szfree_preop = combined_table_updated.PreopScan_ScanFolderName(is_szfree);
szfree_postop = combined_table_updated.LatestScan_ScanFolderName(is_szfree);
nonszfree_preop = combined_table_updated.PreopScan_ScanFolderName(is_nonszfree);
nonszfree_postop = combined_table_updated.LatestScan_ScanFolderName(is_nonszfree);

% rename bc of string mismatches
szfree_postop(1,1) = "Epat13_03";
szfree_preop(1,1) = "Epat13_01";
szfree_preop(3,1) = "Epat38_01";

%% Handedness

% overall
[counts, values] = groupcounts(master_table_updated.Handedness);
fprintf("Handedness counts (all subjects):\n");
for i = 1:numel(values)
    fprintf("  %s: %d\n", string(values(i)), counts(i));
end
fprintf("\n");

% sz-free
[counts_szfree, values_szfree] = groupcounts(combined_table_updated.Handedness(is_szfree));
fprintf("Handedness counts (Sz-Free):\n");
for i = 1:numel(values_szfree)
    fprintf("  %s: %d\n", values_szfree{i}, counts_szfree(i));
end
fprintf("\n");

% non-sz-free
[counts_nonszfree, values_nonszfree] = groupcounts(combined_table_updated.Handedness(is_nonszfree));
fprintf("Handedness counts (Non-Sz-Free):\n");
for i = 1:numel(values_nonszfree)
    fprintf("  %s: %d\n", values_nonszfree{i}, counts_nonszfree(i));
end

%% Age

preopageatscan = master_table_updated.Pre_opAgeAtScan;
preopageatscan_szfree = combined_table_updated.Pre_opAgeAtScan(is_szfree);
preopageatscan_nonszfree = combined_table_updated.Pre_opAgeAtScan(is_nonszfree);

% overall
avg_age_all = mean(preopageatscan, 'omitnan');
std_age_all = std(preopageatscan, 'omitnan');
fprintf("Pre-op Age at Scan (All subjects): %.2f ± %.2f\n", avg_age_all, std_age_all);

% sz-free group
avg_age_szfree = mean(preopageatscan_szfree, 'omitnan');
std_age_szfree = std(preopageatscan_szfree, 'omitnan');
fprintf("Pre-op Age at Scan (Sz-Free): %.2f ± %.2f\n", avg_age_szfree, std_age_szfree);

% non-sz-free group
avg_age_nonszfree = mean(preopageatscan_nonszfree, 'omitnan');
std_age_nonszfree = std(preopageatscan_nonszfree, 'omitnan');
fprintf("Pre-op Age at Scan (Non-Sz-Free): %.2f ± %.2f\n", avg_age_nonszfree, std_age_nonszfree);

[h, p] = adtest(preopageatscan_szfree);
[h1, p1] = adtest(preopageatscan_nonszfree); % h=0, so normally distributed
[h2, p2] = ttest2(preopageatscan_szfree, preopageatscan_nonszfree);

fprintf("Two-sample t-tests (Sz-Free vs Non-Sz-Free, PreOp Age at Scan):\n");
fprintf("  Preop Age: p=%.3f\n", p2);

if p2 < 0.05
    fprintf("\nPreop Age is significantly different between outcome groups: p=%d\n", p2);
else
    fprintf("\nPreop Age is NOT sig. different between outcome groups: p=%d\n", p2);
end

%% Sex

% overall
[counts, values] = groupcounts(master_table_updated.Gender);
fprintf("Gender counts (all subjects):\n");
for i = 1:numel(values)
    fprintf("  %s: %d\n", string(values(i)), counts(i));
end
fprintf("\n");

% sz-free
[counts_szfree, values_szfree] = groupcounts(combined_table_updated.Gender(is_szfree));
fprintf("Gender counts (Sz-Free):\n");
for i = 1:numel(values_szfree)
    fprintf("  %s: %d\n", string(values_szfree(i)), counts_szfree(i));
end
fprintf("\n");

% non-sz-Free
[counts_nonszfree, values_nonszfree] = groupcounts(combined_table_updated.Gender(is_nonszfree));
fprintf("Gender counts (Non-Sz-Free):\n");
for i = 1:numel(values_nonszfree)
    fprintf("  %s: %d\n", string(values_nonszfree(i)), counts_nonszfree(i));
end

%% Epilepsy Duration

epilepsy_dur = combined_full_table.Duration_Yrs;
epilepsy_dur_szfree = combined_full_table.Duration_Yrs(is_szfree);
epilepsy_dur_nonszfree = combined_full_table.Duration_Yrs(is_nonszfree);

% overall
avg_dur_all = mean(epilepsy_dur, 'omitnan');
std_dur_all = std(epilepsy_dur, 'omitnan');
fprintf("Duration (All subjects): %.2f ± %.2f\n", avg_dur_all, std_dur_all);

% sz-free group
avg_dur_szfree = mean(epilepsy_dur_szfree, 'omitnan');
std_dur_szfree = std(epilepsy_dur_szfree, 'omitnan');
fprintf("Duration (Sz-Free): %.2f ± %.2f\n", avg_dur_szfree, std_dur_szfree);

% non-sz-free group
avg_dur_nonszfree = mean(epilepsy_dur_nonszfree, 'omitnan');
std_dur_nonszfree = std(epilepsy_dur_nonszfree, 'omitnan');
fprintf("Duration (Non-Sz-Free): %.2f ± %.2f\n", avg_dur_nonszfree, std_dur_nonszfree);

[h, p] = adtest(epilepsy_dur_szfree);
[h1, p1] = adtest(epilepsy_dur_nonszfree); % h=0, so normally distributed
[h2, p2] = ttest2(epilepsy_dur_szfree, epilepsy_dur_nonszfree);

fprintf("Two-sample t-tests (Sz-Free vs Non-Sz-Free, Duration):\n");
fprintf("  Duration: p=%.3f\n", p2);

if p2 < 0.05
    fprintf("\nEpilepsy Duration is significantly different between outcome groups: p=%d\n", p2);
else
    fprintf("\nEpilepsy Duration. is NOT sig. different between outcome groups: p=%d\n", p2);
end

%% Onset Side

% overall
[counts, values] = groupcounts(master_table_updated.SideofSeizures);
fprintf("Side of Seizures counts (all subjects):\n");
for i = 1:numel(values)
    fprintf("  %s: %d\n", string(values(i)), counts(i));
end
fprintf("\n");

% sz-free
[counts_szfree, values_szfree] = groupcounts(combined_full_table.Side_1_Left_2_Right_3_Bilateral_(is_szfree));
fprintf("Side of Seizures counts (Sz-Free):\n");
for i = 1:numel(values_szfree)
    fprintf("  %s: %d\n", string(values_szfree(i)), counts_szfree(i));
end
fprintf("\n");

% non-sz-free
[counts_nonszfree, values_nonszfree] = groupcounts(combined_full_table.Side_1_Left_2_Right_3_Bilateral_(is_nonszfree));
fprintf("Side of Seizures counts (Non-Sz-Free):\n");
for i = 1:numel(values_nonszfree)
    fprintf("  %s: %d\n", string(values_nonszfree(i)), counts_nonszfree(i));
end

%% Pathology
% [counts, values] = groupcounts(master_table_35.Pathology);
% 1 - 50 - MTS
% 2 - 10 - some gliosis
% 4 - 2 - normal
    % one of them says "generally normal with focal gliosis from electrodes"
% 5 - 6 - none (NA)
% 6 - 6 - other
    % 1) MTS + encephalomalacia
    % 2) L post temp laminal necrosis
    % 2) reactive CNS tissue
    % 3) mild reactive change
    % 4) hippicampus with neuronal loss and R lateral temporal cortex with
    % gliosis and spongiosis
    % 5) Cortical dysplastic and/or malformative process
% 7 - 9

all_codes = 1:7;
all_labels = [ ...
    "MTS", ...
    "Some gliosis", ...
    "Tumor", ...
    "Normal", ...
    "None", ...
    "Other", ...
    "Other2" ...   
];

% function to compute counts for any array of codes
compute_counts = @(x) cellfun(@(c) sum(x == c), num2cell(all_codes));

% overall
x_all = double(master_table_updated.Pathology);  % make sure numeric
counts_all = compute_counts(x_all);
fprintf("Pathology counts (all subjects):\n");
for i = 1:numel(all_codes)
    fprintf("  %d - %s: %d\n", all_codes(i), all_labels(i), counts_all(i));
end
fprintf("\n");

% sz-free group
x_szfree = double(combined_full_table.Pathology_1_MTS_2_SomeGliosis_3_Tumor_4_Normal_5None_6_Other_de(is_szfree));
counts_szfree = compute_counts(x_szfree);
fprintf("Pathology counts (Sz-Free):\n");
for i = 1:numel(all_codes)
    fprintf("  %d - %s: %d\n", all_codes(i), all_labels(i), counts_szfree(i));
end
fprintf("\n");

% non-sz-free group
x_nonszfree = double(combined_full_table.Pathology_1_MTS_2_SomeGliosis_3_Tumor_4_Normal_5None_6_Other_de(is_nonszfree));
counts_nonszfree = compute_counts(x_nonszfree);
fprintf("Pathology counts (Non-Sz-Free):\n");
for i = 1:numel(all_codes)
    fprintf("  %d - %s: %d\n", all_codes(i), all_labels(i), counts_nonszfree(i));
end

%% Seizure Frequency

sf_FAS_all = master_table_updated.FAS_monthly;
sf_FIAS_all = master_table_updated.FIAS_monthly;
sf_FBTC_all = master_table_updated.FBTC_monthly;

sf_FAS = combined_full_table.FAS_Freq_monthly;
sf_FIAS = combined_full_table.FIAS_Freq_monthly;
sf_FBTC = combined_full_table.FBTC_Freq_monthly;

sf_FAS_szfree = sf_FAS(is_szfree);
sf_FIAS_szfree = sf_FIAS(is_szfree);
sf_FBTC_szfree = sf_FBTC(is_szfree);

sf_FAS_nonszfree = sf_FAS(is_nonszfree);
sf_FIAS_nonszfree = sf_FIAS(is_nonszfree);
sf_FBTC_nonszfree = sf_FBTC(is_nonszfree);

% overall
fprintf("Seizure Frequency (All subjects, monthly):\n");
fprintf("  FAS: %.2f ± %.2f\n", mean(sf_FAS_all,'omitnan'), std(sf_FAS_all,'omitnan'));
fprintf("  FIAS: %.2f ± %.2f\n", mean(sf_FIAS_all,'omitnan'), std(sf_FIAS_all,'omitnan'));
fprintf("  FBTC: %.2f ± %.2f\n", mean(sf_FBTC_all,'omitnan'), std(sf_FBTC_all,'omitnan'));
fprintf("\n");

% sz-free group
fprintf("Seizure Frequency (Sz-Free, monthly):\n");
fprintf("  FAS: %.2f ± %.2f\n", mean(sf_FAS_szfree,'omitnan'), std(sf_FAS_szfree,'omitnan'));
fprintf("  FIAS: %.2f ± %.2f\n", mean(sf_FIAS_szfree,'omitnan'), std(sf_FIAS_szfree,'omitnan'));
fprintf("  FBTC: %.2f ± %.2f\n", mean(sf_FBTC_szfree,'omitnan'), std(sf_FBTC_szfree,'omitnan'));
fprintf("\n");

% non-sz-free group
fprintf("Seizure Frequency (Non-Sz-Free, monthly):\n");
fprintf("  FAS: %.2f ± %.2f\n", mean(sf_FAS_nonszfree,'omitnan'), std(sf_FAS_nonszfree,'omitnan'));
fprintf("  FIAS: %.2f ± %.2f\n", mean(sf_FIAS_nonszfree,'omitnan'), std(sf_FIAS_nonszfree,'omitnan'));
fprintf("  FBTC: %.2f ± %.2f\n", mean(sf_FBTC_nonszfree,'omitnan'), std(sf_FBTC_nonszfree,'omitnan'));
fprintf("\n");

% test normality
[h_FAS_szfree,p_FAS_szfree] = adtest(sf_FAS_szfree);
[h_FAS_nonszfree,p_FAS_nonszfree] = adtest(sf_FAS_nonszfree);
[h_FIAS_szfree,p_FIAS_szfree] = adtest(sf_FIAS_szfree);
[h_FIAS_nonszfree,p_FIAS_nonszfree] = adtest(sf_FIAS_nonszfree);
[h_FBTC_szfree,p_FBTC_szfree] = adtest(sf_FBTC_szfree);
[h_FBTC_nonszfree,p_FBTC_nonszfree] = adtest(sf_FBTC_nonszfree);

% compare szfree vs nonszfree
[h_t_FAS,p_t_FAS] = ranksum(sf_FAS_szfree, sf_FAS_nonszfree);
[h_t_FIAS,p_t_FIAS] = ranksum(sf_FIAS_szfree, sf_FIAS_nonszfree);
[h_t_FBTC,p_t_FBTC] = ranksum(sf_FBTC_szfree, sf_FBTC_nonszfree);

fprintf("Two-sample t-tests (Sz-Free vs Non-Sz-Free, monthly):\n");
fprintf("  FAS: p=%.3f\n", p_t_FAS);
fprintf("  FIAS: p=%.3f\n", p_t_FIAS);
fprintf("  FBTC: p=%.3f\n", p_t_FBTC);

if p_t_FAS < 0.05
    fprintf("\nMonthly FAS Freq. is significantly different between outcome groups: p=%d\n", p2);
else
    fprintf("\nMonthly FAS Freq. is NOT sig. different between outcome groups: p=%d\n", p2);
end
if p_t_FIAS < 0.05
    fprintf("\nMonthly FIAS Freq. is significantly different between outcome groups: p=%d\n", p2);
else
    fprintf("\nMonthly FIAS Freq. is NOT sig. different between outcome groups: p=%d\n", p2);
end
if p_t_FBTC < 0.05
    fprintf("\nMonthly FBTC Freq. is significantly different between outcome groups: p=%d\n", p2);
else
    fprintf("\nMonthly FBTC Freq. is NOT sig. different between outcome groups: p=%d\n", p2);
end

%% Previous Surgery 

% overall
[counts, values] = groupcounts(master_table_updated.PreviousSurgery);
fprintf("Previous Surgery counts (all subjects):\n");
for i = 1:numel(values)
    if isnan(values(i))
        valStr = "NaN";         
    else
        valStr = string(values(i));
    end
    fprintf("  %s: %d\n", valStr, counts(i));
end
fprintf("\n");

% sz-free
[counts_szfree, values_szfree] = groupcounts(combined_full_table.PreviousSurgery_0_No_1_Yes_(is_szfree));
fprintf("Previous Surgery counts (Sz-Free):\n");
for i = 1:numel(values_szfree)
    if isnan(values_szfree(i))
        valStr = "NaN";
    else
        valStr = string(values_szfree(i));
    end
    fprintf("  %s: %d\n", valStr, counts_szfree(i));
end
fprintf("\n");

% non-sz-free
[counts_nonszfree, values_nonszfree] = groupcounts(combined_full_table.PreviousSurgery_0_No_1_Yes_(is_nonszfree));
fprintf("Previous Surgery counts (Non-Sz-Free):\n");
for i = 1:numel(values_nonszfree)
    if isnan(values_nonszfree(i))
        valStr = "NaN";
    else
        valStr = string(values_nonszfree(i));
    end
    fprintf("  %s: %d\n", valStr, counts_nonszfree(i));
end
fprintf("\n");

%% Create latest scan date column

combined_full_table.LatestScanDate = NaT(height(combined_full_table),1);  % datetime column
scan_suffix = extractBetween(combined_full_table.LatestScan_ScanFolderName, strlength(combined_full_table.LatestScan_ScanFolderName)-1, strlength(combined_table_updated.LatestScan_ScanFolderName));

is_6mo = scan_suffix == "02";
is_1yr = scan_suffix == "03";
is_2yr = scan_suffix == "04";
is_3yr = scan_suffix == "05";

combined_full_table.LatestScanDate(is_6mo) = combined_full_table.x6MoPost_opScanDate(is_6mo);
combined_full_table.LatestScanDate(is_1yr) = combined_full_table.x1YrPost_opScanDate(is_1yr);
combined_full_table.LatestScanDate(is_2yr) = combined_full_table.x2YrPost_opScanDate(is_2yr);
combined_full_table.LatestScanDate(is_3yr) = combined_full_table.x3YrPost_opScanDate(is_3yr);

combined_full_table.LatestScanDate.Format = 'MM-dd-yy';
combined_full_table.LatestScan_ScanFolderName(5) = "Epat13_03";

head(combined_full_table(:, {'LatestScan_ScanFolderName','LatestScanDate'}))

%% check and fix dates
% Fix LatestScanDate for pat27
idx = find(combined_full_table.SubjectID == "pat27");
combined_full_table.LatestScanDate(idx) = datetime('02-18-2017', 'InputFormat', 'MM-dd-yyyy');
idx_2 = find(combined_full_table.SubjectID == "Epat09");
combined_full_table.LatestScanDate(idx_2) = datetime('03-02-2020', 'InputFormat', 'MM-dd-yyyy');

% first surgery date
d1 = datetime(combined_full_table.First_Surgery_Date, 'InputFormat', 'dd-MMM-yyyy');
yr1 = year(d1);
d1(yr1 < 100) = d1(yr1 < 100) + calyears(2000);  % e.g., 0018 -> 2018
combined_full_table.First_Surgery_Date = d1;

% latest scan date
d2 = datetime(combined_full_table.LatestScanDate, 'InputFormat', 'MM-dd-yyyy');
yr2 = year(d2);
d2(yr2 < 100) = d2(yr2 < 100) + calyears(2000);
combined_full_table.LatestScanDate = d2;

% double check that all first surgery dates come before latest date
is_valid = combined_full_table.First_Surgery_Date < combined_full_table.LatestScanDate;
invalid_idx = find(~is_valid);

if ~isempty(invalid_idx)
    fprintf("These rows have First_Surgery_Date >= LatestScanDate:\n");
    disp(combined_full_table(invalid_idx, {'SubjectID','First_Surgery_Date','LatestScanDate'}));
else
    disp("All First_Surgery_Date values occur before LatestScanDate.");
end

%% add new surgery date col
combined_full_table.SurgeryDate = combined_full_table.First_Surgery_Date;

validRows = ~isnat(d1) & ~isnat(d2);
days_between = days(d2(validRows) - d1(validRows));
combined_full_table.TimeFromSurgerytoLatestPostopScan = NaN(height(combined_full_table), 1);
combined_full_table.TimeFromSurgerytoLatestPostopScan(validRows) = days(d2(validRows) - d1(validRows));

% d1 = first surg date
% d2 = latest scan date
diff_months = calmonths(between(d1, d2, 'months'));

% add to table
combined_full_table.MonthsToLatestScan = diff_months;

%% Time to latest scan 

follow_up = combined_full_table.MonthsToLatestScan;
follow_up_szfree = combined_full_table.MonthsToLatestScan(is_szfree);
follow_up_nonszfree = combined_full_table.MonthsToLatestScan(is_nonszfree);

% overall
avg_dur_all = mean(follow_up, 'omitnan');
std_dur_all = std(follow_up, 'omitnan');
fprintf("Time to latest scan (All subjects): %.2f ± %.2f\n", avg_dur_all, std_dur_all);

% sz-Free group
avg_dur_szfree = mean(follow_up_szfree, 'omitnan');
std_dur_szfree = std(follow_up_szfree, 'omitnan');
fprintf("Duration (Sz-Free): %.2f ± %.2f\n", avg_dur_szfree, std_dur_szfree);

% non-sz-free group
avg_dur_nonszfree = mean(follow_up_nonszfree, 'omitnan');
std_dur_nonszfree = std(follow_up_nonszfree, 'omitnan');
fprintf("Duration (Non-Sz-Free): %.2f ± %.2f\n", avg_dur_nonszfree, std_dur_nonszfree);

[h, p] = adtest(follow_up_szfree);
[h1, p1] = adtest(follow_up_nonszfree); % h=0, so normally distributed
[h2, p2] = ttest2(follow_up_szfree, follow_up_nonszfree);

if p2 < 0.05
    fprintf("Time to latest scan is significantly different between outcome groups: p=%d", p2);
else
    fprintf("\nTime to latest scan is NOT sig. different between outcome groups: p=%d\n", p2);
end

%% History of FBTC
fprintf("\n History of FBTC:");
fprintf("\n n=%d patients out of the whole group have history of FBTC", sum(combined_full_table.FBTC_Freq_monthly ~= 0));
fprintf("\n n=%d patients seizure-free patients have history of FBTC", sum(combined_full_table.FBTC_Freq_monthly(is_szfree) ~= 0));
fprintf("\n n=%d patients NON-seizure-free patients have history of FBTC\n", sum(combined_full_table.FBTC_Freq_monthly(is_nonszfree) ~= 0));

%% Control demographics
control_table = readtable("C:\Users\cavenda\OneDrive - Vanderbilt\Projects\Project_NBM\demographics\NBMproj_control_cohort.xlsx", "Sheet", "Both");

% Age Comparison
[h_age, p_age, ci, stats] = ttest2(control_table.Control_Age, control_table.Patient_Age);

fprintf("Age comparison between Control_Age and Patient_Age:\n");
fprintf("  t-statistic = %.2f, p-value = %.4f\n\n", stats.tstat, p_age);

% gender
Gender = categorical(control_table.Gender_All, [0 1], {'Male','Female'});
[gender_counts, gender_vals] = groupcounts(Gender);
fprintf("Gender distribution:\n");
for i = 1:numel(gender_vals)
    fprintf("  %s: %d\n", string(gender_vals(i)), gender_counts(i));
end
fprintf("\n");

% handedness
Handedness = categorical(control_table.Handedness_All, [0 1], {'LH','RH'});
[hand_counts, hand_vals] = groupcounts(Handedness);
fprintf("Handedness distribution:\n");
for i = 1:numel(hand_vals)
    fprintf("  %s: %d\n", string(hand_vals(i)), hand_counts(i));
end
fprintf("\n");

% group vs gender chi square
Group = categorical(control_table.Group);
[tab_gender, chi2stat_gender, p_chi_gender] = crosstab(Group, Gender);
fprintf("Chi-square test: Group vs Gender:\n");
disp(tab_gender);
fprintf("Chi-square statistic = %.2f, p-value = %.4f\n\n", chi2stat_gender, p_chi_gender);

% group vs handedness chi square
[tab_hand, chi2stat_hand, p_chi_hand] = crosstab(Group, Handedness);
fprintf("Chi-square test: Group vs Handedness:\n");
disp(tab_hand);
fprintf("Chi-square statistic = %.2f, p-value = %.4f\n\n", chi2stat_hand, p_chi_hand);