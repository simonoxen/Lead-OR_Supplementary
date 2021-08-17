load('Figure_5-source_data_1.mat')

%% Normalize distance to target by closest_dtt

displace = closest_dtt - mean(closest_dtt);
nrms_displace = nan(size(nrms));

for j = 1:size(nrms,1)
    % mfr
    [displaced_index, values] = displace_row(nrms(j,:), round(displace(j)));
    nrms_displace(j,displaced_index) = values;
    % single unit
    idx = find([single_unit.trajectory]==j);
    for su_idx = idx
        single_unit(su_idx).dtt = single_unit(su_idx).dtt - displace(j)*diff(query_DTT(1:2));
    end
    stn_entry_exit{j} = stn_entry_exit{j} - displace(j)*diff(query_DTT(1:2));
end

nrms = nrms_displace;

%% Sort

[distance_to_stn_sorted ,IDX] = sort(distance_to_stn);
nrms = nrms(IDX,:);
stn_entry_exit = stn_entry_exit(IDX);

% re-order single unit according to sorted trajectories
single_unit_cp = single_unit;
[single_unit_cp.trajectory] = deal(0);
for j = 1:length(IDX)
    [single_unit_cp([single_unit.trajectory] == IDX(j)).trajectory] = deal(j);
end
single_unit = single_unit_cp;
clear single_unit_cp

%% Panel A, B

% nrms
figure('Position', [-1836 17 583 1153])
ax(1) = subplot(8,1,[1:7]);
imagesc(ax(1),query_DTT,1:size(nrms,1),nrms)
caxis(ax(1),[0 3])
colormap(ax(1), 'summer');

% stn entry exit
for j = 1:length(stn_entry_exit)
    arrayfun(@(x) line(ax(1), 'XData', [x x], 'YData',[-0.5 0.5]+j, 'LineWidth', 2), stn_entry_exit{j});
end

% top bottom
N = round(size(nrms,1)*0.2); % 20 percent
N2 = round(size(nrms,1)*0.8); % 20 percent
line('XData', xlim, 'YData', [N N]);
line('XData', xlim, 'YData', [N2 N2]);

% median
top = nrms(1:N,:);
bot = nrms(end-N+1:end,:);
ax(2) = subplot(8,1,8); hold all
plot(ax(2), query_DTT, movmean(nanmedian(top),5), 'b');
plot(ax(2), query_DTT, movmean(nanmedian(bot),5), 'r');

% quantile
top_quantile = movmean(quantile(top, [0.25 0.75]), 5, 2);
bot_quantile = movmean(quantile(bot, [0.25 0.75]), 5, 2);
patch('XData',[query_DTT fliplr(query_DTT)], 'YData', [top_quantile(1,:) fliplr(top_quantile(2,:))], 'FaceColor', 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
patch('XData',[query_DTT fliplr(query_DTT)], 'YData', [bot_quantile(1,:) fliplr(bot_quantile(2,:))], 'FaceColor', 'r', 'FaceAlpha', 0.1, 'EdgeColor', 'none');

% significance
p = [];
for i = 1:size(top,2)
    try
        p(end+1) = signrank(top(:,i),bot(:,i));
    catch % No data remaining after removal of NaNs.
        p(end+1) = nan;
    end
end
[crit_p, adj_p, h] = fdr_bh(p, 0.01);
plot(query_DTT(~isnan(p)),h)

% stn entry exit
in = [];
out = [];
for j = 1:N
    io = stn_entry_exit{j};
    in = [in io(io>2)];
    out = [out io(io<2)];  
end
in = median(in);
out = median(out);
line('XData', [out out], 'YData',ylim, 'LineWidth', 2)
line('XData', [in in], 'YData',ylim, 'LineWidth', 2)

% link
linkaxes(ax,'x');
xlim([-4 10])

%% Panel C, D

% discard
single_unit([single_unit.proportion3ms] > 0.1) = [];
single_unit([single_unit.spike_count] < 100) = [];

% stn entry exit
figure('Position', [-1836 17 583 1153])
ax(1) = subplot(8,1,[1:7]); hold on;
for j = 1:length(stn_entry_exit)
    arrayfun(@(x) line(ax(1), 'XData', [x x], 'YData',[-0.5 0.5]+j, 'LineWidth', 2), stn_entry_exit{j});
end

% cluster
scatter(ax(1), [single_unit.dtt],[single_unit.trajectory],'.');
ax(1).YDir = 'reverse';
ylim([1 size(nrms,1)]);

% distribution
ax(2) = subplot(8,1,8); hold all
x1 = [single_unit.dtt];
s=swarmchart(x1, ones(size(x1)), 'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
s.XJitter = 'none';
s.YJitter = 'density';

% linkaxes
linkaxes(ax,'x');
xlim([-4 10])

%% Cluster density

clusters_MNI = reshape([single_unit.mni],3,[])';

lead path; % need lead-dbs toolbox

stn_n = ea_load_nii(fullfile(ea_space,'atlases','DISTAL Minimal (Ewert 2017)','lh','STN.nii.gz'));
stn_n.dt = [16,0];
stn_n.img = zeros(stn_n.dim);

[I1,I2,I3] = ind2sub(stn_n.dim, 1:numel(stn_n.img));
q_XYZ = ea_vox2mm([I1',I2',I3'], stn_n.mat);

Idx_clusters = rangesearch(clusters_MNI, q_XYZ, 1);
Idx_all_recordings = rangesearch(all_recordings_MNI, q_XYZ, 1);

n_recordings = stn_n;
n_recordings.img = zeros(n_recordings.dim);
n_recordings.img(:) = cellfun(@(x) length(x),  Idx_all_recordings);

cluster_density_n = stn_n;
cluster_density_n.img = zeros(cluster_density_n.dim);
cluster_density_n.img(:) = cellfun(@(x,y) length(x)/length(y), Idx_clusters, Idx_all_recordings);
cluster_density_n.img = cluster_density_n.img .* (n_recordings.img >=10);
cluster_density_n.img(n_recordings.img==0) = -1;
cluster_density_n.fname = 'Figure_5-source_data_2.nii';
ea_write_nii(cluster_density_n);
gzip(cluster_density_n.fname);
delete(cluster_density_n.fname);

%% Aux

function [displaced_index, values] = displace_row(row, displace_value)
    index = find(row);
    values = row(index);
    displaced_index = index - displace_value;
    out_of_bounds = displaced_index<1 | displaced_index>length(row);
    values(out_of_bounds) = [];
    displaced_index(out_of_bounds) = [];

end

% fdr_bh() - Executes the Benjamini & Hochberg (1995) and the Benjamini &
%            Yekutieli (2001) procedure for controlling the false discovery 
%            rate (FDR) of a family of hypothesis tests. FDR is the expected
%            proportion of rejected hypotheses that are mistakenly rejected 
%            (i.e., the null hypothesis is actually true for those tests). 
%            FDR is a somewhat less conservative/more powerful method for 
%            correcting for multiple comparisons than procedures like Bonferroni
%            correction that provide strong control of the family-wise
%            error rate (i.e., the probability that one or more null
%            hypotheses are mistakenly rejected).
%
% Usage:
%  >> [h, crit_p, adj_p]=fdr_bh(pvals,q,method,report);
%
% Required Input:
%   pvals - A vector or matrix (two dimensions or more) containing the
%           p-value of each individual test in a family of tests.
%
% Optional Inputs:
%   q       - The desired false discovery rate. {default: 0.05}
%   method  - ['pdep' or 'dep'] If 'pdep,' the original Bejnamini & Hochberg
%             FDR procedure is used, which is guaranteed to be accurate if
%             the individual tests are independent or positively dependent
%             (e.g., Gaussian variables that are positively correlated or
%             independent).  If 'dep,' the FDR procedure
%             described in Benjamini & Yekutieli (2001) that is guaranteed
%             to be accurate for any test dependency structure (e.g.,
%             Gaussian variables with any covariance matrix) is used. 'dep'
%             is always appropriate to use but is less powerful than 'pdep.'
%             {default: 'pdep'}
%   report  - ['yes' or 'no'] If 'yes', a brief summary of FDR results are
%             output to the MATLAB command line {default: 'no'}
%
%
% Outputs:
%   h       - A binary vector or matrix of the same size as the input "pvals."
%             If the ith element of h is 1, then the test that produced the 
%             ith p-value in pvals is significant (i.e., the null hypothesis
%             of the test is rejected).
%   crit_p  - All uncorrected p-values less than or equal to crit_p are 
%             significant (i.e., their null hypotheses are rejected).  If 
%             no p-values are significant, crit_p=0.
%   adj_p   - All adjusted p-values less than or equal to q are significant
%             (i.e., their null hypotheses are rejected). Note, adjusted 
%             p-values can be greater than 1.
%
%
% References:
%   Benjamini, Y. & Hochberg, Y. (1995) Controlling the false discovery
%     rate: A practical and powerful approach to multiple testing. Journal
%     of the Royal Statistical Society, Series B (Methodological). 57(1),
%     289-300.
%
%   Benjamini, Y. & Yekutieli, D. (2001) The control of the false discovery
%     rate in multiple testing under dependency. The Annals of Statistics.
%     29(4), 1165-1188.
%
% Example:
%   [dummy p_null]=ttest(randn(12,15)); %15 tests where the null hypothesis
%                                       %is true
%   [dummy p_effect]=ttest(randn(12,5)+1); %5 tests where the null
%                                          %hypothesis is false
%   [h crit_p adj_p]=fdr_bh([p_null p_effect],.05,'pdep','yes');
%
%
% Author:
% David M. Groppe
% Kutaslab
% Dept. of Cognitive Science
% University of California, San Diego
% March 24, 2010

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
%
% 5/7/2010-Added FDR adjusted p-values

function [crit_p adj_p h]=fdr_bh(pvals,q,method,report)

if nargin<1,
    error('You need to provide a vector or matrix of p-values.');
else
    if ~isempty(find(pvals<0,1)),
        error('Some p-values are less than 0.');
    elseif ~isempty(find(pvals>1,1)),
        error('Some p-values are greater than 1.');
    end
end

pvals(isnan(pvals))=[];

if nargin<2,
    q=.05;
end

if nargin<3,
    method='pdep';
end

if nargin<4,
    report='no';
end

s=size(pvals);
if (length(s)>2) || s(1)>1,
    [p_sorted, sort_ids]=sort(reshape(pvals,1,prod(s)));
else
    %p-values are already a row vector
    [p_sorted, sort_ids]=sort(pvals);
end
[dummy, unsort_ids]=sort(sort_ids); %indexes to return p_sorted to pvals order
m=length(p_sorted); %number of tests
adj_p=zeros(1,m)*NaN;

if strcmpi(method,'pdep'),
    %BH procedure for independence or positive dependence
    thresh=[1:m]*q/m;
    wtd_p=m*p_sorted./[1:m];
    %compute adjusted p-values
    for a=1:m,
       adj_p(a)=min(wtd_p(a:end)); 
    end
elseif strcmpi(method,'dep')
    %BH procedure for any dependency structure
    denom=m*sum(1./[1:m]);
    thresh=[1:m]*q/denom;
    wtd_p=denom*p_sorted./[1:m];
    %Note, it can produce adjusted p-values greater than 1!
    %compute adjusted p-values
    for a=1:m,
        adj_p(a)=min(wtd_p(a:end));
    end
else
    error('Argument ''method'' needs to be ''pdep'' or ''dep''.');
end
adj_p=reshape(adj_p(unsort_ids),s);

rej=p_sorted<=thresh;
max_id=find(rej,1,'last'); %find greatest significant pvalue
if isempty(max_id),
    crit_p=0;
    h=pvals*0;
else
    crit_p=p_sorted(max_id);
    h=pvals<=crit_p;
end

if strcmpi(report,'yes'),
    n_sig=sum(p_sorted<=crit_p);
    if n_sig==1,
        fprintf('Out of %d tests, %d is significant using a false discovery rate of %f.\n',m,n_sig,q);
    else
        fprintf('Out of %d tests, %d are significant using a false discovery rate of %f.\n',m,n_sig,q);
    end
    if strcmpi(method,'pdep'),
        fprintf('FDR procedure used is guaranteed valid for independent or positively dependent tests.\n');
    else
        fprintf('FDR procedure used is guaranteed valid for independent or dependent tests.\n');
    end
end


end
