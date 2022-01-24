lead path
load('Figure_5-source_data_1.mat')

%% NRMS and STN distance transformation for xcorr calculation

stnd_threshold = 2;
nrms_threshold = 2;

nrms_proc = nrms;
nrms_proc(nrms<1) = 1;
nrms_proc = (0.5 - atan(nrms_threshold - nrms_proc') / pi)';
stnd_proc = (0.5 - atan(stnd' - stnd_threshold) / pi)';

%% xcorr

r = [];
max_r = [];
lag_at_max = [];
r_at_zero = [];

for i = 1:size(nrms_proc,1)
    
    [x, lags] = xcorr(nrms_proc(i,:), stnd_proc(i,:), 'biased');
    r = [r x];
    
    [m,j] = max(x);
    max_r = [max_r m];
    lag_at_max = [lag_at_max lags(j)];
    r_at_zero = [r_at_zero x(lags==0)];
    
end

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

%% Sort according to max xcorr value

[~ ,IDX] = sort(max_r,'descend');
max_r = max_r(IDX);
distance_to_stn_sorted = distance_to_stn(IDX);
nrms = nrms(IDX,:);
stn_entry_exit = stn_entry_exit(IDX);
brsh = brsh(IDX);
lag_at_max = lag_at_max(IDX);
nrms_proc = nrms_proc(IDX,:);
stnd_proc = stnd_proc(IDX,:);

%% plot sorted trajectories and separate ones with low xcorr value

half = 0.5;
N1 = round(size(nrms,1) * half);

% nrms
f=figure('Position', [-1836 17 650 1153]);
ax(1) = subplot(10,2,1:2:(half*20));
imagesc(ax(1),query_DTT,1:N1,nrms(1:N1,:)); colormap(ax(1), 'summer');caxis(ax(1),[0 3]); %colorbar('northoutside')
% imagesc(ax(1),1, 1:N1, max_r(1:N1)'); colormap(ax(1), ea_redblue); caxis([0 0.5]); colorbar('northoutside')
% imagesc(ax(1),1, 1:N1, 0.1*abs(lag_at_max(1:N1)')); colormap(ax(1), ea_redblue); caxis([0 2]); colorbar('northoutside')

set(gca,'YTickLabel',{})

% stn entry exit
for j = 1:N1
    arrayfun(@(x) line(ax(1), 'XData', [x x], 'YData',[-0.5 0.5]+j, 'LineWidth', 2), stn_entry_exit{j});
end

% nrms
ax(2) = subplot(10,2,(half*20+1):2:20);
imagesc(ax(2),query_DTT,1:(size(nrms,1)-N1),nrms(N1+1:end,:)); colormap(ax(2), 'summer');caxis(ax(2),[0 3])
% imagesc(ax(2),1,1:(size(nrms,1)-N1),max_r(N1+1:end)'); colormap(ax(2), ea_redblue); caxis([0 0.5]); %colorbar('northoutside')
% imagesc(ax(2),1,1:(size(nrms,1)-N1), 0.1*abs(lag_at_max(N1+1:end)')); colormap(ax(2), ea_redblue); caxis([0 2]); %colorbar('northoutside')
set(gca,'YTickLabel',{})
set(gca,'XTickLabel',{})

% stn entry exit
for j = 1:(size(nrms,1)-N1)
    arrayfun(@(x) line(ax(2), 'XData', [x x], 'YData',[-0.5 0.5]+j, 'LineWidth', 2), stn_entry_exit{j+N1});
end

%% Same but sorting with lag value

% sort top with lag + max_r
lag_and_r = lag_at_max;
max_r_refactor = max_r(lag_at_max == 0);
max_r_refactor = (max_r_refactor - min(max_r_refactor)) ./ (max(max_r_refactor) - min(max_r_refactor));
max_r_refactor = 1 ./ (1 + max_r_refactor);
max_r_refactor(1:2:end) = -max_r_refactor(1:2:end);
lag_and_r(lag_at_max == 0) = max_r_refactor;

[~,IDX] = sort(lag_and_r(1:N1));
lag_at_max2 = lag_at_max(IDX); 
nrms2 = nrms(IDX,:);
stn_entry_exit2 = stn_entry_exit(IDX);
brsh2 = brsh(IDX);
nrms_proc2 = nrms_proc(IDX,:);
stnd_proc2 = stnd_proc(IDX,:);

% nrms
ax(3) = subplot(10,2,2:2:(half*20));
imagesc(ax(3),query_DTT,1:size(nrms2,1),nrms2); colormap(ax(3), 'summer'); caxis(ax(3),[0 3])
% imagesc(ax(3),1,1:size(nrms2,1),0.1*abs(lag_at_max2)'); colormap(ax(3), ea_redblue); caxis([0 2]);% colorbar('northoutside')

% stn entry exit
for j = 1:length(stn_entry_exit2)
    arrayfun(@(x) line(ax(3), 'XData', [x x], 'YData',[-0.5 0.5]+j, 'LineWidth', 2), stn_entry_exit2{j});
end

% link
linkaxes(ax,'x');
xlim([-4 10])

%% brain shift analysis
lag_at_max2 = abs(0.1*lag_at_max2); % rescale to physical units
N2 = round(sum(lag_at_max2 > std(lag_at_max2)) / 2);

best_match = brsh2((-(N2-1):N2)+round(size(nrms2,1)/2));
worst_match = [brsh2(1:N2) brsh2(end-(N2-1):end)];

line(ax(3), 'XData', [-4 10], 'YData', [round(size(nrms2,1)/2)-N2 round(size(nrms2,1)/2)-N2], 'Color', [ 0    0.4470    0.7410], 'LineWidth', 2)
line(ax(3), 'XData', [-4 10], 'YData', [round(size(nrms2,1)/2)+N2 round(size(nrms2,1)/2)+N2], 'Color', [ 0    0.4470    0.7410], 'LineWidth', 2)
line(ax(3), 'XData', [-4 10], 'YData', [N2 N2], 'Color', [0.8500    0.3250    0.0980], 'LineWidth', 2)
line(ax(3), 'XData', [-4 10], 'YData', [size(nrms2,1)-N2 size(nrms2,1)-N2], 'Color', [0.8500    0.3250    0.0980], 'LineWidth', 2)

% ax(4)=subplot(10,2,(half*20+2):2:20); hold all % to plot in same figure

% rain clud
figure;
h1 = ea_raincloud_plot(best_match, 'box_on', 1, 'color', [ 0    0.4470    0.7410], 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .35,...
    'box_col_match', 1);
h2 = ea_raincloud_plot(worst_match, 'box_on', 1, 'color', [0.8500    0.3250    0.0980], 'alpha', 0.5,...
    'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .75,...
    'box_col_match', 1);
p = signrank(best_match, worst_match);
xlabel('Brainshift (mm)')
title(sprintf('Wilcoxon p = %.4f', p))
xlim([-0.2 1.2])

% corr_plot
[f2,R,p,g]=ea_corrplot(worst_match', [lag_at_max2(1:N2) lag_at_max2(end-(N2-1):end)]', {'', 'Brainshift (mm)', 'Lag (mm)'}, 'spearman');
% copyobj(allchild(f2.Children(3)), ax(4));
% xlabel(ax(4),{'brain shift', sprintf('R = %.3f p = %.3f', R, p)})
% ylabel(ax(4),'lag at max(xcorr)')
% close(f2)

%% Aux

function [displaced_index, values] = displace_row(row, displace_value)
    index = find(row);
    values = row(index);
    displaced_index = index - displace_value;
    out_of_bounds = displaced_index<1 | displaced_index>length(row);
    values(out_of_bounds) = [];
    displaced_index(out_of_bounds) = [];

end
