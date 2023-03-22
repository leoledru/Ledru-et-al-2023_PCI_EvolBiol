%% CODE POUR ANALYSE PRCC ET PRODUCTION DE FIGURES

clear all; close all;
% load('PRCC.mat');
% load('PRCC_b.mat');
% load('PRCC_c.mat');
% load('PRCC_all.mat');
load('../Data/PRCC_e.mat');
foraging_trait = linspace(0,1,11);

%
% Partial Rank correlation coefficient %

% Initialization Outcomes %
parameters_all = [];
for i = 1:length(Output_a)
    params = Param{i};
    parameters_all = [parameters_all;params];
end

proxy = []; % average z
A = [];
for i = 1:length(Output_a)
% for i = 1:20
    a = Output_a{i};
    
    % take the 11 time steps
%     for j = 1:size(test,1)
%         b = a(j,:);
%         b = reshape(b,11,11);
%         b_sum = sum(b,1);
%         A = [A;b_sum];
%     end
%     A_sum = sum(A,1);

    % takes only the last time step
    a = a(end,:);
    a = reshape(a,11,11);
    a_sum = sum(a,1);
    a_norm = a_sum./sum(a_sum);
    proxy = [proxy;sum(a_norm.*foraging_trait)];
end
% selon les params il est possible que des simus vont a l'extinction
% dans ce cas retirer ces simus de l'analyse
idx_crash = find(isnan(proxy));
proxy(idx_crash) = [];
parameters_all(idx_crash,:) = [];
Output_a(idx_crash) = [];

nbr_parameters = size(parameters_all,2);
spearman_r = zeros(nbr_parameters,1);
spearman_p_value = zeros(nbr_parameters,1);

for is = 1:nbr_parameters
    param_new = parameters_all;
    param_new(:,is) = []; % remove the focus parameter 
    Z = param_new;
    % Partial correlation
    Wx = Z\parameters_all(:,is); % compare the focus parameter with all others (to remove the linear correlation)
    Wy = Z\proxy; % compare the focus parameter with the proxy to remove the linear correlation
    Res_x = parameters_all(:,is) - Z*Wx;
    Res_y = proxy - Z*Wy;
    % spearman
    [rhoi,pi] = spearman(Res_x,Res_y,0.05);
    spearman_r(is) = rhoi;
    spearman_p_value(is) = pi;
end


%% Figure Sensitivity
% figure
% clf

% subplot(1,2,1)
plot(1:nbr_parameters,spearman_r,'ko','MarkerFaceColor','k','MarkerSize',8)
xticks(1:6);
% sigma, sigmaK, hmax_loop, animal_intrinsic_growth, animal_intcraspe_compet, plant_intrinsic_growth;
xticklabels({'$\sigma$','$\sigma_K$','$s_{max}$','$d$','$I_c$','$g$'});
set(gca,'TickLabelInterpreter','latex')
set(gca,'TickLabelInterpreter','latex')
ax = gca;
ax.XAxis.FontSize = 18;
ax.YAxis.FontSize = 12;
ylabel('Correlation Coefficient','interpreter','latex','FontSize',20)

hold on
plot(0:nbr_parameters+1,zeros(nbr_parameters+2,1),'k--')
% plot([.8 1 1.2],[spearman_r(1)+.04 spearman_r(1)+.04 spearman_r(1)+.04],'k*')
% plot([1.8 2 2.2],[spearman_r(2)+.04 spearman_r(2)+.04 spearman_r(2)+.04],'k*')
% plot([2.8 3 3.2],[spearman_r(3)+.04 spearman_r(3)+.04 spearman_r(3)+.04],'k*')
% plot([3.9 4.1],[spearman_r(4)+.04 spearman_r(4)+.04],'k*')
% plot([4.9 5.1],[spearman_r(5)+.04 spearman_r(5)+.04],'k*')
% plot([5.8 6 6.2],[spearman_r(6)+.04 spearman_r(6)+.04 spearman_r(6)+.04],'k*')

xlim([0 nbr_parameters+1]);
ylim([-.8 .5])
hold off
% 
% subplot(1,2,2)
% plot(1:nbr_parameters,spearman_p_value,'ko','MarkerFaceColor','k','MarkerSize',8)
% xticks(1:6);
% % sigma, sigmaK, hmax_loop, animal_intrinsic_growth;
% xticklabels({'$\sigma$','$\sigma_K$','$s_{max}$','$d$','$I_c$','$g$'});
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% ax = gca;
% ax.XAxis.FontSize = 18;
% ax.YAxis.FontSize = 12;
% ylabel('p-value','interpreter','latex','FontSize',20)
% hold on
% plot(0:nbr_parameters+1,0.05.*ones(nbr_parameters+2,1),'k--')
% xlim([0 nbr_parameters+1]);
% hold off

%% Z average
% boxplot(proxy)

idx_min = find(proxy==min(proxy));
idx_max = find(proxy==max(proxy));
target = 0.35;
[closest_target,idx_target] = min(abs(proxy - target)); 
param_min = parameters_all(idx_min,:);
param_max = parameters_all(idx_max,:);
param_target = parameters_all(idx_target,:);

% plot(param_min,'bo','MarkerFaceColor','b')
% hold on
% plot(param_max,'ro','MarkerFaceColor','r')
% hold off
% xticks(1:4);
% % sigma, sigmaK, hmax_loop, animal_intrinsic_growth;
% xticklabels({'$\sigma$','$\sigma_K$','$s_{max}$','$d$'});
% set(gca,'TickLabelInterpreter','latex')
% set(gca,'TickLabelInterpreter','latex')
% ax = gca;
% ax.XAxis.FontSize = 18;

% parameters_all
% sigma, sigmaK, hmax_loop, animal_intrinsic_growth, animal_intcraspe_compet, plant_intrinsic_growth;

high_af = Output_a{idx_max};
high_af = high_af(end,:);
high_af = reshape(high_af,11,11);
low_af = Output_a{idx_min};
low_af = low_af(end,:);
low_af = reshape(low_af,11,11);
intermediate_af = Output_a{idx_target};
intermediate_af = intermediate_af(end,:);
intermediate_af = reshape(intermediate_af,11,11);
colormap pink
subplot(1,3,1)
imagesc(low_af)
set(gca,'YTick',1:1:11,'YTickLabel',-5:5)
set(gca,'XTick',1:2:11,'XTickLabel',0:.2:1)
xlabel('foraging trait','interpreter','latex','Fontsize',20)
ylabel('niche position trait','interpreter','latex','Fontsize',20)
title({['$\sigma = \\$' num2str(round(parameters_all(idx_min,1),3))],...
    ['$\sigma_K = \\$' num2str(round(parameters_all(idx_min,2),3))],...
    ['$s_{max} = \\$' num2str(round(parameters_all(idx_min,3),3))],...
    ['$d = \\$' num2str(round(parameters_all(idx_min,4),3))],...
    ['$I_c = \\$' num2str(round(parameters_all(idx_min,5),3))],...
    ['$g = \\$' num2str(round(parameters_all(idx_min,6),3))]},'interpreter','latex','FontSize',15);
c = colorbar;
ylabel(c,'abundance','interpreter','latex','FontSize',20);

subplot(1,3,2)
imagesc(intermediate_af)
set(gca,'YTick',1:1:11,'YTickLabel',-5:5)
set(gca,'XTick',1:2:11,'XTickLabel',0:.2:1)
xlabel('foraging trait','interpreter','latex','Fontsize',20)
ylabel('niche position trait','interpreter','latex','Fontsize',20)
title({['$\sigma = \\$' num2str(round(parameters_all(idx_target,1),3))],...
    ['$\sigma_K = \\$' num2str(round(parameters_all(idx_target,2),3))],...
    ['$s_{max} = \\$' num2str(round(parameters_all(idx_target,3),3))],...
    ['$d = \\$' num2str(round(parameters_all(idx_target,4),3))],...
    ['$I_c = \\$' num2str(round(parameters_all(idx_target,5),3))],...
    ['$g = \\$' num2str(round(parameters_all(idx_target,6),3))]},'interpreter','latex','FontSize',15);
c = colorbar;
ylabel(c,'abundance','interpreter','latex','FontSize',20);

subplot(1,3,3)
imagesc(high_af)
set(gca,'YTick',1:1:11,'YTickLabel',-5:5)
set(gca,'XTick',1:2:11,'XTickLabel',0:.2:1)
xlabel('foraging trait','interpreter','latex','Fontsize',20)
ylabel('niche position trait','interpreter','latex','Fontsize',20)
title({['$\sigma = \\$' num2str(round(parameters_all(idx_max,1),3))],...
    ['$\sigma_K = \\$' num2str(round(parameters_all(idx_max,2),3))],...
    ['$s_{max} = \\$' num2str(round(parameters_all(idx_max,3),3))],...
    ['$d = \\$' num2str(round(parameters_all(idx_max,4),3))],...
    ['$I_c = \\$' num2str(round(parameters_all(idx_max,5),3))],...
    ['$g = \\$' num2str(round(parameters_all(idx_max,6),3))]},'interpreter','latex','FontSize',15);
c = colorbar;
ylabel(c,'abundance','interpreter','latex','FontSize',20);


%%


