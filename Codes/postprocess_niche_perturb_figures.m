%% FIGURES FOR NICHE AND MORTALITY PERTURBATION
clear all; close all; clc
load('../Data/prep.mat')
foraging_trait = linspace(0,1,11);

%% NICHE SPEED (iteration each 101 time steps)
increment_loop = [0,linspace(.05,.2,20)]; % increment of perturbation
tol_a_t = [zeros(101,1);prep{1, 1}.niche_speed.tolerance.biomass.tolerance_a];
tol_a_f = [zeros(101,1);prep{1, 2}.niche_speed.tolerance.biomass.tolerance_a];
animal_density_t = [prep{1, 1}.niche_speed.animal_density;prep{1, 1}.niche_speed.animal_density_perturb];
animal_density_f = [prep{1, 2}.niche_speed.animal_density;prep{1, 2}.niche_speed.animal_density_perturb];

window = 20;
idx_t = 101:101:length(tol_a_t);
tol_a_t_end = [];
for i = 1:length(idx_t)
    tol_a = tol_a_t(idx_t(i)-window:idx_t(i));
    tol_a_t_end = [tol_a_t_end; mean(tol_a)];
end

window = 20;
idx_f = 101:101:length(tol_a_f);
tol_a_f_end = [];
for i = 1:length(idx_f)
    tol_a = tol_a_f(idx_f(i)-window:idx_f(i));
    tol_a_f_end = [tol_a_f_end; mean(tol_a)];
end

% compute z_mean for f communities
Z_mean = [];
animal_f_end = animal_density_f(idx_f,:);
for i = 1:size(animal_f_end,1)
    a = animal_f_end(i,:);
    a = reshape(a,77,11);
    a_z = sum(a,1);
    a_norm = a_z./sum(a_z);
    z_mean = sum(a_norm.*foraging_trait);
    Z_mean = [Z_mean, z_mean];
end

figure(1)
clf
hold on
l = line([increment_loop(1),increment_loop(length(idx_t))+0.01],[0,0],'color','k','LineWidth',1);
set(get(get(l,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
map = jet(256);
map(1,:) = [0,0,0];
Z_mean(isnan(Z_mean)) = -0.01;
s = scatter(increment_loop(1:length(idx_f)),tol_a_f_end,150,Z_mean,'filled','s');
set(get(get(s,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

s = scatter(increment_loop(1:length(idx_t)),tol_a_t_end,60,[zeros(length(idx_t)-1,1);-0.01],'filled');
set(get(get(s,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
colormap(map);
c = colorbar;
ylabel(c,'mean foraging trait','interpreter','latex','FontSize',20)
plot(increment_loop(1:length(idx_f)),tol_a_f_end,'--k','LineWidth',2)
plot(increment_loop(1:length(idx_t)),tol_a_t_end,'-k','LineWidth',2)

xlabel('pace of environmental change $c$','interpreter','latex','FontSize',20)
ylabel({'proportion of';'consumer biomass lost'},'interpreter','latex','FontSize',20)
xlim([increment_loop(1),increment_loop(length(idx_t))+0.01])
xlim([0,0.13])
caxis([0,1])
% legend('community with AF evolution','community without AF evolution','interpreter','latex','FontSize',15)
hold off


%% MORTALITY INCREASE 
% data stop at the last iteration before extinction
% we add a final data for each type of commu (t and f) corresponding to the
% mortality leading to extinction (to be coherent with the niche speed
% where we plot the data until the data of extinction)

%%% Threshold
increase_morta_f = [0.1;prep{1, 2}.press_mortality_threshold.animal_growth_dynamic];
increase_morta_f = [increase_morta_f;increase_morta_f(end) + increase_morta_f(end)*0.10]-0.1; % add the last mortality leading to extinction
increase_morta_t = prep{1, 1}.press_mortality_threshold.animal_growth_dynamic;
increase_morta_t = [0.1;increase_morta_t;increase_morta_t(end) + increase_morta_t(end)*0.10]-0.1; % add the last mortality leading to extinction
tol_a_t = prep{1, 1}.press_mortality_threshold.tolerance.biomass.tolerance_a;
tol_a_t = [0;tol_a_t;-1];
tol_a_f = prep{1, 2}.press_mortality_threshold.tolerance.biomass.tolerance_a;
tol_a_f = [0;tol_a_f;-1];
animal_density_t = prep{1, 1}.press_mortality_threshold.animal_density_perturb;
animal_density_f = prep{1, 2}.press_mortality_threshold.animal_density_perturb;

% %%% Increase
% increase_morta_f = prep{1, 2}.press_mortality_increase.animal_growth_dynamic;
% increase_morta_f = [increase_morta_f;increase_morta_f(end) + increase_morta_f(end)*0.10]; % add the last mortality leading to extinction
% increase_morta_t = prep{1, 1}.press_mortality_increase.animal_growth_dynamic;
% increase_morta_t = [increase_morta_t;increase_morta_f(end-1)]; % add the last mortality leading to extinction
% tol_a_t = prep{1, 1}.press_mortality_increase.tolerance.biomass.tolerance_a;
% tol_a_t = [tol_a_t;-1];
% tol_a_f = prep{1, 2}.press_mortality_increase.tolerance.biomass.tolerance_a;
% tol_a_f = [tol_a_f;-1];
% animal_density_t = prep{1, 1}.press_mortality_increase.animal_density_perturb;
% animal_density_f = prep{1, 2}.press_mortality_increase.animal_density_perturb;


% window = 20;
% idx_t = 101:101:length(tol_a_t);
% tol_a_t_end = [];
% for i = 1:length(idx_t)
%     tol_a = tol_a_t(idx_t(i)-window:idx_t(i));
%     tol_a_t_end = [tol_a_t_end; mean(tol_a)];
% end
% 
% window = 20;
% idx_f = 101:101:length(tol_a_f);
% tol_a_f_end = [];
% for i = 1:length(idx_f)
%     tol_a = tol_a_f(idx_f(i)-window:idx_f(i));
%     tol_a_f_end = [tol_a_f_end; mean(tol_a)];
% end

% compute z_mean for f communities
load('../Data/morta_with_timelaps.mat')
animal_density_f = [animal_density_perturb(1,:);animal_density_perturb(101:101:end,:)];
Z_mean = [];
for i =  1:size(animal_density_f,1)
    a = animal_density_f(i,:);
    a = reshape(a,11,11);
    a_z = sum(a,1);
    a_norm = a_z./sum(a_z);
    z_mean = sum(a_norm.*foraging_trait);
    Z_mean = [Z_mean, z_mean];
end
%
map = jet(256);
map(1,:) = [0,0,0];
Z_mean = [Z_mean, -0.01];

% data goes until extinction of communauty;
% to plot only until last sustainable mortality value;
% it must take the data until end-1th value (1:end-1);
figure(2)
clf
hold on
s = scatter(increase_morta_f,tol_a_f,150,Z_mean,'filled','s');
% s = scatter(increase_morta_f(1:end-1),tol_a_f(1:end-1),150,Z_mean,'filled','s');
set(get(get(s,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

s = scatter(increase_morta_t,tol_a_t,100,[zeros(length(increase_morta_t)-1,1);-0.01],'filled');
% s = scatter(increase_morta_t(1:end-1),tol_a_t(1:end-1),100,[zeros(length(increase_morta_t(1:end-1))-1,1);-0.01],'filled');
set(get(get(s,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
s.MarkerFaceAlpha = .7;
colormap(map);
% colormap jet
c = colorbar;
ylabel(c,'mean foraging trait','interpreter','latex','FontSize',20)
plot(increase_morta_f,tol_a_f,'--k','LineWidth',2)
plot(increase_morta_t,tol_a_t,'-k','LineWidth',2)
% plot(increase_morta_f(1:end-1),tol_a_f(1:end-1),'--k','LineWidth',2)
% plot(increase_morta_t(1:end-1),tol_a_t(1:end-1),'-k','LineWidth',2)
p = plot([increase_morta_f(1),increase_morta_f(end)],[0,0],'-k','LineWidth',.5);
set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
ylabel({'proportion of';'consumer biomass lost'},'interpreter','latex','FontSize',20)
xlabel('Additional consumers intrinsic mortality $\Delta d$','interpreter','latex','FontSize',20)
% legend('community with AF evolution','community without AF evolution','interpreter','latex','FontSize',15)
xlim([0 .6])
caxis([0,1])
hold off

%% NICHE JUMP (iteration each 101 time steps)
y0_dyn_t = 25+[-25;prep{1,1}.niche_displacement_range.y0_dyn];
y0_dyn_f = 25+[-25;prep{1,2}.niche_displacement_range.y0_dyn];

tol_a_t = [zeros(101,1);prep{1, 1}.niche_displacement_range.tolerance.biomass.tolerance_a];
tol_a_f = [zeros(101,1);prep{1, 2}.niche_displacement_range.tolerance.biomass.tolerance_a];
animal_density_t = [prep{1, 1}.niche_displacement_range.animal_density;prep{1, 1}.niche_displacement_range.animal_density_perturb];
animal_density_f = [prep{1, 2}.niche_displacement_range.animal_density;prep{1, 2}.niche_displacement_range.animal_density_perturb];

window = 20;
idx_t = 101:101:length(tol_a_t);
tol_a_t_end = [];
for i = 1:length(idx_t)
    tol_a = tol_a_t(idx_t(i)-window:idx_t(i));
    tol_a_t_end = [tol_a_t_end; mean(tol_a)];
end

window = 20;
idx_f = 101:101:length(tol_a_f);
tol_a_f_end = [];
for i = 1:length(idx_f)
    tol_a = tol_a_f(idx_f(i)-window:idx_f(i));
    tol_a_f_end = [tol_a_f_end; mean(tol_a)];
end

% compute z_mean for f communities
Z_mean = [];
animal_f_end = animal_density_f(idx_f,:);
for i = 1:size(animal_f_end,1)
    a = animal_f_end(i,:);
    a = reshape(a,77,11);
    a_z = sum(a,1);
    a_norm = a_z./sum(a_z);
    z_mean = sum(a_norm.*foraging_trait);
    Z_mean = [Z_mean, z_mean];
end

figure(1)
clf
hold on
l = line([y0_dyn_f(1),y0_dyn_f(end)],[0,0],'color','k','LineWidth',1);
set(get(get(l,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
map = jet(256);
map(1,:) = [0,0,0];
Z_mean(isnan(Z_mean)) = -0.01;
s = scatter(y0_dyn_f,tol_a_f_end,150,Z_mean,'filled','s');
set(get(get(s,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

s = scatter(y0_dyn_t,tol_a_t_end,60,[zeros(length(idx_t)-1,1);-0.01],'filled');
set(get(get(s,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
colormap(map);
c = colorbar;
ylabel(c,'mean foraging trait','interpreter','latex','FontSize',20)
plot(y0_dyn_f,tol_a_f_end,'--k','LineWidth',2)
plot(y0_dyn_t,tol_a_t_end,'-k','LineWidth',2)

xlabel('instantaneous niche shift $\Delta y$','interpreter','latex','FontSize',20)
ylabel({'proportion of';'consumer biomass lost'},'interpreter','latex','FontSize',20)
% xlim([increment_loop(1),increment_loop(length(idx_t))+0.01])
% xlim([0,0.13])
caxis([0,1])
% legend('community with AF evolution','community without AF evolution','interpreter','latex','FontSize',15)
hold off


%% NICHE SPEED INDIVIDUAL EFFECT OF Z
% load('niche_speed_with_fixed_z.mat')
load('../Data/niche_speed_with_fixed_z_101.mat')
increment_loop = linspace(.05,.2,20);
% foraging_trait = linspace(0,1,11);
foraging_trait = linspace(0,1,101);
% tspan = 101 time-steps per iteration :
for i = 1:length(y0_dyn_all)
    nbr_iteration_for_each_z(i) = length(y0_dyn_all{i})/101;
end
for i = 1:length(nbr_iteration_for_each_z)
    increment_for_each_z(i) = increment_loop(nbr_iteration_for_each_z(i));
end
plot(foraging_trait,increment_for_each_z,'ok','MarkerFaceColor','k','MarkerSize',10)
% scatter(foraging_trait,increment_for_each_z,'filled')
xlabel('foraging trait $z$','interpreter','latex','FontSize',20)
ylabel({'maximum sustainable pace'; 'of environmental change'},'interpreter','latex','FontSize',20)

%% MORTALITY INDIVIDUAL EFFECT OF Z
% load('../Data/morta_increase_z_fixed.mat')
load('../Data/increase_morta_z_fixed_101.mat')
foraging_trait = linspace(0,1,101);
morta_z = increase_morta_f;
morta_increase_z_fixed = morta_z(count_all+1);
% plot(foraging_trait,morta_increase_z_fixed,'--ok','MarkerFaceColor','k','MarkerSize',10)
figure(5)
clf
plot(foraging_trait,morta_increase_z_fixed,'--ok','MarkerFaceColor','k','MarkerSize',10)
xlabel('foraging trait $z$','interpreter','latex','FontSize',20)
ylabel({'maximum sustainable'; 'mortality increase'},'interpreter','latex','FontSize',20)


%% NICHE JUMP INDIVIDUAL EFFET OF Z
% load('jump_with_fixed_z.mat')
load('../Data/jump_with_fixed_z_101.mat')

jump_size = [];
for i = 1:length(y0_dyn_all)
    a = y0_dyn_all{i};
    jump_size = [jump_size,length(a)];
end
foraging_trait = linspace(0,1,101);
plot(foraging_trait,jump_size,'ok','MarkerFaceColor','k','MarkerSize',10)
% scatter(foraging_trait,jump_size,'filled','b')
xlabel('foraging trait $z$','interpreter','latex','FontSize',20)
ylabel({'maximum sustainable' ;'niche shift'},'interpreter','latex','FontSize',20)



