%% POSTPROCESS DE LA COMPARAISON ENTRE COMMUNAUTES AVEC/SANS AF

clear all; clc;
close all


Color = get(gca,'colororder');
choice = 2; % 1 = MEAN +CI 5%
            % 2 = MEDIAN + QUANTILE 25%-75%

number_of_animals = 11;
number_of_plants  = 11;
number_of_foraging = 11;

xmin = -5;
xmax = 5;
xx = linspace(xmin,xmax,number_of_animals);
dx = xx(2)-xx(1);
dx = 1;
traits_of_animals = xx';
foraging_trait = linspace(0,1,number_of_foraging);

% load('community_comparison6.mat');
% load('community_comparison_CK.mat');

% load('community_comparison.mat');
% load('community_comparison_b.mat');
load('../Data/community_comparison_all.mat'); % community_comparison.mat + community_comparison_b.mat

% each column is a simulation, row 1 = animal simu AF; row 2 = plant simu AF;
% row 3 = Effort simu AF; row 4 = animal simu tondeuse; row 5 = plant simu
% tondeuse; row 6 = Effort tondeuse; row 7 = params

% postprocess = cell(7,size(Param,2));
% load('postprocess_densities')

for ii = 1:size(Param,2)
    % for ii = 1:1
    postprocess{1,ii} = Output_a_f(:,(ii-1)*number_of_animals*number_of_foraging+1:ii*number_of_animals*number_of_foraging);
    postprocess{2,ii} = Output_p_f(:,(ii-1)*number_of_plants+1:ii*number_of_plants);
%     postprocess{3,ii} = Effort_f(:,(ii-1)*number_of_animals*number_of_foraging+1:ii*number_of_animals*number_of_foraging*number_of_plants);
    postprocess{4,ii} = Output_a_t(:,(ii-1)*number_of_animals*number_of_foraging+1:ii*number_of_animals*number_of_foraging);
    postprocess{5,ii} = Output_p_t(:,(ii-1)*number_of_plants+1:ii*number_of_plants);
    %   postprocess{6,ii} = Effort_t(:,(ii-1)*number_of_animals*number_of_foraging+1:ii*number_of_animals*number_of_foraging*number_of_plants);
    postprocess{7,ii} = Param(:,ii);
end


%% Compute this section for agregation of densities over z axis (necessary
%% to further compute the diversity indices which only consider the
%% diversity along x axis).

% Reshape animal density (postprocess rows 1 and 4)
Out_f = cell(1,size(postprocess,2));
Out_t = cell(1,size(postprocess,2));
Out_F = cell(1,size(postprocess,2));
Out_T = cell(1,size(postprocess,2));
for ii = 1:size(postprocess,2)
    Output_a_f = zeros(101*11,number_of_plants);
    Output_a_t = zeros(101*11,number_of_plants);
    output_a_f = postprocess{1,ii};
    output_a_t = postprocess{4,ii};
    for jj = 1:size(output_a_f,1)
        Output_a_f(11*(jj-1)+1:11*jj,:) = reshape(output_a_f(jj,:),number_of_animals,number_of_foraging);
        Output_a_t(11*(jj-1)+1:11*jj,:) = reshape(output_a_t(jj,:),number_of_animals,number_of_foraging);
    end
    Out_F{1,ii} = Output_a_f;
    Out_T{1,ii} = Output_a_t;
end
%
% Sum animal densities on foraging trait

for ii = 1:size(postprocess,2)
    animal_niche_f = zeros(101,number_of_plants);
    animal_f = Out_F{:,ii};
    animal_niche_t = zeros(101,number_of_plants);
    animal_t = Out_T{:,ii};
    for jj = 1:101
        animal_niche_f(jj,:) = sum(animal_f((jj-1)*number_of_foraging+1:jj*number_of_foraging,:),2)';
        animal_niche_t(jj,:) = sum(animal_t((jj-1)*number_of_foraging+1:jj*number_of_foraging,:),2)';
    end
    Out_f{:,ii} = animal_niche_f;
    Out_t{:,ii} = animal_niche_t;
end

%%
% BOXPLOT (average on the last 100 time steps)

biomass_animal_f = zeros(size(Param,2),size(postprocess{1,1},1));
biomass_animal_t = zeros(size(Param,2),size(postprocess{1,1},1));
biomass_plant_f  = zeros(size(Param,2),size(postprocess{1,1},1));
biomass_plant_t  = zeros(size(Param,2),size(postprocess{1,1},1));
for ii = 1:size(Param,2)
    biomass_animal_f(ii,:) = sum(postprocess{1,ii},2)';
    biomass_animal_t(ii,:) = sum(postprocess{4,ii},2)';
    biomass_plant_f(ii,:)  = sum(postprocess{2,ii},2)';
    biomass_plant_t(ii,:)  = sum(postprocess{5,ii},2)';
end

Biomass_plant_f = mean(biomass_plant_f(:,end-10:end),2);
Biomass_plant_t = mean(biomass_plant_t(:,end-10:end),2);
Biomass_plant_boxplot = [Biomass_plant_f,Biomass_plant_t];

Biomass_animal_f = mean(biomass_animal_f(:,end-10:end),2);
Biomass_animal_t = mean(biomass_animal_t(:,end-10:end),2);
Biomass_animal_boxplot = [Biomass_animal_f,Biomass_animal_t];

figure(10)
clf
subplot(1,2,1)
boxplot(Biomass_animal_boxplot)
title('Total consumers biomass','interpreter','latex','FontSize',18)
xticks([1 2 3 4])
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTickLabel', {'with foraging','without'});
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)

subplot(1,2,2)
boxplot(Biomass_plant_boxplot)
title('Total resources biomass','interpreter','latex','FontSize',18)
xticks([1 2 3 4])
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTickLabel', {'with foraging','without'});
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)

%% Compute z average for color in scatters below
Z_mean = zeros(1,size(postprocess,2));
Z_var  = zeros(1,size(postprocess,2));
for ii = 1:size(postprocess,2)
    animal = postprocess{1,ii};
    z_mean = [];
    z_var =  [];
    for i = 0:10
        a = animal(end-i,:);
        a_reshape = reshape(a,number_of_animals,number_of_foraging);
        a_sum_and_norm = sum(a_reshape,1)./sum(sum(a_reshape,1));
        z_m    = sum(a_sum_and_norm.*foraging_trait);
        z_mean = [z_mean,z_m];
        z_var  = [z_var,sum(a_sum_and_norm.*(foraging_trait-z_m).^2)];
    end
    Z_mean(ii) = mean(z_mean);
    Z_var(ii)  = mean(z_var);
end


%% VARIANTE PLUS CONCISE : PLOT DE LA DISTANCE ENTRE BIOMASSE_TONDEUSE ET BIOMASSE_AF EN FONCTION DU Z_MOYEN
figure(1)
clf
hold on
l = line([1,0],[0,0],'color','k','linewidth',2,'linestyle','-');
set(get(get(l,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

%%% Absolute differnece
% dist_biomasse_animal = Biomass_animal_boxplot(:,1) - Biomass_animal_boxplot(:,2);
%%% Percentage of differnece rate
dist_biomasse_animal = (Biomass_animal_boxplot(:,1) - Biomass_animal_boxplot(:,2))./Biomass_animal_boxplot(:,2)*100;

% yyaxis left
s = scatter(Z_mean,dist_biomasse_animal,5,Color(1,:),'filled');
s.MarkerFaceAlpha = .5;
%%% Absolute differnece
% dist_biomasse_resource = Biomass_plant_boxplot(:,1) - Biomass_plant_boxplot(:,2);
%%% Percentage of differnece rate
dist_biomasse_resource = (Biomass_plant_boxplot(:,1) - Biomass_plant_boxplot(:,2))./Biomass_plant_boxplot(:,2)*100;

% yyaxis right
s = scatter(Z_mean,dist_biomasse_resource,5,'d','filled','MarkerFaceColor',Color(5,:));
s.MarkerFaceAlpha = .5;

legend('consumers','resources','interpreter','latex','FontSize',15,'location','northwest')
% ax = gca;
% ax.YAxis(1).Color = Color(1,:); %'b';
% ax.YAxis(2).Color = Color(5,:); %[.2 .8 .2];

% AJOUT MEAN + SHADE STD
Z_mean_round = round(Z_mean,1); % round pour réduire le nombre de Z différents
% et avoir plusieurs communautés (réplicas)
% pour chaque Z
Z_unique = unique(Z_mean_round);
% RESOURCE
dist_biomasse_resource_replica = {};
for i = 1:length(Z_unique)
    value = Z_unique(i);
    dist_biomasse_resource_replica{i} = dist_biomasse_resource(Z_mean_round==value);
end
dist_biomasse_resource_median = cellfun(@median, dist_biomasse_resource_replica);
dist_biomasse_resource_quantile_inf = cellfun(@(x) quantile(x,0.25), dist_biomasse_resource_replica);
dist_biomasse_resource_quantile_sup = cellfun(@(x) quantile(x,0.75), dist_biomasse_resource_replica);
dist_biomasse_resource_mean = cellfun(@mean, dist_biomasse_resource_replica);
dist_biomasse_resource_std = cellfun(@std, dist_biomasse_resource_replica);
% CONSUMER
dist_biomasse_animal_replica = {};
for i = 1:length(Z_unique)
    value = Z_unique(i);
    dist_biomasse_animal_replica{i} = dist_biomasse_animal(Z_mean_round==value);
end
dist_biomasse_animal_median = cellfun(@median, dist_biomasse_animal_replica);
dist_biomasse_animal_quantile_inf = cellfun(@(x) quantile(x,0.25), dist_biomasse_animal_replica);
dist_biomasse_animal_quantile_sup = cellfun(@(x) quantile(x,0.75), dist_biomasse_animal_replica);

dist_biomasse_animal_mean = cellfun(@mean, dist_biomasse_animal_replica);
dist_biomasse_animal_std = cellfun(@std, dist_biomasse_animal_replica);

x = 0:.1:1;
x2 = [x, fliplr(x)];

dist_biomasse_resource_mean = movmean(dist_biomasse_resource_mean,2);
dist_biomasse_resource_std = movmean(dist_biomasse_resource_std,2);
% yyaxis right
hold on

if (choice == 1)
    pr = plot(Z_unique,dist_biomasse_resource_mean,'--','LineWidth',3,'Color',Color(5,:));
    inBetween = [dist_biomasse_resource_mean-dist_biomasse_resource_std,...
        fliplr(dist_biomasse_resource_mean+dist_biomasse_resource_std)];
    set(get(get(pr,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    h = fill(x2, inBetween, Color(5,:));
    set(h,'facealpha',.2,'EdgeColor',Color(5,:));
    ymax_r = max(dist_biomasse_resource_mean+dist_biomasse_resource_std)*1.1;
    ymin_r = min(dist_biomasse_resource_mean-dist_biomasse_resource_std)*1.1;
%     ylim([ymin,ymax])
elseif (choice == 2 )
    pr = plot(Z_unique,dist_biomasse_resource_median,'--','LineWidth',3,'Color',Color(5,:));
    set(get(get(pr,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    inBetween = [dist_biomasse_resource_quantile_inf,fliplr(dist_biomasse_resource_quantile_sup)];
    h = fill(x2, inBetween, Color(5,:));
    set(h,'facealpha',.2,'EdgeColor',Color(5,:));
    ymax_r = max(dist_biomasse_resource_quantile_sup)*1.1;
    ymin_r = min(dist_biomasse_resource_quantile_inf)*1.1;
%     ylim([ymin,ymax])
end

set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hold off

dist_biomasse_animal_mean = movmean(dist_biomasse_animal_mean,2);
dist_biomasse_animal_std = movmean(dist_biomasse_animal_std,2);
% yyaxis left
hold on

if (choice == 1)
    pa = plot(Z_unique,dist_biomasse_animal_mean,'--','color',Color(1,:),'LineWidth',3);
    set(get(get(pa,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    inBetween = [dist_biomasse_animal_mean-dist_biomasse_animal_std,...
        fliplr(dist_biomasse_animal_mean+dist_biomasse_animal_std)];
    h = fill(x2, inBetween, Color(1,:));
    set(h,'facealpha',.2,'EdgeColor',Color(1,:));
    ymax_a = max(dist_biomasse_animal_mean+dist_biomasse_animal_std)*1.1;
    ymin_a= min(dist_biomasse_animal_mean-dist_biomasse_animal_std)*1.1;
%     ylim([ymin,ymax])
elseif (choice ==2)
    pa = plot(Z_unique,dist_biomasse_animal_median,'--','LineWidth',3,'Color',Color(1,:));
    set(get(get(pa,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    inBetween = [dist_biomasse_animal_quantile_inf,fliplr(dist_biomasse_animal_quantile_sup)];
    h = fill(x2, inBetween, Color(1,:));
    set(h,'facealpha',.2,'EdgeColor',Color(1,:));
    ymax_a = max(dist_biomasse_animal_quantile_sup)*1.1;
    ymin_a = min(dist_biomasse_animal_quantile_inf)*1.1;
%     ylim([ymin,ymax])
end
ymin = min(ymin_a,ymin_r);
ymax = max(ymax_r,ymax_a);
ylim([ymin,ymax])
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hold off

ax = gca;
ax.FontSize = 16;
yticks(ax,[-75 -50 -25 0 25 50 75])
yticklabels(ax,{'-75','-50','-25','0','25','50','75'})
grid on
ax.GridLineStyle = '--';
xlabel({'mean foraging trait of','the system with AF evolution'},'interpreter','latex','FontSize',20)
ylabel({'Percent difference'; 'in biomass'},'interpreter','latex','FontSize',20,'color','k')
% yyaxis right
% ylabel({'difference in','resource biomass'},'interpreter','latex','FontSize',20)
% Make the axis align
% align_yyaxis_zero(ax)

%%

%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONAL DIVERSITY %
% Laliberté & Legendre 2010 %
%%%%%%%%%%%%%%%%%%%%%%%%

FDis_animal_f = zeros(101,size(Param,2));
FDis_animal_t = zeros(101,size(Param,2));
FDis_plant_f = zeros(101,size(Param,2));
FDis_plant_t = zeros(101,size(Param,2));
for ii = 1:size(postprocess,2)
    % Densities
    animal_density_f = Out_f{:,ii};
    plant_density_f = postprocess{2,ii};
    animal_density_t = Out_t{:,ii};
    plant_density_t = postprocess{5,ii};
    
    % Centroid
    C_animal_f = sum(animal_density_f.*xx,2)./sum(animal_density_f,2);
    C_animal_t = sum(animal_density_t.*xx,2)./sum(animal_density_t,2);
    C_plant_f = sum(plant_density_f.*xx,2)./sum(plant_density_f,2);
    C_plant_t = sum(plant_density_t.*xx,2)./sum(plant_density_t,2);
    
    % Functional dispersion
    z_animal_f = abs(xx-C_animal_f);
    FDis_animal_f(:,ii) = sum(animal_density_f.*z_animal_f,2)./sum(animal_density_f,2);
    z_animal_t = abs(xx-C_animal_t);
    FDis_animal_t(:,ii) = sum(animal_density_t.*z_animal_t,2)./sum(animal_density_t,2);
    z_plant_f = abs(xx-C_plant_f);
    FDis_plant_f(:,ii) = sum(plant_density_f.*z_plant_f,2)./sum(plant_density_f,2);
    z_plant_t = abs(xx-C_plant_t);
    FDis_plant_t(:,ii) = sum(plant_density_t.*z_plant_t,2)./sum(plant_density_t,2);
end

%%
%%%%%%%%%%%
% Boxplot %
% average on the last 100 time steps
FDis_animal_boxplot = [mean(FDis_animal_f(end-10:end,:),1);mean(FDis_animal_t(end-10:end,:),1)];
FDis_plant_boxplot = [mean(FDis_plant_f(end-10:end,:),1);mean(FDis_plant_t(end-10:end,:),1)];

figure(11)
clf
subplot(1,2,1)
boxplot(FDis_animal_boxplot')
title('Functional dispersion of consumers','interpreter','latex','FontSize',18)
xticks([1 2 3 4])
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTickLabel', {'with foraging','without'});
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)

subplot(1,2,2)
boxplot(FDis_plant_boxplot')
title('Functional dispersion of resources','interpreter','latex','FontSize',18)
xticks([1 2 3 4])
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTickLabel', {'with foraging','without'});
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)

%%  Compute z average for color in scatters below
Z_mean = [];
for ii = 1:size(postprocess,2)
    animal = postprocess{1,ii};
    z_mean = [];
    for i = 0:10
        a = animal(end-i,:);
        a_reshape = reshape(a,number_of_animals,number_of_foraging);
        a_sum_and_norm = sum(a_reshape,1)./sum(sum(a_reshape,1));
        z_mean = [z_mean,sum(a_sum_and_norm.*foraging_trait)];
    end
    Z_mean = [Z_mean, mean(z_mean)];
end

%% VARIANTE PLUS CONCISE : PLOT DE LA DISTANCE ENTRE BIOMASSE_TONDEUSE ET BIOMASSE_AF EN FONCTION DU Z_MOYEN
figure(2)
clf
hold on
l = line([1,0],[0,0],'color','k','linewidth',2,'linestyle','-');
set(get(get(l,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

%%% Absolute difference
% dist_FDis_animal = FDis_animal_boxplot(1,:) - FDis_animal_boxplot(2,:);
%%% Relative difference
dist_FDis_animal = (FDis_animal_boxplot(1,:) - FDis_animal_boxplot(2,:))./FDis_animal_boxplot(2,:)*100;

s = scatter(Z_mean,dist_FDis_animal,5,Color(1,:),'filled');
s.MarkerFaceAlpha = .5;
hold on

%%% Absolute differnece
% dist_FDis_resource = FDis_plant_boxplot(1,:) - FDis_plant_boxplot(2,:);
%%% Relative difference
dist_FDis_resource = (FDis_plant_boxplot(1,:) - FDis_plant_boxplot(2,:) )./FDis_plant_boxplot(2,:)*100;

s = scatter(Z_mean,dist_FDis_resource,5,'d','filled','MarkerFaceColor',Color(5,:));
s.MarkerFaceAlpha = .5;

legend('consumers','resources','interpreter','latex','FontSize',15,'location','northwest')
xlabel({'mean foraging trait','of community with AF evolution'},'interpreter','latex','FontSize',20)

% AJOUT MEAN + SHADE STD
Z_mean_round = round(Z_mean,1); % round pour réduire le nombre de Z différents
% et avoir plusieurs communautés (réplicas)
% pour chaque Z
Z_unique = unique(Z_mean_round);
% RESOURCE
dist_FDis_resource_replica = {};
for i = 1:length(Z_unique)
    value = Z_unique(i);
    dist_FDis_resource_replica{i} = dist_FDis_resource(Z_mean_round==value);
end
dist_FDis_resource_median = cellfun(@median, dist_FDis_resource_replica);
dist_FDis_resource_quantile_inf = cellfun(@(x) quantile(x,0.25), dist_FDis_resource_replica);
dist_FDis_resource_quantile_sup = cellfun(@(x) quantile(x,0.75), dist_FDis_resource_replica);
dist_FDis_resource_mean = cellfun(@mean, dist_FDis_resource_replica);
dist_FDis_resource_std = cellfun(@std, dist_FDis_resource_replica);
% CONSUMER
dist_FDis_animal_replica = {};
for i = 1:length(Z_unique)
    value = Z_unique(i);
    dist_FDis_animal_replica{i} = dist_FDis_animal(Z_mean_round==value);
end
dist_FDis_animal_median = cellfun(@median, dist_FDis_animal_replica);
dist_FDis_animal_quantile_inf = cellfun(@(x) quantile(x,0.25), dist_FDis_animal_replica);
dist_FDis_animal_quantile_sup = cellfun(@(x) quantile(x,0.75), dist_FDis_animal_replica);
dist_FDis_animal_mean = cellfun(@mean, dist_FDis_animal_replica);
dist_FDis_animal_std = cellfun(@std, dist_FDis_animal_replica);

x = 0:.1:1;
x2 = [x, fliplr(x)];

dist_FDis_resource_mean = movmean(dist_FDis_resource_mean,2);
dist_FDis_resource_std = movmean(dist_FDis_resource_std,2);
hold on
if (choice==1)
    p = plot(Z_unique,dist_FDis_resource_mean,'--','LineWidth',3,'Color',Color(5,:));
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    inBetween = [dist_FDis_resource_mean-dist_FDis_resource_std,...
        fliplr(dist_FDis_resource_mean+dist_FDis_resource_std)];
    h = fill(x2, inBetween, Color(5,:));
    set(h,'facealpha',.2,'EdgeColor',Color(5,:));
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    ymax_r = max(dist_FDis_resource_mean+dist_FDis_resource_std)*1.1;
    ymin_r = min(dist_FDis_resource_mean-dist_FDis_resource_std)*1.1;
    
    dist_FDis_animal_mean = movmean(dist_FDis_animal_mean,2);
    dist_FDis_animal_std = movmean(dist_FDis_animal_std,2);
    p = plot(Z_unique,dist_FDis_animal_mean,'--','LineWidth',3,'Color',Color(1,:));
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    inBetween = [dist_FDis_animal_mean-dist_FDis_animal_std,...
        fliplr(dist_FDis_animal_mean+dist_FDis_animal_std)];
    h = fill(x2, inBetween, Color(1,:));
    set(h,'facealpha',.2,'EdgeColor',Color(1,:));
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    ymax_a = max(dist_FDis_animal_mean+dist_FDis_animal_std)*1.1;
    ymin_a = min(dist_FDis_animal_mean-dist_FDis_animal_std)*1.1;
    
elseif(choice==2)
    p = plot(Z_unique,dist_FDis_resource_median,'--','LineWidth',3,'Color',Color(5,:));
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    inBetween = [dist_FDis_resource_quantile_inf,...
        fliplr(dist_FDis_resource_quantile_sup)];
    h = fill(x2, inBetween, Color(5,:));
    set(h,'facealpha',.2,'EdgeColor',Color(5,:));
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    ymax_r = max(dist_FDis_resource_quantile_sup)*1.1;
    ymin_r = min(dist_FDis_resource_quantile_inf)*1.1;
    
    p = plot(Z_unique,dist_FDis_animal_median,'--','LineWidth',3,'Color',Color(1,:));
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    inBetween = [dist_FDis_animal_quantile_inf,...
        fliplr(dist_FDis_animal_quantile_sup)];
    h = fill(x2, inBetween, Color(1,:));
    set(h,'facealpha',.2,'EdgeColor',Color(1,:));
    set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    ymax_a = max(dist_FDis_animal_quantile_sup)*1.1;
    ymin_a = min(dist_FDis_animal_quantile_inf)*1.1;
end
ymin = min(ymin_r,ymin_a)-0.05;
ymax = max(ymax_r,ymax_a);
ylim([ymin,ymax])
hold off

ax = gca;
ax.FontSize = 16;
grid on
ax.GridLineStyle = '--';
xlabel({'mean foraging trait of','the system with AF evolution'},'interpreter','latex','FontSize',20)
ylabel({'Percent difference in','functional diversity'},'interpreter','latex','FontSize',20)
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PRODUCTIVITY %%%
% Sum of all consumers functional responses weigthed by the consumers abundances
% = flux from resources to consumers

load('../Data/community_comparison_all.mat'); % community_comparison.mat + community_comparison_b.mat

number_of_animals = 11;
number_of_plants  = 11;
number_of_foraging = 11;

% FONCTION K(y)
K0 = 50; % maximal carrying capacity
y0 = 0; % optimal trait
sigmaK = 2.5;
Kf = @(trait_plant) K0*exp(-(trait_plant - y0).^2./(2*sigmaK^2));

% FONCTION C(y-y0)
Beta = 0; % Beta = 0 : symetrical competition
sigmaC = sigmaK-1;
C = @(trait_plant_y,trait_plant_yi) exp(sigmaC^2*Beta^2/2).*exp(-(trait_plant_y-trait_plant_yi + sigmaC^2*Beta^2).^2./(2*sigmaC^2));

% TRAITS %
% Animal traits
xmin = -5;
xmax = 5;
xx = linspace(xmin,xmax,number_of_animals);
dx = xx(2)-xx(1);
traits_of_animals = xx';

foraging_trait = linspace(0,1,number_of_foraging);
dz = foraging_trait(2) - foraging_trait(1);
Foraging_trait = reshape(foraging_trait,1,1,number_of_foraging);

% Plants traits
ymin = -5;
ymax = 5;
yy = linspace(ymin,ymax,number_of_plants);
dy =  yy(2)-yy(1);
traits_of_plants = yy';

[XX,YY] = meshgrid(xx,yy);
SIGMA = .5;
delta_ij = complementary_traits(SIGMA,XX,YY);
delta_ij = delta_ij';

effort_speed_of_change = 1;

% Competition %
Compet = C(traits_of_plants,traits_of_plants');
% Carrying capacity %
K = Kf(traits_of_plants);

% HOMOGENEOUS INITIALISATION %
% r_j = .5;                      % .*ones(number_of_plants,1);
% handling_time     = .01;
extraction_coeff = .8;         % .*ones(number_of_animals,1);
conversion_coeff = .3;         % .*ones(number_of_animals,1);
animal_intrinsic_growth = .1; % .*ones(number_of_animals,1);
animal_intraspe_compet  = .01; % .*ones(number_of_animals,1);
plant_intrinsic_growth  = .8; % .5.*ones(number_of_plants,1);
seuil_abondance = 1e-5;
seuil_effort = seuil_abondance;

% HANDLING TIME TRADE-OFF
% alpha_h = 5; % expo
alpha_h = 1; % lineaire
% alpha_h = .5; % racine
hmin = .1;
hmax = 0.55;
h = @(z) hmin + (hmax - hmin).*z.^alpha_h;
handling_time = h(foraging_trait);
handling_time = reshape(handling_time,1,1,number_of_foraging);


Productivity_f = zeros(size(Param,2),101);
Productivity_t = zeros(size(Param,2),101);
for jj=1:size(Param,2)
    b_animal_f = Output_a_f(:,(jj-1)*number_of_animals*number_of_foraging+1:jj*number_of_animals*number_of_foraging);
    b_plant_f = Output_p_f(:,(jj-1)*number_of_plants+1:jj*number_of_plants);
    effort_ij_f = Effort_f(:,(jj-1)*number_of_animals*number_of_foraging*number_of_plants+1:jj*number_of_animals*number_of_foraging*number_of_plants);
    b_animal_t = Output_a_t(:,(jj-1)*number_of_animals*number_of_foraging+1:jj*number_of_animals*number_of_foraging);
    b_plant_t = Output_p_t(:,(jj-1)*number_of_plants+1:jj*number_of_plants);
    effort_ij_t = Effort_t(:,(jj-1)*number_of_animals*number_of_foraging*number_of_plants+1:jj*number_of_animals*number_of_foraging*number_of_plants);
    
    for ii=1:101
        % AF
        B_plant_f = b_plant_f(ii,:)'; % B_plant must be a vertical vector
        B_animal_f = reshape(b_animal_f(ii,:),number_of_animals,number_of_foraging);
        % reshape effort AF
        Effort_ij_f = reshape(effort_ij_f(ii,:),number_of_animals,number_of_plants,number_of_foraging);
        % compute effort mower
        Effort_sans_of_f = zeros(number_of_animals,number_of_plants);
        Effort_sans_of_f(:,B_plant_f>seuil_effort) = Effort_sans_of_f(:,B_plant_f>seuil_effort) + B_plant_f(B_plant_f>seuil_effort)';
        Effort_sans_of_f = Effort_sans_of_f./(sum(Effort_sans_of_f,2)+(sum(Effort_sans_of_f,2)==0));
        % compute real effort
        Effort_real_f(:,:,:) = Effort_ij_f(:,:,:).*Foraging_trait + (1-Foraging_trait).*Effort_sans_of_f(:,:);
        
        dc_ij_f = (1 + handling_time.*extraction_coeff.*sum((Effort_real_f.*delta_ij.*B_plant_f'),2).*dy);
        c_ij_f = (extraction_coeff.*delta_ij.*B_plant_f')./dc_ij_f;
        
        functional_response_animal_f = sum(Effort_real_f.*c_ij_f,2).*dy;
        functional_response_animal_f = reshape(functional_response_animal_f,number_of_animals,number_of_foraging);
        functional_response_animal_f = functional_response_animal_f.*B_animal_f;
        
        Productivity_f(jj,ii) = sum(functional_response_animal_f,'all')*dz*dx;
        
        % MOWER
        B_plant_t = b_plant_t(ii,:)'; % B_plant must be a vertical vector
        B_animal_t = reshape(b_animal_t(ii,:),number_of_animals,number_of_foraging);
        % reshape effort AF
        Effort_ij_t = reshape(effort_ij_t(ii,:),number_of_animals,number_of_plants,number_of_foraging);
        % compute effort mower
        Effort_sans_of_t = zeros(number_of_animals,number_of_plants);
        Effort_sans_of_t(:,B_plant_t>seuil_effort) = Effort_sans_of_t(:,B_plant_t>seuil_effort) + B_plant_t(B_plant_t>seuil_effort)';
        Effort_sans_of_t = Effort_sans_of_t./(sum(Effort_sans_of_t,2)+(sum(Effort_sans_of_t,2)==0));
        % compute real effort
        Effort_real_t(:,:,:) = Effort_ij_t(:,:,:).*Foraging_trait + (1-Foraging_trait).*Effort_sans_of_t(:,:);
        
        dc_ij_t = (1 + handling_time.*extraction_coeff.*sum((Effort_real_t.*delta_ij.*B_plant_t'),2).*dy);
        c_ij_t = (extraction_coeff.*delta_ij.*B_plant_t')./dc_ij_t;
        
        functional_response_animal_t = sum(Effort_real_t.*c_ij_t,2).*dy;
        functional_response_animal_t = reshape(functional_response_animal_t,number_of_animals,number_of_foraging);
        functional_response_animal_t = functional_response_animal_t.*B_animal_t;
        
        Productivity_t(jj,ii) = sum(functional_response_animal_t,'all')*dz*dx;
    end
end

%% BOXPLOT
% average productivity (on the last 100 time steps)
Productivity_boxplot_t = mean(Productivity_t(:,end-10:end),2);
Productivity_boxplot_f = mean(Productivity_f(:,end-10:end),2);
Productivity_boxplot_mean = [Productivity_boxplot_f,Productivity_boxplot_t];
figure(12)
clf
boxplot(Productivity_boxplot_mean);
title('Productivity','interpreter','latex','FontSize',18)
%% PRODUCTIVITY IN FUNCTION OF Z_MEAN
Z_mean = [];
for ii = 1:size(postprocess,2)
    animal = postprocess{1,ii};
    z_mean = [];
    for i = 0:10
        a = animal(end-i,:);
        a_reshape = reshape(a,number_of_animals,number_of_foraging);
        a_sum_and_norm = sum(a_reshape,1)./sum(sum(a_reshape,1));
        z_mean = [z_mean,sum(a_sum_and_norm.*foraging_trait)];
    end
    Z_mean = [Z_mean, mean(z_mean)];
end

%% PLOT DE LA DISTANCE ENTRE BIOMASSE_TONDEUSE ET BIOMASSE_AF EN FONCTION DU Z_MOYEN

figure(3)
clf
hold on
l = line([1,0],[0,0],'color','k','linewidth',2,'linestyle','-');
set(get(get(l,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%%% Absolute diffenrence
% dist_productivity_animal = Productivity_boxplot_f - Productivity_boxplot_t;
%%% Relative differnece
dist_productivity_animal = (Productivity_boxplot_f - Productivity_boxplot_t)./Productivity_boxplot_t*100;

s = scatter(Z_mean,dist_productivity_animal,5,Color(1,:),'filled');
s.MarkerFaceAlpha = .5;
hold on

% AJOUT MEAN + SHADE STD
Z_mean_round = round(Z_mean,1); % round pour réduire le nombre de Z différents
% et avoir plusieurs communautés (réplicas)
% pour chaque Z
Z_unique = unique(Z_mean_round);
% CONSUMER
dist_productivity_animal_replica = {};
for i = 1:length(Z_unique)
    value = Z_unique(i);
    dist_productivity_animal_replica{i} = dist_productivity_animal(Z_mean_round==value);
end
dist_productivity_animal_median = cellfun(@median, dist_productivity_animal_replica);
dist_productivity_animal_quantile_inf = cellfun(@(x) quantile(x,0.25), dist_productivity_animal_replica);
dist_productivity_animal_quantile_sup = cellfun(@(x) quantile(x,0.75), dist_productivity_animal_replica);
dist_productivity_animal_mean = cellfun(@mean, dist_productivity_animal_replica);
dist_productivity_animal_std = cellfun(@std, dist_productivity_animal_replica);

x = 0:.1:1;
x2 = [x, fliplr(x)];

if (choice == 1)
    p = plot(Z_unique,dist_productivity_animal_mean,'--','LineWidth',3,'color',Color(1,:));
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    inBetween = [dist_productivity_animal_mean-dist_productivity_animal_std,...
        fliplr(dist_productivity_animal_mean+dist_productivity_animal_std)];
    h = fill(x2, inBetween, Color(1,:));
    set(h,'facealpha',.2,'EdgeColor',Color(1,:));
    ymin = min(dist_productivity_animal_mean-dist_productivity_animal_std)*1.1;
    ymax = max(dist_productivity_animal_mean+dist_productivity_animal_std)*1.1;
    ylim([ymin,ymax])
elseif(choice == 2)
    p = plot(Z_unique,dist_productivity_animal_median,'--','LineWidth',3,'color',Color(1,:));
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    inBetween = [dist_productivity_animal_quantile_inf,...
        fliplr(dist_productivity_animal_quantile_sup)];
    h = fill(x2, inBetween, Color(1,:));
    set(h,'facealpha',.2,'EdgeColor',Color(1,:));
    ymin = min(dist_productivity_animal_quantile_inf)*1.1;
    ymax = max(dist_productivity_animal_quantile_sup)*1.1;
    ylim([ymin,ymax])
end
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');


hold off

ax = gca;
ax.FontSize = 16;
grid on
ax.GridLineStyle = '--';
xlabel({'mean foraging trait of','the system with AF evolution'},'interpreter','latex','FontSize',20)
ylabel({'Percent difference'; 'in productivity'},'interpreter','latex','FontSize',20)

% %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% IMAGESC OF ONE COMMUNITY
%
% % simu_animal = Out_t{1,1};
% % simu_plant = Output_p_t(1:50,1:11);
% simu_animal = Out_f{1,1};
% simu_plant = Output_p_f(1:50,1:11);
%
% subplot(1,2,1)
% % imAlpha=ones(size(simu_animal(1:50,:)));
% % imAlpha(simu_animal(1:50,:)==0)=0;
% % imagesc(heatmap,'AlphaData',imAlpha);
% cmap = jet(256);
% cmap(1,:) = 0;
% imagesc(simu_animal(1:50,:));%,'AlphaData',imAlpha);
% colormap(cmap);
% h = colorbar;
% xlabel('Niche trait $x$','interpreter','latex','FontSize',18);
% set(gca,'XTick',1:11,'XTickLabel',xx)
% ylabel(h,'Abundance','interpreter','latex','FontSize',18);
% ylabel('Time','interpreter','latex','FontSize',18)
% title('Consumers','interpreter','latex','FontSize',18)
% subplot(1,2,2)
% imagesc(simu_plant)
% xlabel('Niche trait $y$','interpreter','latex','FontSize',18);
% set(gca,'XTick',1:11,'XTickLabel',xx)
% ylabel(h,'Abundance','interpreter','latex','FontSize',18);
% ylabel('Time ','interpreter','latex','FontSize',18)
% title('Resources','interpreter','latex','FontSize',18)
% h = colorbar;
%
% % sgtitle('Mower community','interpreter','latex','FontSize',18);
% sgtitle('Foraging community','interpreter','latex','FontSize',18);
%
% % xlabel('niche trait $x$','interpreter','latex','FontSize',20);
% % set(gca,'XTick',1:11,'XTickLabel',xx)
% % caxis([0 1]);
% %
% % % lsline;
% % % colormap(myColorMap)
% % h = colorbar;
% % ylabel(h,'foraging trait $z$','interpreter','latex','FontSize',20);
%
%
% %% FIND PARAMETERS LEADING TO DYNAMIC COMMUNITY
%
% biomass_animal = cell(1,length(Out_f));
%
% for ii = 1:length(Out_f)
%     biomass_animal{ii} = sum(Out_f{ii},2);
%     biomass_eq = biomass_animal{ii};
%     biomass_eq = biomass_eq(end-50:end);
%     biomass_var(ii) = var(biomass_eq);
% end
%
% %% check the stationnary and dynamic communities
% index_statio = find(biomass_var<1);
% index_dyn = find(biomass_var>=1);
%
% for ii = 1:length(index_statio)
%     subplot(2,1,1)
%     plot(biomass_animal{index_statio(ii)});
%     hold on
% end
% title('total animal biomass for stationnary communities','interpreter','latex','FontSize',20)
% hold off
%
% for ii = 1:length(index_dyn)
%     subplot(2,1,2)
%     plot(biomass_animal{index_dyn(ii)});
%     hold on
% end
% title('total animal biomass for dynamic communities','interpreter','latex','FontSize',20)
% hold off
%
% %% find if AF evolve in the stationnary communities (not so much, AF is stronger in dynamic communities)
% close
%
% Z_final = [];
% for ii = index_statio
%     final = Out_F{ii};
%     final = reshape(final(end-10:end,:),number_of_animals,number_of_foraging);
%     final = sum(final,1);
%     z_final = final./(sum(final));
%     Z_final = [Z_final;z_final];
% end
%
% Z_final_dyn = [];
% for ii = index_dyn
%     final = Out_F{ii};
%     final = reshape(final(end-10:end,:),number_of_animals,number_of_foraging);
%     final = sum(final,1);
%     z_final = final./(sum(final));
%     Z_final_dyn = [Z_final_dyn;z_final];
% end
%
% plot(foraging_trait,sum(Z_final,1),'-*','LineWidth',2)
% hold on
% plot(foraging_trait,sum(Z_final_dyn,1),'r-*','LineWidth',2)
% hold off
% legend('stationnary communities','dynamic communities','interpreter','latex','FontSize',15)
% xlabel('foraging trait','interpreter','latex','FontSize',20)
% ylabel('mean densities over all runs','interpreter','latex','FontSize',20)
%
% set(gcf, 'Position', get(0, 'Screensize'));
%
% %% parameters leading to stationnary or dynamic communities
% % Param = [SIGMA,sigmaK,hmax,alpha_h]];
% close
%
% subplot(1,4,1)
% plot(Param(1,index_statio),'b*')
% hold on
% plot(Param(1,index_dyn),'r*')
% title('$\sigma$','interpreter','latex','FontSize',15)
%
% subplot(1,4,2)
% plot(Param(2,index_statio),'b*')
% hold on
% plot(Param(2,index_dyn),'r*')
% title('$\sigma_K$','interpreter','latex','FontSize',15)
%
% subplot(1,4,3)
% plot(Param(3,index_statio),'b*')
% hold on
% plot(Param(3,index_dyn),'r*')
% title('$h_{max}$','interpreter','latex','FontSize',15)
%
% subplot(1,4,4)
% plot(Param(4,index_statio),'b*')
% hold on
% plot(Param(4,index_dyn),'r*')
% title('$\alpha_h$','interpreter','latex','FontSize',15)
% legend('stationnary communities','dynamical communities','interpreter','latex','Location','northwest','FontSize',15)
%
% set(gcf, 'Position', get(0, 'Screensize'));
%
% %% Test of visualization of community-proxy
%
% % dominant trait
% close
%
% simu = Out_F{5};
%
% dominant_x = cell(1,101);
% for ii = 1:101
%     matrix = simu((ii-1)*number_of_animals+1:ii*number_of_animals,:);
%     sum_on_z = sum(matrix,2);
%     dominant_x{ii} = find(sum_on_z==max(sum_on_z));
% end
%
% for ii = 1:length(dominant_x)
%     plot(ii,xx(dominant_x{ii}),'b-*')
%     hold on
% end
%
%
%
%
%
%
%
%
