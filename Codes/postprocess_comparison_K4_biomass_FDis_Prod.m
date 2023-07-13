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

% load('../Data/community_comparison_K4_0.mat'); % community_comparison.mat + community_comparison_b.mat
data = [1:23,31:66]; 
Ndata = length(data);
Iloop = 20;
Nt = 10;
loop = Ndata*Iloop;
Output_a_f = zeros(loop,Nt+1,number_of_animals,number_of_foraging);
Output_p_f = zeros(loop,Nt+1,number_of_plants);
Output_a_t = zeros(loop,Nt+1,number_of_animals,number_of_foraging);
Output_p_t = zeros(loop,Nt+1,number_of_plants);
Effort_f = zeros(loop,Nt+1,number_of_animals,number_of_plants,number_of_foraging);
Effort_t = zeros(loop,Nt+1,number_of_animals,number_of_plants,number_of_foraging);
PARAM    = zeros(Iloop*Ndata,4);
for idata = 1:Ndata
    load(['../Data/community_comparison_K4_',num2str(data(idata)),'.mat']);
    Output_a_f((idata-1)*Iloop+1:idata*Iloop,:,:,:) = output_a_f; 
    Output_a_t((idata-1)*Iloop+1:idata*Iloop,:,:,:) = output_a_t ;
    Output_p_f((idata-1)*Iloop+1:idata*Iloop,:,:) = output_p_f; 
    Output_p_t((idata-1)*Iloop+1:idata*Iloop,:,:) = output_p_t; 
    Effort_f((idata-1)*Iloop+1:idata*Iloop,:,:,:,:) = effort_f;
    Effort_t((idata-1)*Iloop+1:idata*Iloop,:,:,:,:) = effort_t;
    PARAM((idata-1)*Iloop+1:idata*Iloop,:) = Param(end-Iloop+1:end,:);
end

%% Biomass
Biomass_plant_f = mean(sum(Output_p_f,3),2);
Biomass_plant_t = mean(sum(Output_p_t,3),2);
Biomass_plant_boxplot = [Biomass_plant_f,Biomass_plant_t];

Biomass_animal_f = mean(sum(Output_a_f,[3,4]),2);
Biomass_animal_t = mean(sum(Output_a_t,[3,4]),2);
Biomass_animal_boxplot = [Biomass_animal_f,Biomass_animal_t];
%%% Percentage of differnece rate
dist_biomasse_animal = (Biomass_animal_boxplot(:,1) - Biomass_animal_boxplot(:,2))./Biomass_animal_boxplot(:,2)*100;
dist_biomasse_resource = (Biomass_plant_boxplot(:,1) - Biomass_plant_boxplot(:,2))./Biomass_plant_boxplot(:,2)*100;

Z_mean = mean(sum(sum(Output_a_f,3)./sum(Output_a_f,[3,4]).*reshape(foraging_trait,[1,1,1,number_of_foraging]),4),2);

%% FIGURE BIOMASS
figure(1)
clf
hold on
l = line([1,0],[0,0],'color','k','linewidth',2,'linestyle','-');
set(get(get(l,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% yyaxis left
s = scatter(Z_mean,dist_biomasse_animal,5,Color(1,:),'filled');
s.MarkerFaceAlpha = .5;

% yyaxis right
s = scatter(Z_mean,dist_biomasse_resource,5,'d','filled','MarkerFaceColor',Color(5,:));
s.MarkerFaceAlpha = .5;
axis([0,1,-50,100])
legend('consumers','resources','interpreter','latex','FontSize',15,'location','northwest')


% AJOUT MEAN + SHADE STD
Z_mean_round = round(Z_mean,1); % round pour réduire le nombre de Z différents
% et avoir plusieurs communautés (réplicas)
% pour chaque 
x =  0:.1:1;
Z_unique = x';
% Z_unique = unique(Z_mean_round);
% RESOURCE / CONSUMER
dist_biomasse_resource_replica = {};
dist_biomasse_animal_replica = {};

for i = 1:length(Z_unique)
    value = Z_unique(i);
    Iz = logical((Z_mean_round<value+0.01).*(Z_mean_round>value-0.01));
    dist_biomasse_resource_replica{i} = dist_biomasse_resource(Iz);
    dist_biomasse_animal_replica{i}   = dist_biomasse_animal(Iz);
end
% RESOURCE
dist_biomasse_resource_median = cellfun(@median, dist_biomasse_resource_replica);
dist_biomasse_resource_quantile_inf = cellfun(@(x) quantile(x,0.25), dist_biomasse_resource_replica);
dist_biomasse_resource_quantile_sup = cellfun(@(x) quantile(x,0.75), dist_biomasse_resource_replica);
dist_biomasse_resource_mean = cellfun(@mean, dist_biomasse_resource_replica);
dist_biomasse_resource_std = cellfun(@std, dist_biomasse_resource_replica);
% CONSUMER
dist_biomasse_animal_median = cellfun(@median, dist_biomasse_animal_replica);
dist_biomasse_animal_quantile_inf = cellfun(@(x) quantile(x,0.25), dist_biomasse_animal_replica);
dist_biomasse_animal_quantile_sup = cellfun(@(x) quantile(x,0.75), dist_biomasse_animal_replica);
dist_biomasse_animal_mean = cellfun(@mean, dist_biomasse_animal_replica);
dist_biomasse_animal_std = cellfun(@std, dist_biomasse_animal_replica);

x2 = [x, fliplr(x)];

dist_biomasse_resource_mean = movmean(dist_biomasse_resource_mean,2);
dist_biomasse_resource_std  = movmean(dist_biomasse_resource_std,2);
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
% yticks(ax,[-75 -50 -25 0 25 50 75])
% yticklabels(ax,{'-75','-50','-25','0','25','50','75'})
grid on
ax.GridLineStyle = '--';
xlabel({'mean foraging trait of','the system with AF evolution'},'interpreter','latex','FontSize',20)
ylabel({'Percent difference'; 'in biomass'},'interpreter','latex','FontSize',20,'color','k')


%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONAL DIVERSITY %
%% Laliberté & Legendre 2010 %
%%%%%%%%%%%%%%%%%%%%%%%%

FDis_animal_f = zeros(loop,Nt+1);
FDis_animal_t = zeros(loop,Nt+1);
FDis_plant_f = zeros(loop,Nt+1);
FDis_plant_t = zeros(loop,Nt+1);
% Densities
animal_density_f = Output_a_f;
plant_density_f = Output_p_f;
animal_density_t = Output_a_t;
plant_density_t = Output_p_t;

% Centroid
C_animal_f = sum(animal_density_f.*reshape(xx,[1,1,number_of_animals,1]),[3,4])./sum(animal_density_f,[3,4]);
C_animal_t = sum(animal_density_t.*reshape(xx,[1,1,number_of_animals,1]),[3,4])./sum(animal_density_t,[3,4]);
C_plant_f  = sum(plant_density_f.*reshape(xx,[1,1,number_of_plants]),3)./sum(plant_density_f,3);
C_plant_t  = sum(plant_density_t.*reshape(xx,[1,1,number_of_plants]),3)./sum(plant_density_t,3);


% Functional dispersion
z_animal_f = abs(xx-C_animal_f);
FDis_animal_f = sum(animal_density_f.*reshape(z_animal_f,[loop,1,number_of_animals,1]),[3,4])./sum(animal_density_f,[3,4]);
z_animal_t    = abs(xx-C_animal_t);
FDis_animal_t = sum(animal_density_t.*reshape(z_animal_t,[loop,1,number_of_animals,1]),[3,4])./sum(animal_density_t,[3,4]);
z_plant_f     = abs(xx-C_plant_f);
FDis_plant_f  = sum(plant_density_f.*reshape(z_plant_f,[loop,1,number_of_plants]),3)./sum(plant_density_f,3);
z_plant_t     = abs(xx-C_plant_t);
FDis_plant_t  = sum(plant_density_t.*reshape(z_plant_t,[loop,1,number_of_plants]),3)./sum(plant_density_t,3);

FDis_animal_boxplot = [mean(FDis_animal_f,2)';mean(FDis_animal_t,2)'];
FDis_plant_boxplot = [mean(FDis_plant_f,2)';mean(FDis_plant_t,2)'];


%% FIGURE FDIS
figure(2)
clf
hold on
l = line([1,0],[0,0],'color','k','linewidth',2,'linestyle','-');
set(get(get(l,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

%%% Relative difference
dist_FDis_animal = (FDis_animal_boxplot(1,:) - FDis_animal_boxplot(2,:))./FDis_animal_boxplot(2,:)*100;

s = scatter(Z_mean,dist_FDis_animal,5,Color(1,:),'filled');
s.MarkerFaceAlpha = .5;
hold on

%%% Relative difference
dist_FDis_resource = (FDis_plant_boxplot(1,:) - FDis_plant_boxplot(2,:) )./FDis_plant_boxplot(2,:)*100;

s = scatter(Z_mean,dist_FDis_resource,5,'d','filled','MarkerFaceColor',Color(5,:));
s.MarkerFaceAlpha = .5;

legend('consumers','resources','interpreter','latex','FontSize',15,'location','northwest')
xlabel({'mean foraging trait','of community with AF evolution'},'interpreter','latex','FontSize',20)

% AJOUT MEAN + SHADE STD
% Z_mean_round = round(Z_mean,1); % round pour réduire le nombre de Z différents
% et avoir plusieurs communautés (réplicas)
% pour chaque Z
% x =  0:.1:1;
% Z_unique = x';
% Z_unique = unique(Z_mean_round);
% RESOURCE / CONSUMER
dist_FDis_resource_replica = {};
dist_FDis_animal_replica = {};

for i = 1:length(Z_unique)
    value = Z_unique(i);
    Iz = logical((Z_mean_round<value+0.01).*(Z_mean_round>value-0.01));
    dist_FDis_resource_replica{i} = dist_FDis_resource(Iz);
    dist_FDis_animal_replica{i}   = dist_FDis_animal(Iz);
end
% RESOURCE
dist_FDis_resource_median = cellfun(@median, dist_FDis_resource_replica);
dist_FDis_resource_quantile_inf = cellfun(@(x) quantile(x,0.25), dist_FDis_resource_replica);
dist_FDis_resource_quantile_sup = cellfun(@(x) quantile(x,0.75), dist_FDis_resource_replica);
dist_FDis_resource_mean = cellfun(@mean, dist_FDis_resource_replica);
dist_FDis_resource_std = cellfun(@std, dist_FDis_resource_replica);
% CONSUMER
dist_FDis_animal_median = cellfun(@nanmedian, dist_FDis_animal_replica);
dist_FDis_animal_quantile_inf = cellfun(@(x) quantile(x,0.25), dist_FDis_animal_replica);
dist_FDis_animal_quantile_sup = cellfun(@(x) quantile(x,0.75), dist_FDis_animal_replica);
dist_FDis_animal_mean = cellfun(@mean, dist_FDis_animal_replica);
dist_FDis_animal_std = cellfun(@std, dist_FDis_animal_replica);


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PRODUCTIVITY %%%
% Sum of all consumers functional responses weigthed by the consumers abundances
% = flux from resources to consumers

% load('../Data/community_comparison_all.mat'); % community_comparison.mat + community_comparison_b.mat
% load('../Data/community_comparison_K4.mat'); % community_comparison.mat + community_comparison_b.mat

number_of_animals = 11;
number_of_plants  = 11;
number_of_foraging = 11;

% FONCTION K(y)
K0 = 50; % maximal carrying capacity
y0 = 0; % optimal trait
sigmaK = 2.5;
% Kf = @(trait_plant) K0*exp(-(trait_plant - y0).^2./(2*sigmaK^2));
Kf = @(trait_plant) K0*exp(-(trait_plant - y0).^4./(12*sigmaK^4));


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


Productivity_f = zeros(loop,Nt+1);
Productivity_t = zeros(loop,Nt+1);
for jj=1:loop 
    sigma_loop  = PARAM(:,1);
    sigmaK_loop = PARAM(:,2);
    hmax_loop   = PARAM(:,3);
    animal_intrinsic_growth_loop = PARAM(:,4);
    %% SIGMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SIGMA = sigma_loop(jj);
    delta_ij = complementary_traits(SIGMA,XX,YY);
    delta_ij = delta_ij';
    
    %% SIGMA_K %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sigmaK = sigmaK_loop(jj);
%   Kf = @(trait_plant) K0*exp(-(trait_plant - y0).^2./(2*sigmaK^2));
    Kf = @(trait_plant) K0*exp(-(trait_plant - y0).^4./(12*sigmaK^4));
    
    sigmaC = sigmaK-1;
    C = @(trait_plant_y,trait_plant_yi) exp(sigmaC^2*Beta^2/2).*exp(-(trait_plant_y-trait_plant_yi + sigmaC^2*Beta^2).^2./(2*sigmaC^2));
    Compet = C(traits_of_plants,traits_of_plants');
    
    %% h_max & alpha_h %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hmax = hmax_loop(jj);
    h = @(z) hmin + (hmax - hmin).*z.^alpha_h;
    handling_time = h(foraging_trait);
    handling_time = reshape(handling_time,1,1,number_of_foraging);
    
    % animal_growth
    animal_intrinsic_growth = animal_intrinsic_growth_loop(jj);
    
    for ii=1:Nt+1
        % AF
        B_plant_f = permute(Output_p_f(jj,ii,:),[3,1,2]); % B_plant must be a vertical vector
        B_animal_f = permute(Output_a_f(jj,ii,:,:),[3,4,1,2]);
        % reshape effort AF
        Effort_ij_f = permute(Effort_f(jj,ii,:,:,:),[3,4,5,1,2]);
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
        B_plant_t  = permute(Output_p_t(jj,ii,:),[3,1,2]); % B_plant must be a vertical vector
        B_animal_t = permute(Output_a_t(jj,ii,:,:),[3,4,1,2]);
        % reshape effort AF
        Effort_ij_t = permute(Effort_t(jj,ii,:,:,:),[3,4,5,1,2]); % compute effort mower
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

%% FIGURE PRODUCTIVITY
Productivity_boxplot_t = mean(Productivity_t,2);
Productivity_boxplot_f = mean(Productivity_f,2);
Productivity_boxplot_mean = [Productivity_boxplot_f,Productivity_boxplot_t];

figure(3)
clf
hold on
l = line([1,0],[0,0],'color','k','linewidth',2,'linestyle','-');
set(get(get(l,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%%% Absolute differnece
% dist_productivity_animal = (Productivity_boxplot_f - Productivity_boxplot_t);

%%% Relative differnece
dist_productivity_animal = (Productivity_boxplot_f - Productivity_boxplot_t)./Productivity_boxplot_t*100;

s = scatter(Z_mean,dist_productivity_animal,5,Color(1,:),'filled');
s.MarkerFaceAlpha = .5;
hold on

% AJOUT MEAN + SHADE STD
dist_productivity_animal_replica = {};
% CONSUMER
for i = 1:length(Z_unique)
    value = Z_unique(i);
    Iz = logical((Z_mean_round<value+0.01).*(Z_mean_round>value-0.01));
    dist_productivity_animal_replica{i} = dist_productivity_animal(Iz);
end

dist_productivity_animal_median = cellfun(@nanmedian, dist_productivity_animal_replica);
dist_productivity_animal_quantile_inf = cellfun(@(x) quantile(x,0.25), dist_productivity_animal_replica);
dist_productivity_animal_quantile_sup = cellfun(@(x) quantile(x,0.75), dist_productivity_animal_replica);
dist_productivity_animal_mean = cellfun(@mean, dist_productivity_animal_replica);
dist_productivity_animal_std = cellfun(@std, dist_productivity_animal_replica);

% x = 0:.1:1;
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

