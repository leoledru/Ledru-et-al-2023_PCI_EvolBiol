clear all;clc

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

load('../community_comparison_all.mat');

% each column is a simulation, row 1 = consumers (animal) biomass AF community; row 2 = resources (plant) biomass AF community;
% row 3 = foraging efforts AF community; row 4 = consumers biomass no-AF community; row 5 = resources biomass no-AF community
% row 6 = foraging efforts no-AF community; row 7 = parameters for each run

for ii = 1:size(Param,2)
  postprocess{1,ii} = Output_a_f(:,(ii-1)*number_of_animals*number_of_foraging+1:ii*number_of_animals*number_of_foraging);
  postprocess{2,ii} = Output_p_f(:,(ii-1)*number_of_plants+1:ii*number_of_plants);
%   postprocess{3,ii} = Effort_f(:,(ii-1)*number_of_animals*number_of_foraging+1:ii*number_of_animals*number_of_foraging*number_of_plants);
  postprocess{4,ii} = Output_a_t(:,(ii-1)*number_of_animals*number_of_foraging+1:ii*number_of_animals*number_of_foraging);
  postprocess{5,ii} = Output_p_t(:,(ii-1)*number_of_plants+1:ii*number_of_plants);
%   postprocess{6,ii} = Effort_t(:,(ii-1)*number_of_animals*number_of_foraging+1:ii*number_of_animals*number_of_foraging*number_of_plants);
  postprocess{7,ii} = Param(:,ii); 
end

%% Compute this section for agregation of densities over z axis (necessary
%% to further compute the diversity indices which only consider the diversity along x axis).

% Reshape animal density (postprocess rows 1 and 4)
Out_f = cell(1,size(postprocess,2));
Out_t = cell(1,size(postprocess,2));
Out_F = cell(1,size(postprocess,2));
Out_T = cell(1,size(postprocess,2));
for ii = 1:size(postprocess,2)
    Output_a_f = [];
    Output_a_t = [];
    output_a_f = postprocess{1,ii};
    output_a_t = postprocess{4,ii};
    for jj = 1:size(output_a_f,1)
        Output_a_f = [Output_a_f;reshape(output_a_f(jj,:),number_of_animals,number_of_foraging)];
        Output_a_t = [Output_a_t;reshape(output_a_t(jj,:),number_of_animals,number_of_foraging)];    
    end
    Out_F{1,ii} = Output_a_f;
    Out_T{1,ii} = Output_a_t;
end

% Sum animal densities on foraging trait
animal_niche_f = [];
for ii = 1:size(postprocess,2)
    animal_niche_f = [];
    animal_f = Out_F{:,ii};
    animal_niche_t = [];
    animal_t = Out_T{:,ii};
    for jj = 1:101
    animal_niche_f = [animal_niche_f;sum(animal_f((jj-1)*number_of_foraging+1:jj*number_of_foraging,:),2)'];
    animal_niche_t = [animal_niche_t;sum(animal_t((jj-1)*number_of_foraging+1:jj*number_of_foraging,:),2)'];
    end
    Out_f{:,ii} = animal_niche_f;
    Out_t{:,ii} = animal_niche_t;
end

%% Plot of distance between biomass of AF community and biomass of no-AF community, 
%% in function of the mean foraging trait of AF community

biomass_animal_f = [];
biomass_animal_t = [];
biomass_plant_f = [];
biomass_plant_t = [];
for ii = 1:size(Param,2)
    biomass_animal_f = [biomass_animal_f;sum(postprocess{1,ii},2)'];
    biomass_animal_t = [biomass_animal_t;sum(postprocess{4,ii},2)'];
    biomass_plant_f = [biomass_plant_f;sum(postprocess{2,ii},2)'];
    biomass_plant_t = [biomass_plant_t;sum(postprocess{5,ii},2)'];
end

Biomass_plant_f = mean(biomass_plant_f(:,end-10:end),2);
Biomass_plant_t = mean(biomass_plant_t(:,end-10:end),2);
Biomass_plant_boxplot = [Biomass_plant_f,Biomass_plant_t];

Biomass_animal_f = mean(biomass_animal_f(:,end-10:end),2);
Biomass_animal_t = mean(biomass_animal_t(:,end-10:end),2);
Biomass_animal_boxplot = [Biomass_animal_f,Biomass_animal_t];

% Compute z average for color in scatters below
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

dist_biomasse_animal = Biomass_animal_boxplot(:,1) - Biomass_animal_boxplot(:,2); 
yyaxis left
s = scatter(Z_mean,dist_biomasse_animal,5,'b','filled');
s.MarkerFaceAlpha = .5;

dist_biomasse_resource = Biomass_plant_boxplot(:,1) - Biomass_plant_boxplot(:,2); 
yyaxis right
s = scatter(Z_mean,dist_biomasse_resource,5,'d','filled','MarkerFaceColor',[.2 .8 .2]);
s.MarkerFaceAlpha = .5;

legend('consumers','resources','interpreter','latex','FontSize',15,'location','southeast')
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = [.2 .8 .2];

% add mean + std
Z_mean_round = round(Z_mean,1); % round to reduce the number of different mean foraging traits and have several communities (replicas) for each
Z_unique = unique(Z_mean_round); 
% resources
dist_biomasse_resource_replica = {};
for i = 1:length(Z_unique)
    value = Z_unique(i);
    dist_biomasse_resource_replica{i} = dist_biomasse_resource(Z_mean_round==value);
end
dist_biomasse_resource_median = cellfun(@median, dist_biomasse_resource_replica);
dist_biomasse_resource_mean = cellfun(@mean, dist_biomasse_resource_replica);
dist_biomasse_resource_std = cellfun(@std, dist_biomasse_resource_replica);
% consumers
dist_biomasse_animal_replica = {};
for i = 1:length(Z_unique)
    value = Z_unique(i);
    dist_biomasse_animal_replica{i} = dist_biomasse_animal(Z_mean_round==value);
end
dist_biomasse_animal_median = cellfun(@median, dist_biomasse_animal_replica);
dist_biomasse_animal_mean = cellfun(@mean, dist_biomasse_animal_replica);
dist_biomasse_animal_std = cellfun(@std, dist_biomasse_animal_replica);

x = 0:.1:1;
x2 = [x, fliplr(x)];

dist_biomasse_resource_mean = movmean(dist_biomasse_resource_mean,2);
dist_biomasse_resource_std = movmean(dist_biomasse_resource_std,2);
yyaxis right
hold on
p = plot(Z_unique,dist_biomasse_resource_mean,'--','LineWidth',3,'Color',[.2 .8 .2]);
set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
inBetween = [dist_biomasse_resource_mean-dist_biomasse_resource_std,...
    fliplr(dist_biomasse_resource_mean+dist_biomasse_resource_std)];
h = fill(x2, inBetween, [.2 .8 .2]);
set(h,'facealpha',.2,'EdgeColor',[.2 .8 .2]);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hold off

dist_biomasse_animal_mean = movmean(dist_biomasse_animal_mean,2);
dist_biomasse_animal_std = movmean(dist_biomasse_animal_std,2);
yyaxis left
hold on
p = plot(Z_unique,dist_biomasse_animal_mean,'--b','LineWidth',3);
set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
inBetween = [dist_biomasse_animal_mean-dist_biomasse_animal_std,...
    fliplr(dist_biomasse_animal_mean+dist_biomasse_animal_std)];
h = fill(x2, inBetween, 'b');
set(h,'facealpha',.2,'EdgeColor','b');
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hold off

ax = gca;
ax.FontSize = 16;
xlabel({'mean foraging trait of','community with AF evolution'},'interpreter','latex','FontSize',20)
ylabel({'difference in','resource biomass'},'interpreter','latex','FontSize',20)
yyaxis right
ylabel({'difference in','resource biomass'},'interpreter','latex','FontSize',20)

%%
x = Z_mean;
y = dist_biomasse_animal;
mdl = fitlm(x',y');
Xnew = linspace(min(x), max(x), 1000)';
[ypred,yci] = predict(mdl, Xnew);
figure
plot(x, y, 'p')
hold on
plot(Xnew, ypred, '--g')
plot(Xnew, yci, '--r')
hold off
grid

%% FUNCTIONAL DIVERSITY %
% Laliberté & Legendre 2010 %

FDis_animal_f = [];
FDis_animal_t = [];
FDis_plant_f = [];
FDis_plant_t = [];
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
    FDis_animal_f = [FDis_animal_f,sum(animal_density_f.*z_animal_f,2)./sum(animal_density_f,2)]; 
    z_animal_t = abs(xx-C_animal_t);
    FDis_animal_t = [FDis_animal_t,sum(animal_density_t.*z_animal_t,2)./sum(animal_density_t,2)];
    z_plant_f = abs(xx-C_plant_f);
    FDis_plant_f = [FDis_plant_f,sum(plant_density_f.*z_plant_f,2)./sum(plant_density_f,2)]; 
    z_plant_t = abs(xx-C_plant_t);
    FDis_plant_t = [FDis_plant_t,sum(plant_density_t.*z_plant_t,2)./sum(plant_density_t,2)];   
end

% average on the last 100 time steps
FDis_animal_boxplot = [mean(FDis_animal_f(end-10:end,:),1);mean(FDis_animal_t(end-10:end,:),1)];
FDis_plant_boxplot = [mean(FDis_plant_f(end-10:end,:),1);mean(FDis_plant_t(end-10:end,:),1)];

% Compute z average for color in scatters below
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

%% Plot of distance between FDis of AF community and FDis of no-AF community, 
%% in function of the mean foraging trait of AF community
dist_FDis_animal = FDis_animal_boxplot(1,:) - FDis_animal_boxplot(2,:); 
s = scatter(Z_mean,dist_FDis_animal,5,'b','filled');
s.MarkerFaceAlpha = .5;
hold on

dist_FDis_resource = FDis_plant_boxplot(1,:) - FDis_plant_boxplot(2,:); 
s = scatter(Z_mean,dist_FDis_resource,5,'d','filled','MarkerFaceColor',[.2 .8 .2]);
s.MarkerFaceAlpha = .5;

legend('consumers','resources','interpreter','latex','FontSize',15,'location','southeast')
xlabel({'mean foraging trait','of community with AF evolution'},'interpreter','latex','FontSize',20)

% add mean + std
Z_mean_round = round(Z_mean,1);
Z_unique = unique(Z_mean_round); 
% resources
dist_FDis_resource_replica = {};
for i = 1:length(Z_unique)
    value = Z_unique(i);
    dist_FDis_resource_replica{i} = dist_FDis_resource(Z_mean_round==value);
end
dist_FDis_resource_median = cellfun(@median, dist_FDis_resource_replica);
dist_FDis_resource_mean = cellfun(@mean, dist_FDis_resource_replica);
dist_FDis_resource_std = cellfun(@std, dist_FDis_resource_replica);
% consumers
dist_FDis_animal_replica = {};
for i = 1:length(Z_unique)
    value = Z_unique(i);
    dist_FDis_animal_replica{i} = dist_FDis_animal(Z_mean_round==value);
end
dist_FDis_animal_median = cellfun(@median, dist_FDis_animal_replica);
dist_FDis_animal_mean = cellfun(@mean, dist_FDis_animal_replica);
dist_FDis_animal_std = cellfun(@std, dist_FDis_animal_replica);

x = 0:.1:1;
x2 = [x, fliplr(x)];

dist_FDis_resource_mean = movmean(dist_FDis_resource_mean,2);
dist_FDis_resource_std = movmean(dist_FDis_resource_std,2);
hold on
p = plot(Z_unique,dist_FDis_resource_mean,'--','LineWidth',3,'Color',[.2 .8 .2]);
set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
inBetween = [dist_FDis_resource_mean-dist_FDis_resource_std,...
    fliplr(dist_FDis_resource_mean+dist_FDis_resource_std)];
h = fill(x2, inBetween, [.2 .8 .2]);
set(h,'facealpha',.2,'EdgeColor',[.2 .8 .2]);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

dist_FDis_animal_mean = movmean(dist_FDis_animal_mean,2);
dist_FDis_animal_std = movmean(dist_FDis_animal_std,2);
p = plot(Z_unique,dist_FDis_animal_mean,'--b','LineWidth',3);
set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
inBetween = [dist_FDis_animal_mean-dist_FDis_animal_std,...
    fliplr(dist_FDis_animal_mean+dist_FDis_animal_std)];
h = fill(x2, inBetween, 'b');
set(h,'facealpha',.2,'EdgeColor','b');
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hold off

ax = gca;
ax.FontSize = 16;
xlabel({'mean foraging trait of','community with AF evolution'},'interpreter','latex','FontSize',20)
ylabel({'difference in','functional diversity FDis'},'interpreter','latex','FontSize',20)

%% FRO : evenness %

animal_density_f_end = [];
animal_density_t_end = [];
plant_density_f_end = [];
plant_density_t_end = [];

for ii = 1:size(postprocess,2)  
%     Densities
    animal_density_f = Out_f{:,ii};
    plant_density_f = postprocess{2,ii};
    animal_density_t = Out_t{:,ii};
    plant_density_t = postprocess{5,ii};
    
    % average on the last 100 time step
    animal_density_f_end = [animal_density_f_end;mean(animal_density_f(end-10:end,:),1)];
    animal_density_t_end = [animal_density_t_end;mean(animal_density_t(end-10:end,:),1)];
    plant_density_f_end = [plant_density_f_end;mean(plant_density_f(end-10:end,:),1)];
    plant_density_t_end = [plant_density_t_end;mean(plant_density_t(end-10:end,:),1)];
end

EW_animal_f = [];
EW_animal_t = [];
EW_plant_f = [];
EW_plant_t = [];
for jj = 1:(number_of_animals-1)
    EW_animal_f = [EW_animal_f,abs(xx(jj+1) - xx(jj))./(animal_density_f_end(:,jj+1)+animal_density_f_end(:,jj))];
    EW_animal_t = [EW_animal_t,abs(xx(jj+1) - xx(jj))./(animal_density_t_end(:,jj+1)+animal_density_t_end(:,jj))];
    EW_plant_f = [EW_plant_f,abs(xx(jj+1) - xx(jj))./(plant_density_f_end(:,jj+1)+plant_density_f_end(:,jj))];
    EW_plant_t = [EW_plant_t,abs(xx(jj+1) - xx(jj))./(plant_density_t_end(:,jj+1)+plant_density_t_end(:,jj))];
end

PEW_animal_f = EW_animal_f./sum(EW_animal_f,2);
PEW_animal_t = EW_animal_t./sum(EW_animal_t,2);
PEW_plant_f = EW_plant_f./sum(EW_plant_f,2);
PEW_plant_t = EW_plant_t./sum(EW_plant_t,2);

FRO_animal_f = [];
FRO_animal_t = [];
FRO_plant_f = [];
FRO_plant_t = [];
for k = 1:size(postprocess,2) 
FRO_animal_f = [FRO_animal_f;sum(min(PEW_animal_f(k,:),1/(number_of_animals-1)))];
FRO_animal_t = [FRO_animal_t;sum(min(PEW_animal_t(k,:),1/(number_of_animals-1)))];
FRO_plant_f = [FRO_plant_f;sum(min(PEW_plant_f(k,:),1/(number_of_animals-1)))];
FRO_plant_t = [FRO_plant_t;sum(min(PEW_plant_t(k,:),1/(number_of_animals-1)))];
end

FRO_boxplot_animal = [FRO_animal_f,FRO_animal_t];
FRO_boxplot_plant = [FRO_plant_f,FRO_plant_t];

% Compute z average for color in scatters below
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

%% Plot of distance between FRO of AF community and FRO of no-AF community, 
%% in function of the mean foraging trait of AF community
dist_FRO_animal = FRO_boxplot_animal(:,1) - FRO_boxplot_animal(:,2); 
s = scatter(Z_mean,dist_FRO_animal,5,'b','filled');
s.MarkerFaceAlpha = .5;
hold on
dist_FRO_resource = FRO_boxplot_plant(:,1) - FRO_boxplot_plant(:,2); 
s = scatter(Z_mean,dist_FRO_resource,5,'d','filled','MarkerFaceColor',[.2 .8 .2]);
s.MarkerFaceAlpha = .5;

legend('consumers','resources','interpreter','latex','FontSize',15,'location','southeast')

% add mean + std
Z_mean_round = round(Z_mean,1);
Z_unique = unique(Z_mean_round); 
% resources
dist_FRO_resource_replica = {};
for i = 1:length(Z_unique)
    value = Z_unique(i);
    dist_FRO_resource_replica{i} = dist_FRO_resource(Z_mean_round==value);
end
dist_FRO_resource_median = cellfun(@median, dist_FRO_resource_replica);
dist_FRO_resource_mean = cellfun(@mean, dist_FRO_resource_replica);
dist_FRO_resource_std = cellfun(@std, dist_FRO_resource_replica);
% consumers
dist_FRO_animal_replica = {};
for i = 1:length(Z_unique)
    value = Z_unique(i);
    dist_FRO_animal_replica{i} = dist_FRO_animal(Z_mean_round==value);
end
dist_FRO_animal_median = cellfun(@median, dist_FRO_animal_replica);
dist_FRO_animal_mean = cellfun(@mean, dist_FRO_animal_replica);
dist_FRO_animal_std = cellfun(@std, dist_FRO_animal_replica);

x = 0:.1:1;
x2 = [x, fliplr(x)];

dist_FRO_resource_mean = movmean(dist_FRO_resource_mean,2);
dist_FRO_resource_std = movmean(dist_FRO_resource_std,2);
hold on
p = plot(Z_unique,dist_FRO_resource_mean,'--','LineWidth',3,'Color',[.2 .8 .2]);
set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
inBetween = [dist_FRO_resource_mean-dist_FRO_resource_std,...
    fliplr(dist_FRO_resource_mean+dist_FRO_resource_std)];
h = fill(x2, inBetween, [.2 .8 .2]);
set(h,'facealpha',.2,'EdgeColor',[.2 .8 .2]);
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

dist_FRO_animal_mean = movmean(dist_FRO_animal_mean,2);
dist_FRO_animal_std = movmean(dist_FRO_animal_std,2);
p = plot(Z_unique,dist_FRO_animal_mean,'--b','LineWidth',3);
set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
inBetween = [dist_FRO_animal_mean-dist_FRO_animal_std,...
    fliplr(dist_FRO_animal_mean+dist_FRO_animal_std)];
h = fill(x2, inBetween, 'b');
set(h,'facealpha',.2,'EdgeColor','b');
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hold off

ax = gca;
ax.FontSize = 16;
xlabel({'mean foraging trait of','community with AF evolution'},'interpreter','latex','FontSize',20)
ylabel({'difference in','functional regularity FRO'},'interpreter','latex','FontSize',20)

%% PRODUCTIVITY + Regime-morpho relationship
% Sum of all consumers functional responses weigthed by the consumers abundances = flux from resources to consumers

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
niche_plants = 1:number_of_plants; % redéfini la position de niche comme allant de 1 à nb pour le calcul du régime
niche_consumers = 1:number_of_animals;

Productivity_f = [];
Productivity_t = [];
Regime_morpho_cor_t = [];
Regime_morpho_cor_f = [];
for jj=1:size(Param,2)
    productivity_f = [];
    productivity_t = []; 
    regime_morpho_cor_f = [];
    regime_morpho_cor_t = [];
    
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
        
        % calcul du matching entre régime réalisé et morphologie (x)
        animal_vector = reshape(B_animal_f, number_of_animals, 1, number_of_foraging);
        animal_vector_norm = animal_vector./sum(animal_vector, "all");
        regime_f = Effort_real_f.*c_ij_f.*animal_vector_norm; % régime -> multiplié par le vecteur normalisé de biomasse pour ne garder que les régimes *réalisés*
        regime_norm = regime_f./sum(regime_f,1); 
        regime_norm(isnan(regime_norm)) = 0;
        regime_y_mean = sum(regime_norm.*niche_plants,2); % régime rapporté en fonction de la position de niche de la ressource
        matching_x_y = abs(niche_consumers' - regime_y_mean); % distance entre position de niche du consommateur et de la "position" de son régime
        mean_matching_x_y = nanmean(matching_x_y, "all"); % mean distance pour toutes la communauté

        productivity_f = [productivity_f,sum(functional_response_animal_f,'all')*dz*dx];
        regime_morpho_cor_f = [regime_morpho_cor_f, mean_matching_x_y];
        
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
        
        % calcul du matching entre régime réalisé et morphologie (x)
        animal_vector = reshape(B_animal_t, number_of_animals, 1, number_of_foraging);
        animal_vector_norm = animal_vector./sum(animal_vector, "all");
        regime_t = Effort_real_t.*c_ij_t.*animal_vector_norm;
        regime_norm = regime_t./sum(regime_t,1);
        regime_norm(isnan(regime_norm)) = 0;
        regime_y_mean = sum(regime_norm.*niche_plants,2);
        matching_x_y = abs(niche_consumers' - regime_y_mean);
        mean_matching_x_y = nanmean(matching_x_y, "all");        
        
        productivity_t = [productivity_t,sum(functional_response_animal_t,'all')*dz*dx];
        regime_morpho_cor_t = [regime_morpho_cor_t, mean_matching_x_y];
    end
    Productivity_t = [Productivity_t;productivity_t];
    Productivity_f = [Productivity_f;productivity_f];
    Regime_morpho_cor_t = [Regime_morpho_cor_t;regime_morpho_cor_t];
    Regime_morpho_cor_f = [Regime_morpho_cor_f;regime_morpho_cor_f];
end

%% BOXPLOT 
% average productivity (on the last 100 time steps)
Productivity_boxplot_t = mean(Productivity_t(:,end-10:end),2);
Productivity_boxplot_f = mean(Productivity_f(:,end-10:end),2);
Productivity_boxplot_mean = [Productivity_boxplot_f,Productivity_boxplot_t];
boxplot(Productivity_boxplot_mean);
title('Productivity','interpreter','latex','FontSize',18)
%% average regime-morpho corr
Regmorphocor_boxplot_t = mean(Regime_morpho_cor_t(:,end-10:end),2);
Regmorphocor_boxplot_f = mean(Regime_morpho_cor_f(:,end-10:end),2);
Regmorphocor_boxplot_mean = [Regmorphocor_boxplot_f,Regmorphocor_boxplot_t];
boxplot(Regmorphocor_boxplot_mean);
title('Cor reg-morpho','interpreter','latex','FontSize',18)

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

close
dist_productivity_animal = Productivity_boxplot_f - Productivity_boxplot_t; 
s = scatter(Z_mean,dist_productivity_animal,5,'b','filled');
s.MarkerFaceAlpha = .5;
hold on

% add mean + std
Z_mean_round = round(Z_mean,1);
Z_unique = unique(Z_mean_round); animal_vector_norm = animal_vector./sum(animal_vector, "all");

dist_productivity_animal_replica = {};
for i = 1:length(Z_unique)
    value = Z_unique(i);
    dist_productivity_animal_replica{i} = dist_productivity_animal(Z_mean_round==value);
end
dist_productivity_animal_median = cellfun(@median, dist_productivity_animal_replica);
dist_productivity_animal_mean = cellfun(@mean, dist_productivity_animal_replica);
dist_productivity_animal_std = cellfun(@std, dist_productivity_animal_replica);

x = 0:.1:1;
x2 = [x, fliplr(x)];

p = plot(Z_unique,dist_productivity_animal_mean,'--b','LineWidth',3);
set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
inBetween = [dist_productivity_animal_mean-dist_productivity_animal_std,...
    fliplr(dist_productivity_animal_mean+dist_productivity_animal_std)];
h = fill(x2, inBetween, 'b');
set(h,'facealpha',.2,'EdgeColor','b');
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hold off

ax = gca;
ax.FontSize = 16;
xlabel({'mean foraging trait of','community with AF evolution'},'interpreter','latex','FontSize',20)
ylabel({'difference in productivity'},'interpreter','latex','FontSize',20)

%% Corr regime-morpho in function of Z mean
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

close
dist_reg_morpho_cor = Regmorphocor_boxplot_f - Regmorphocor_boxplot_t; 
s = scatter(Z_mean, dist_reg_morpho_cor, 5, 'b', 'filled');
s.MarkerFaceAlpha = .5;
hold on

% add mean + std
Z_mean_round = round(Z_mean,1);
Z_unique = unique(Z_mean_round); 

dist_reg_morpho_replica = {};
for i = 1:length(Z_unique)
    value = Z_unique(i);
    dist_reg_morpho_replica{i} = dist_reg_morpho_cor(Z_mean_round==value);
end
dist_reg_morpho_cor_median = cellfun(@median, dist_reg_morpho_replica);
dist_reg_morpho_cor_mean = cellfun(@mean, dist_reg_morpho_replica);
dist_reg_morpho_cor_std = cellfun(@std, dist_reg_morpho_replica);

x = 0:.1:1;
x2 = [x, fliplr(x)];

p = plot(Z_unique,dist_reg_morpho_cor_mean,'--b','LineWidth',3);
set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
inBetween = [dist_reg_morpho_cor_mean-dist_reg_morpho_cor_std,...
    fliplr(dist_reg_morpho_cor_mean+dist_reg_morpho_cor_std)];
h = fill(x2, inBetween, 'b');
set(h,'facealpha',.2,'EdgeColor','b');
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hold off

ax = gca;
ax.FontSize = 16;
xlabel({'mean foraging trait of','community with AF evolution'},'interpreter','latex','FontSize',20)
ylabel({'difference in distance', 'between mean consumed y and x'},'interpreter','latex','FontSize',20)

