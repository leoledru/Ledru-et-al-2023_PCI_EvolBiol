%% POSTPROCESS NICHE OVERLAP Z FIXED

clear all; clc; close all;
Color = get(gca,'colororder');

% load('../Data/effect_of_z.mat')
load('../Data/community_comparison_all.mat','Effort_f','Param','Output_a_f','Output_p_f','Output_p_t'); 
% load('../Data/community_comparison_K4.mat','Effort_f','Param','Output_a_f','Output_p_f','Output_p_t'); 

%% INIT
number_of_animals = 11 ; % 21; 
number_of_plants  = 11 ; %21; 
number_of_foraging = 11; %11
foraging_trait = linspace(0,1,number_of_foraging);

np = size(Param,2);
ntt = 101;
Plant_density  = reshape(Output_p_f,ntt,number_of_plants,np);
Effort         = reshape(Effort_f,ntt,number_of_plants,number_of_animals*number_of_foraging,size(Param,2));
Animal_density = reshape(Output_a_f,ntt,number_of_animals,number_of_foraging,np);
Plant_density_t  = reshape(Output_p_t,ntt,number_of_plants,np);
%% Compute z average for color in scatters below
nt = 10;
Z_mean = zeros(1,np);
Z_var  = zeros(1,np);
parfor ip = 1:np
    animal = Animal_density(end-nt+1:end,:,:,ip);
    a_sum_and_norm = permute(sum(animal,2)./sum(animal,[2,3]),[1,3,2]);
    z_mean    = sum(a_sum_and_norm.*foraging_trait,2);
    z_var  = sum(a_sum_and_norm.*(foraging_trait-z_mean).^2,2);
    Z_mean(ip) = mean(z_mean);
    Z_var(ip)  = mean(z_var);
end

%% INITIALISATION %%%    
% FONCTION K(y)
K0 = 50; % maximal carrying capacity
y0 = 0 ; %optimal trait
sigmaK = 2.5;
Kf = @(trait_plant) K0*exp(-(trait_plant - y0).^2./(2*sigmaK^2));

% FONCTION C(y-y0)
Beta = 0; % Beta = 0 : symetrical competition
sigmaC = sigmaK-1;
C = @(trait_plant_y,trait_plant_yi) exp(sigmaC^2*Beta^2/2).*exp(-(trait_plant_y-trait_plant_yi + sigmaC^2*Beta^2).^2./(2*sigmaC^2));

% TRAITS %
% Animal traits
xmin = -5 ;
xmax = 5 ;
xx = linspace(xmin,xmax,number_of_animals);
dx = 1;
traits_of_animals = xx';

% foraging_trait = linspace(0,1,number_of_foraging);
%dz = foraging_trait(2) - foraging_trait(1);
dz = .1;
Foraging_trait = reshape(foraging_trait,1,1,number_of_foraging);

% Plants traits
ymin = -5 ;
ymax = 5 ;
yy = linspace(ymin,ymax,number_of_plants);
dy = 1;
traits_of_plants = yy';

[YY,XX] = meshgrid(yy,xx);
SIGMA = .9;
delta_ij = complementary_traits(SIGMA,XX,YY);
delta_ij = delta_ij';
Delta_ij =repmat(delta_ij,1,number_of_foraging);

% HOMOGENEOUS INITIALISATION %
extraction_coeff = .5; %.8;         % .*ones(number_of_animals,1);
conversion_coeff = .3;         % .*ones(number_of_animals,1);
animal_intrinsic_growth = .1; % .*ones(number_of_animals,1);
animal_intraspe_compet  = .01; % .*ones(number_of_animals,1);
plant_intrinsic_growth  = .8; % .5.*ones(number_of_plants,1);
seuil_abondance = 1e-5;
seuil_effort = seuil_abondance;
plant_intraspe_compet = 0.01;

% HANDLING TIME TRADE-OFF
% alpha_h = 5; % expo
alpha_h = 1; % lineaire
% alpha_h = .5; % racine
hmin = .1;
hmax = 0.55;
h = @(z) hmin + (hmax - hmin).*z.^alpha_h;
hz = h(foraging_trait);
Hz = reshape(repmat(hz,number_of_animals,1),1,number_of_foraging*number_of_animals);
zz = repmat(foraging_trait,number_of_animals,1);
zz = zz(:)';
%% COMPUTE RHO
RHO   = zeros(1,np);
RHO_t = zeros(1,np);
% ZZ = zeros(1,1,length(zz));
% ZZ(1,1,:) = zz;
% D_ij = permute(repmat(Delta_ij,1,1,nt),[3,1,2]);
sigma_loop  = Param(1,:);
sigmaK_loop = Param(2,:);
hmax_loop   = Param(3,:);
animal_intrinsic_growth_loop = Param(4,:);
parfor ip = 1:np
    %% SIGMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    SIGMA = sigma_loop(ip);
    delta_ij = complementary_traits(SIGMA,XX,YY);
    delta_ij = delta_ij';
    Delta_ij =repmat(delta_ij,1,number_of_foraging);
    
    %% SIGMA_K %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sigmaK = sigmaK_loop(ip);
    Kf = @(trait_plant) K0*exp(-(trait_plant - y0).^2./(2*sigmaK^2));
    % Kf = @(trait_plant) K0*exp(-(trait_plant - y0).^4./(12*sigmaK^4));
    
    sigmaC = sigmaK-1;
    C = @(trait_plant_y,trait_plant_yi) exp(sigmaC^2*Beta^2/2).*exp(-(trait_plant_y-trait_plant_yi + sigmaC^2*Beta^2).^2./(2*sigmaC^2));
    Compet = C(traits_of_plants,traits_of_plants');
    
    %% h_max & alpha_h %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hmax = hmax_loop(ip);
    h = @(z) hmin + (hmax - hmin).*z.^alpha_h;
    handling_time = h(foraging_trait);
    handling_time = reshape(handling_time,1,1,number_of_foraging);
    hz = h(foraging_trait);
    Hz = reshape(repmat(hz,number_of_animals,1),1,number_of_foraging*number_of_animals);

    % animal_growth
    animal_intrinsic_growth = animal_intrinsic_growth_loop(ip);
    
    %% AF evolution
    animal_density = Animal_density(end,:,:,ip);
    plant_density  = permute(Plant_density(end,:,ip),[2,3,1]);
    effort_AF  = permute(Effort(end,:,:,ip),[2,3,4,1]);
    effort_t   = plant_density./(sum(plant_density,1) + (sum(plant_density,1)==0) );
    effort_iAF = effort_AF.*zz+(1-zz).*effort_t;
    
    ci  = effort_iAF.*Delta_ij;
    cci = ci./(1+Hz.*extraction_coeff.*sum(ci.*plant_density,1));
    Cci = sum(reshape(cci,number_of_plants,number_of_animals,number_of_foraging).*animal_density,3)./sum(animal_density,3);

    al = Compet*plant_density./(plant_density.*Kf(traits_of_plants));
    
    % cci_mean = sum(ci.*plant_density,1)./(sum(plant_density,1) + (sum(plant_density,1)==0) );
    % Ci  = ((cci-cci_mean).*plant_density)'*(cci-cci_mean);
    % Ci  = (Cci./(plant_intrinsic_growth*al))'*(Cci);
    Ci  = (Cci)'*(Cci);

    norm_Ci = diag(Ci);
    N_Ci = sqrt(norm_Ci*norm_Ci');
    CCi = (Ci-diag(norm_Ci))./N_Ci;
    Rho = mean(CCi(CCi>0));
 
    % ui = effort_iAF.* plant_density.*Delta_ij;
    % dui = 1+Hz.*extraction_coeff.*sum(ui,1);
    % uui = ui./dui;
    % Ui = uui'*uui;
    % norm_Ui = diag(Ui);
    % N_Ui = sqrt(norm_Ui*norm_Ui');
    % UUi = (Ui-diag(norm_Ui))./N_Ui;
    % Rho = mean(UUi(UUi>0));
    RHO(ip) = Rho;
    
    %% Random Foraging
    plant_density = permute(Plant_density_t(end,:,ip),[2,3,1]);
    effort_t      = plant_density./(sum(plant_density,1) + (sum(plant_density,1)==0) );
    effort_iAF    = effort_t;
    ui = effort_iAF.*Delta_ij;
    dui = 1+Hz.*extraction_coeff.*sum(ui.* plant_density,1);
    uui = ui./dui;
    al = Compet*plant_density./(plant_density.*Kf(traits_of_plants));

    % uui_mean = sum(ui.*plant_density,1)./(sum(plant_density,1) + (sum(plant_density,1)==0) );
    % Ui = ((uui-uui_mean).* plant_density)'*(uui-uui_mean);
    % Ui = (uui./(plant_intrinsic_growth.*al))'*(uui);
    Ui = (uui)'*(uui);

    norm_Ui = diag(Ui);
    N_Ui = sqrt(norm_Ui*norm_Ui');
    UUi = (Ui-diag(norm_Ui))./N_Ui;
    Rho = mean(UUi(UUi>0));
    RHO_t(ip) = Rho;
end

DRHO = (RHO-RHO_t)./RHO_t*100;

%% FIGURE DIFFERENCE NICHE OVERLAP
% AJOUT MEAN + SHADE STD
Z_mean_round = round(Z_mean,1); % round pour réduire le nombre de Z différents
% et avoir plusieurs communautés (réplicas)
% pour chaque Z
Z_unique = unique(Z_mean_round);
% RESOURCE
DRHO_replica = {};
for i = 1:length(Z_unique)
    value = Z_unique(i);
    DRHO_replica{i} = DRHO(Z_mean_round==value);
end
%% 
DRHO_replica_med   = cellfun(@median, DRHO_replica);
DRHO_replica_mean   = cellfun(@mean, DRHO_replica);
DRHO_replica_q_inf = cellfun(@(x) quantile(x,0.25), DRHO_replica);
DRHO_replica_q_sup = cellfun(@(x) quantile(x,0.75), DRHO_replica);

Z = Z_unique;
ZZ = [Z, fliplr(Z)];

%%% Linear regression
[Z_sort,Isort] = sort(Z_mean);
% X = [ones(length(Z_sort),1),Z_sort'];
% b = X\(RHO(Isort)');
% Rho_lr = X*b;
[p,s] = polyfit(Z_sort,DRHO(Isort),1);
[Rho_lr,dRho] = polyval(p,Z_sort,s);
Rho_inf = Rho_lr-2*dRho;
Rho_sup = Rho_lr+2*dRho;

figure(1)
clf
hold on
sc = scatter(Z_mean,DRHO,5,'d','filled','MarkerFaceColor',Color(1,:));
sc.MarkerFaceAlpha = .5;

line([0,1],[0,0],'color','k','linewidth',2)
%%% LINEAR REGRESSION
% pr = plot(Z_sort,Rho_lr,'-','LineWidth',2,'Color',Color(1,:))%[0.6,0.6,0.6]);
% set(get(get(pr,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% inBetween = [Rho_inf,fliplr(Rho_sup)];
% h = fill([Z_sort,fliplr(Z_sort)], inBetween, Color(2,:));
% set(h,'facealpha',.2,'EdgeColor',Color(2,:));

%%% MEDIAN
pmed = plot(Z_unique,DRHO_replica_med,'--','LineWidth',3,'Color',Color(1,:));
    set(get(get(pmed,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    inBetween = [DRHO_replica_q_inf,fliplr(DRHO_replica_q_sup)];
    h = fill(ZZ, inBetween, Color(1,:));
    set(h,'facealpha',.2,'EdgeColor',Color(1,:));
    
%%% MEAN    
% pm = plot(Z_unique,DRHO_replica_mean,'--','LineWidth',3,'Color',Color(3,:));

ymax = max(DRHO_replica_q_sup)*1.1;
ymin = max(min(DRHO_replica_q_inf)*1.1,-100);
    ylim([ymin,ymax])

set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
hold off
ax = gca;
ax.FontSize = 16;
grid on
ax.GridLineStyle = '--';
xlabel({'mean foraging trait of','the system with PF evolution'},'interpreter','latex','FontSize',20)
ylabel({'Percent difference'; 'in niche overlap'},'interpreter','latex','FontSize',20,'color','k')


%% FIGURE NICHE OVERLAP
% % AJOUT MEAN + SHADE STD
% Z_mean_round = round(Z_mean,1); % round pour réduire le nombre de Z différents
% % et avoir plusieurs communautés (réplicas)
% % pour chaque Z
% Z_unique = unique(Z_mean_round);
% % RESOURCE
% RHO_replica = {};
% for i = 1:length(Z_unique)
%     value = Z_unique(i);
%     RHO_replica{i} = RHO(Z_mean_round==value);
% end
% %% 
% RHO_replica_med   = cellfun(@median, RHO_replica);
% RHO_replica_mean   = cellfun(@mean, RHO_replica);
% RHO_replica_q_inf = cellfun(@(x) quantile(x,0.05), RHO_replica);
% RHO_replica_q_sup = cellfun(@(x) quantile(x,0.95), RHO_replica);
% 
% Z = Z_unique;
% ZZ = [Z, fliplr(Z)];
% 
% %%% Linear regression
% [Z_sort,Isort] = sort(Z_mean);
% % X = [ones(length(Z_sort),1),Z_sort'];
% % b = X\(RHO(Isort)');
% % Rho_lr = X*b;
% [p,s] = polyfit(Z_sort,RHO(Isort),1);
% [Rho_lr,dRho] = polyval(p,Z_sort,s);
% Rho_inf = Rho_lr-2*dRho;
% Rho_sup = Rho_lr+2*dRho;
% 
% figure(1)
% clf
% hold on
% sc = scatter(Z_mean,RHO,5,'d','filled','MarkerFaceColor',Color(1,:));
% sc.MarkerFaceAlpha = .5;
% 
% % plot(Z_sort,Rho_lr,'linewidth',3)
% %%% LINEAR REGRESSION
% pr = plot(Z_sort,Rho_lr,'-','LineWidth',2,'Color',Color(1,:))%[0.6,0.6,0.6]);
% % set(get(get(pr,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% % inBetween = [Rho_inf,fliplr(Rho_sup)];
% % h = fill([Z_sort,fliplr(Z_sort)], inBetween, Color(2,:));
% % set(h,'facealpha',.2,'EdgeColor',Color(2,:));
% 
% %%% MEDIAN
% pmed = plot(Z_unique,RHO_replica_med,'--','LineWidth',3,'Color',Color(1,:));
%     set(get(get(pmed,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%     inBetween = [RHO_replica_q_inf,fliplr(RHO_replica_q_sup)];
%     h = fill(ZZ, inBetween, Color(1,:));
%     set(h,'facealpha',.2,'EdgeColor',Color(1,:));
%     
% %%% MEAN    
% % pm = plot(Z_unique,RHO_replica_mean,'--','LineWidth',3,'Color',Color(3,:));
% 
% %     ymax = max(dist_biomasse_resource_quantile_sup)*1.1;
% %     ymin = min(dist_biomasse_resource_quantile_inf)*1.1;
% %     ylim([ymin,ymax])
% 
% set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
% hold off
% ax = gca;
% ax.FontSize = 16;
% xlabel({'foraging trait $z$'},'interpreter','latex','FontSize',20)
% ylabel('Niche overlap $\rho$','interpreter','latex','FontSize',20,'color','k')
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
