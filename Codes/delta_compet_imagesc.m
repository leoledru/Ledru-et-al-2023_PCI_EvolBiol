%% VISUALISATION DES FONCTIONS DE COMPLEMENTARITE (DELTA) ET DE COMPETITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

number_of_animals = 100 ; % 21; 
number_of_plants  = 100 ; %21; 
number_of_foraging = 11; %11

% FONCTION K(y)
K0 = 50; % maximal carrying capacity
y0 = 0 ; %optimal trait
sigmaK = 1.5;
Kf = @(trait_plant) K0*exp(-(trait_plant - y0).^2./(2*sigmaK^2));

% FONCTION C(y-y0)
Beta = 0; % Beta = 0 : symetrical competition
sigmaC = sigmaK-1;
C = @(trait_plant_y,trait_plant_yi) exp(sigmaC^2*Beta^2/2).*exp(-(trait_plant_y-trait_plant_yi + sigmaC^2*Beta^2).^2./(2*sigmaC^2));

% TRAITS %
% Animal traits
xmin = -5 ;
xmax = 5 ;
% xmin = -35; %pour deplacement niche
% xmax = 35;
xx = linspace(xmin,xmax,number_of_animals);
%dx = xx(2)-xx(1);
dx = 1;
traits_of_animals = xx';

foraging_trait = linspace(0,1,number_of_foraging);
%dz = foraging_trait(2) - foraging_trait(1);
dz = 1;
Foraging_trait = reshape(foraging_trait,1,1,number_of_foraging);

% Plants traits
ymin = -5 ;
ymax = 5 ;
% ymin = -35; % pour deplacement niche
% ymax = 35;
yy = linspace(ymin,ymax,number_of_plants);
%dy =  yy(2)-yy(1);
dy = 1;
traits_of_plants = yy';

[YY,XX] = meshgrid(yy,xx);
SIGMA = .9;
delta_ij = complementary_traits(SIGMA,XX,YY);
delta_ij = delta_ij';

effort_speed_of_change = .5;

% Competition %
Compet = C(traits_of_plants,traits_of_plants');
% Carrying capacity %
K = Kf(traits_of_plants);

% HOMOGENEOUS INITIALISATION %
extraction_coeff = .8;         % .*ones(number_of_animals,1);
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
handling_time = h(foraging_trait);
handling_time = reshape(handling_time,1,1,number_of_foraging);

% K
subplot(2,2,1)
plot(K,'k','LineWidth',3)
xlabel('ressource niche trait y','FontSize',20,'interpreter','latex')
ylabel('carrying capacity','FontSize',20,'interpreter','latex')
xticks(0:10:100)
xticklabels(linspace(-5,5,11))
hold on
plot([50,50],[0,50],'--r','LineWidth',2)
text(51,4,'$y_0$','interpreter','latex','FontSize',20,'Color','red')
hold off

% Compet
subplot(2,2,3)
imagesc(Compet);
colormap jet
c = colorbar;
ylabel(c,'ressource competition','FontSize',20,'interpreter','latex')
xticks([1,10:10:100])
xticklabels(linspace(-5,5,11))
yticks([1,10:10:100])
yticklabels(linspace(-5,5,11))
xlabel('ressource niche trait y','FontSize',20,'interpreter','latex')
ylabel("ressource niche trait y'",'FontSize',20,'interpreter','latex')

% Interaction
subplot(2,2,4)
imagesc(delta_ij);
colormap jet
c = colorbar;
ylabel(c,{'maximum','interaction strength'},'FontSize',20,'interpreter','latex')
xticks([1,10:10:100])
xticklabels(linspace(-5,5,11))
yticks([1,10:10:100])
yticklabels(linspace(-5,5,11))
xlabel('consumer niche trait x','FontSize',20,'interpreter','latex')
ylabel("ressource niche trait y",'FontSize',20,'interpreter','latex')


%%

% K
subplot(2,2,1)

sigmaK = 2;
Kf = @(trait_plant) K0*exp(-(trait_plant - y0).^2./(2*sigmaK^2));
K = Kf(traits_of_plants);
plot(K,'k','LineWidth',2)
xlabel('ressource niche trait y','FontSize',20,'interpreter','latex')
ylabel('carrying capacity','FontSize',20,'interpreter','latex')
xticks(0:10:100)
xticklabels(linspace(-5,5,11))
hold on
plot([50,50],[0,50],'--r','LineWidth',2)
text(51,4,'$y_0$','interpreter','latex','FontSize',20,'Color','red')
text(70,40,'$\sigma_K = 2$','interpreter','latex','FontSize',20)

sigmaK = 1;
Kf = @(trait_plant) K0*exp(-(trait_plant - y0).^2./(2*sigmaK^2));
K = Kf(traits_of_plants);
plot(K,'--k','LineWidth',2)
xlabel('ressource niche trait y','FontSize',20,'interpreter','latex')
ylabel('carrying capacity','FontSize',20,'interpreter','latex')
xticks(0:10:100)
xticklabels(linspace(-5,5,11))

sigmaK = .5;
Kf = @(trait_plant) K0*exp(-(trait_plant - y0).^2./(2*sigmaK^2));
K = Kf(traits_of_plants);
plot(K,'-.k','LineWidth',2)
xlabel('resource niche trait y','FontSize',20,'interpreter','latex')
ylabel('carrying capacity','FontSize',20,'interpreter','latex')
xticks(0:10:100)
xticklabels(linspace(-5,5,11))

hold off

% Compet
subplot(2,2,2)
sigmaC = 2;
C = @(trait_plant_y,trait_plant_yi) exp(sigmaC^2*Beta^2/2).*exp(-(trait_plant_y-trait_plant_yi + sigmaC^2*Beta^2).^2./(2*sigmaC^2));
Compet = C(traits_of_plants,traits_of_plants');
plot(Compet(60,:),'k','LineWidth',2)
xticks(0:10:100)
xticklabels(linspace(-5,5,11))
xlabel('resource niche trait y','FontSize',20,'interpreter','latex')
ylabel({'competition experienced','by resource y = 1'},'FontSize',20,'interpreter','latex')
text(80,.8,'$\sigma_C = 2$','interpreter','latex','FontSize',20)
hold on 

sigmaC = 1;
C = @(trait_plant_y,trait_plant_yi) exp(sigmaC^2*Beta^2/2).*exp(-(trait_plant_y-trait_plant_yi + sigmaC^2*Beta^2).^2./(2*sigmaC^2));
Compet = C(traits_of_plants,traits_of_plants');
plot(Compet(60,:),'--k','LineWidth',2)
xticks(0:10:100)
xticklabels(linspace(-5,5,11))
xlabel('resource niche trait y','FontSize',20,'interpreter','latex')
ylabel({'competition experienced','by resource y = 1'},'FontSize',20,'interpreter','latex')

sigmaC = .5;
C = @(trait_plant_y,trait_plant_yi) exp(sigmaC^2*Beta^2/2).*exp(-(trait_plant_y-trait_plant_yi + sigmaC^2*Beta^2).^2./(2*sigmaC^2));
Compet = C(traits_of_plants,traits_of_plants');
plot(Compet(60,:),'-.k','LineWidth',2)
xticks(0:10:100)
xticklabels(linspace(-5,5,11))
xlabel('resource niche trait y','FontSize',20,'interpreter','latex')
ylabel({'competition experienced',"by resource y' = 1"},'FontSize',20,'interpreter','latex')
hold off

% Interaction
subplot(2,2,3)
SIGMA = 2;
delta_ij = complementary_traits(SIGMA,XX,YY);
delta_ij = delta_ij';
plot(delta_ij(60,:),'k','LineWidth',2)
xticks(0:10:100)
xticklabels(linspace(-5,5,11))
xlabel('ressource niche trait y','FontSize',20,'interpreter','latex')
ylabel({'maximum interaction','strength of consumer x = 1'},'FontSize',20,'interpreter','latex')
text(80,.6,'$\sigma = 2$','interpreter','latex','FontSize',20)
hold on

SIGMA = 1;
delta_ij = complementary_traits(SIGMA,XX,YY);
delta_ij = delta_ij';
plot(delta_ij(60,:),'--k','LineWidth',2)
xticks(0:10:100)
xticklabels(linspace(-5,5,11))
xlabel('ressource niche trait y','FontSize',20,'interpreter','latex')
ylabel({'maximum interaction','strength of consumer x = 1'},'FontSize',20,'interpreter','latex')

SIGMA = .5;
delta_ij = complementary_traits(SIGMA,XX,YY);
delta_ij = delta_ij';
plot(delta_ij(60,:),'-.k','LineWidth',2)
xticks(0:10:100)
xticklabels(linspace(-5,5,11))
xlabel('ressource niche trait y','FontSize',20,'interpreter','latex')
ylabel({'maximum interaction','strength of consumer x = 1'},'FontSize',20,'interpreter','latex')

hold off

% trade-off
alpha_h = 1; % lineaire
% alpha_h = .5; % racine
hmin = .1;
hmax = 0.55;
h = @(z) hmin + (hmax - hmin).*z.^alpha_h;
handling_time = h(foraging_trait);

subplot(2,2,4)
plot(foraging_trait,handling_time,'k','LineWidth',2)
ylim([0 .7])
xlabel('foraging trait z','FontSize',20,'interpreter','latex')
ylabel('searching time','FontSize',20,'interpreter','latex')
hold on
plot([0 1],[hmax hmax],'--k','LineWidth',1)
text(0.05,.61,'$s_{max} = 0.55$','interpreter','latex','FontSize',20)
hold off












