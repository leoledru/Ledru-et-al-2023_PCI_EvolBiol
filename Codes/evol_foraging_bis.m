%% CODE PERMETTANT DE FAIRE UNE SIMU UNIQUE POUR DIVERS TESTS
%% UN BLOCK DE CODE PERMET LA VISUALISATION DYNAMIQUE (POUR CRÉER UN GIF PAR EXEMPLE)

% clear all; clc;

%%% Version avec la méthode de calcul des efforts de Jimmy (matrice de
% transition Q). Ici pour les efforts et delta, les plantes sont en lignes
% et les animaux en colonnes
tic

number_of_animals = 31;
number_of_plants  = 31;
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
dx = 1;
traits_of_animals = xx';

foraging_trait = linspace(0,1,number_of_foraging);
dz = foraging_trait(2) - foraging_trait(1);
dz = 1;
Foraging_trait = reshape(foraging_trait,1,1,number_of_foraging);

% Plants traits
ymin = -5;
ymax = 5;
yy = linspace(ymin,ymax,number_of_plants);
dy =  yy(2)-yy(1);
dy = 1;
traits_of_plants = yy';

[YY,XX] = meshgrid(yy,xx);
SIGMA = .3;
delta_ij = complementary_traits(SIGMA,XX,YY);
delta_ij = delta_ij';

effort_speed_of_change = .5;

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
plant_intraspe_compet = 0.01;

% HANDLING TIME TRADE-OFF
% alpha_h = 5; % expo
alpha_h = 1; % lineaire
% alpha_h = .5; % racine
hmin = .1;
hmax = 0.7;
h = @(z) hmin + (hmax - hmin).*z.^alpha_h;
handling_time = h(foraging_trait);
handling_time = reshape(handling_time,1,1,number_of_foraging);

% DENSITIES

% animal_density = (xx<3*dx).*(xx>-3*dx);
% plant_density  = (yy<3*dy).*(yy>-3*dy);
animal_density = zeros(number_of_animals,number_of_foraging);
plant_density = zeros(1,number_of_plants);
animal_density((number_of_animals+1)/2,1) = .1;
plant_density((number_of_plants+1)/2) = .1;

% animal_density(:,:) = 1; % animal everywhere
% plant_density(:) = 1; % plant everywhere
% animal_density(15,:) = 1;
% plant_density(:) = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EFFORT

effort_ij0 = delta_ij;
effort_ij0 = effort_ij0./sum(effort_ij0,2);
effort_ij0 = reshape(effort_ij0,[],1);
effort_ij0 = repmat(effort_ij0,number_of_foraging,1);
% effort_ij = reshape(effort_ij,number_of_animals,number_of_plants,number_of_foraging);
% effort_ij = zeros(number_of_animals,number_of_plants);
% effort_ij(animal_density>0,plant_density>0) = 1;
% effort_ij(:,:) = effort_ij(:,:)./sum(effort_ij(:,:),2);
% effort_ij(isnan(effort_ij)) = 0;

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIFFUSION MATRIX

D = 1e-3; % Diffusion coefficient
%         D = 0;
e = ones(number_of_animals,1);
A_animal = spdiags([e,-2*e,e],-1:1,number_of_animals,number_of_animals);
A_animal(1,1) = -1; A_animal(end,end) = -1;
%         A_animal(1,2) = -1; A_animal(end,end-1) = -1;
A_animal = D*A_animal/(dx^2);

D = 1e-3; % Diffusion coefficient
%         D = 0;
e = ones(number_of_plants,1);
A_plant = spdiags([e,-2*e,e],-1:1,number_of_plants,number_of_plants);
A_plant(1,1) = -1; A_plant(end,end) = -1;
%         A_plant(1,2) = -1; A_plant(end,end-1) = -1;
A_plant = D*A_plant/(dx^2);

D = 1e-3; % Diffusion coefficient
%         D = 0;
e = ones(number_of_foraging,1);
A_foraging = spdiags([e,-2*e,e],-1:1,number_of_foraging,number_of_foraging);
A_foraging(1,1) = -1; A_foraging(end,end) = -1;
%         A_foraging(1,2) = -1; A_foraging(end,end-1) = -1;
A_foraging = D*A_foraging/(dx^2);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    % ODE %
                    
% Continuer une simu avec EFFORT %
% animal_density = reshape(animal_density(end,:),number_of_animals,number_of_foraging);
% plant_density = plant_density(end,:);
% effort_ij0 = effort_dynamic(end,:);                            
% effort_ij0 = reshape(effort_ij0,[],1);

% couper les mutations 
% A_foraging = A_foraging.*zeros(1);
% A_plant = A_plant.*zeros(1);
% A_animal = A_animal.*zeros(1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% AVEC DYNAMIQUE DES EFFORTS

% options = odeset('RelTol',1e-3,'AbsTol',1e-4);

% Effort_ij0 = reshape(effort_ij,[],1);
animal_density = reshape(animal_density,number_of_animals*number_of_foraging,1);
B = [animal_density ; plant_density' ; effort_ij0];

tmax = 5000;
tspan = [0:10:tmax];

% bis
options = odeset('RelTol',1e-3,'AbsTol',1e-4);

[t,dB]=ode45(@(t,B) demographic_system_evol_foraging_bis(t, B, number_of_animals, delta_ij, extraction_coeff, conversion_coeff,...
animal_intrinsic_growth, animal_intraspe_compet, plant_intrinsic_growth, seuil_effort, A_foraging,...
handling_time, number_of_plants, A_animal, A_plant, dx, dy, dz, Compet, K, seuil_abondance, effort_speed_of_change,...
number_of_foraging, Foraging_trait, plant_intraspe_compet), tspan , B, options); % tspan


animal_density = dB(:,1:number_of_animals*number_of_foraging);
plant_density = dB(:,number_of_animals*number_of_foraging+1:number_of_animals*number_of_foraging+number_of_plants);
effort_dynamic = dB(:,number_of_animals*number_of_foraging+number_of_plants+1:end);

animal_density = animal_density.*(animal_density>seuil_abondance);
plant_density = plant_density.*(plant_density>seuil_abondance);

elapsedTime = toc;

% save(['new_method_with_index','.mat'],'animal_density','plant_density','effort_dynamic','elapsedTime')

%% FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%
%% VISUALISATION DYNAMIQUE 
animal = [];

for i=1:size(plant_density,1)
    if mod(i,10)==0
        drawnow
    
        subplot(1,2,1)
        plot(plant_density(i,:),'g');
    %     plot(yy,plant_density_oldmethod(i,:),'g');

    
        subplot(1,2,2)
        imagesc(reshape(animal_density(i,:),number_of_animals,number_of_foraging));
        colorbar
    %     imagesc(reshape(animal_density_oldmethod(i,:),number_of_animals,number_of_foraging));
    end
end

%% Figure for example of community emergence

animal_niche = [];
animal_z = [];
for i = 1:size(animal_density,1)
    animal = reshape(animal_density(i,:),31,11);
    animal_niche = [animal_niche;sum(animal,2)'];
    
    animal = reshape(animal_density(i,:),31,11);
    animal_z = [animal_z;sum(animal,1)];
end

subplot(3,1,1)
imagesc(plant_density(1:100,:)')
colormap pink
c = colorbar;
ylabel(c,{'resources','density'},'interpreter','latex','FontSize',20)
yticks(1:3:31)
yticklabels(linspace(-5,5,11))
xlabel('time','interpreter','latex','FontSize',20)
ylabel('niche trait','interpreter','latex','FontSize',20)

subplot(3,1,2)
imagesc(animal_niche(1:100,:)')
colormap pink
c = colorbar;
ylabel(c,{'consumers','density'},'interpreter','latex','FontSize',20)
yticks(1:3:31)
yticklabels(linspace(-5,5,11))
xlabel('time','interpreter','latex','FontSize',20)
ylabel('niche trait','interpreter','latex','FontSize',20)

subplot(3,1,3)
imagesc(animal_z(1:100,:)')
colormap pink
c = colorbar;
ylabel(c,{'consumers','density'},'interpreter','latex','FontSize',20)
yticks(1:1:11)
yticklabels(linspace(0,1,11))
xlabel('time','interpreter','latex','FontSize',20)
ylabel('foraging trait','interpreter','latex','FontSize',20)

























