clear all;clc;
%% Fait émerger un certains nombres (loop) de communautés foraging et tondeuse
%% Pour chaque itération une combinaison de paramètres est tirée dans une gamme
%% Enregistre les densités pour chaque communauté pour pouvoir ensuite les comparer (biomasse, diversité)

loop = 1000;
% sigma=.5 [.1;1] / sigmaK=2.5 ]1;4] / hmax=1 [hmin;2] / animal_intrinsic_growth=.1 [.1 .5] / sigmaC=sigmaK-1 
% hypercube latin
[X_scaled,X_normalized]=lhsdesign_modified(loop,[.1 1 .1 .1],[1 4 2 .5]);
parametres = X_scaled;

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
% dx = xx(2)-xx(1);
dx = 1;
traits_of_animals = xx';

foraging_trait = linspace(0,1,number_of_foraging);
% dz = foraging_trait(2) - foraging_trait(1);
dz = 1;
Foraging_trait = reshape(foraging_trait,1,1,number_of_foraging);

% Plants traits
ymin = -5;
ymax = 5;
yy = linspace(ymin,ymax,number_of_plants);
% dy =  yy(2)-yy(1);
dy = 1;
traits_of_plants = yy';

[YY,XX] = meshgrid(yy,xx);
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
plant_intraspe_compet = 0;

% HANDLING TIME TRADE-OFF
% alpha_h = 5; % expo
alpha_h = 1; % lineaire
% alpha_h = .5; % racine
hmin = .1;
hmax = .55;
h = @(z) hmin + (hmax - hmin).*z.^alpha_h;
handling_time = h(foraging_trait);
handling_time = reshape(handling_time,1,1,number_of_foraging);

%%%%%%%%
% PARFOR

Output_a_f = [];
Output_p_f = [];
Output_a_t = [];
Output_p_t = [];
Effort_f = [];
Effort_t = [];
Param = [];
sigma_loop = parametres(:,1);
sigmaK_loop = parametres(:,2);
hmax_loop = parametres(:,3);
animal_intrinsic_growth_loop = parametres(:,4);

for j = 1:loop/20
    idx = (j-1)*20 + (1:20);
    sigma_parfor = sigma_loop(idx);
    sigmaK_parfor = sigmaK_loop(idx);
    hmax_parfor = hmax_loop(idx);
    animal_intrinsic_growth_parfor = animal_intrinsic_growth_loop(idx);
    parfor i = 1:20
        % SIGMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         SIGMA = sigma_parfor(i);
         delta_ij = complementary_traits(SIGMA,XX,YY);
         delta_ij = delta_ij';

        % SIGMA_K %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         sigmaK = sigmaK_parfor(i);
         Kf = @(trait_plant) K0*exp(-(trait_plant - y0).^2./(2*sigmaK^2));
         sigmaC = sigmaK-1;
         C = @(trait_plant_y,trait_plant_yi) exp(sigmaC^2*Beta^2/2).*exp(-(trait_plant_y-trait_plant_yi + sigmaC^2*Beta^2).^2./(2*sigmaC^2));
         Compet = C(traits_of_plants,traits_of_plants');

        % h_max & alpha_h %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hmax = hmax_parfor(i);
        h = @(z) hmin + (hmax - hmin).*z.^alpha_h
        handling_time = h(foraging_trait);
        handling_time = reshape(handling_time,1,1,number_of_foraging);

        % animal_growth
        animal_intrinsic_growth = animal_intrinsic_growth_parfor(i);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % EFFORT
        effort_ij0 = delta_ij;
        effort_ij0 = effort_ij0./sum(effort_ij0,2);
        effort_ij0 = reshape(effort_ij0,[],1);
        effort_ij0 = repmat(effort_ij0,number_of_foraging,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % DIFFUSION MATRIX
        D =1e-3;
        e = ones(number_of_animals,1);
        A_animal = spdiags([e,-2*e,e],-1:1,number_of_animals,number_of_animals);
        A_animal(1,1) = -1; A_animal(end,end) = -1;
    %     A_animal(1,2) = -1; A_animal(end,end-1) = -1;
        A_animal = D*A_animal/(dx^2);

        D = 1e-3;
        e = ones(number_of_plants,1);
        A_plant = spdiags([e,-2*e,e],-1:1,number_of_plants,number_of_plants);
        A_plant(1,1) = -1; A_plant(end,end) = -1;
    %     A_plant(1,2) = -1; A_plant(end,end-1) = -1;
        A_plant = D*A_plant/(dx^2);

        D = 1e-3;
        e = ones(number_of_foraging,1);
        A_foraging = spdiags([e,-2*e,e],-1:1,number_of_foraging,number_of_foraging);
        A_foraging(1,1) = -1; A_foraging(end,end) = -1;
        A_foraging = D*A_foraging/(dx^2);

        % SIMU AF
        animal_density = zeros(number_of_animals,number_of_foraging);
        plant_density = zeros(1,number_of_plants);
        plant_density((number_of_plants+1)/2) = 1.1;
        animal_density((number_of_animals+1)/2,:) = .1;

        options = odeset('RelTol',1e-3,'AbsTol',1e-5);
        animal_density = reshape(animal_density,number_of_animals*number_of_foraging,1);
        B = [animal_density ; plant_density' ; effort_ij0];
        tmax = 1e3;
        tspan = 0:10:tmax;

        [t,dB]=ode45(@(t,B) demographic_system_evol_foraging_bis(t, B, number_of_animals, delta_ij, extraction_coeff, conversion_coeff,...
            animal_intrinsic_growth, animal_intraspe_compet, plant_intrinsic_growth, seuil_effort, A_foraging,...
            handling_time, number_of_plants, A_animal, A_plant, dx, dy, dz, Compet, K, seuil_abondance, effort_speed_of_change,...
            number_of_foraging, Foraging_trait, plant_intraspe_compet),tspan, B, options);

        animal_density = dB(:,1:number_of_animals*number_of_foraging);
        plant_density = dB(:,number_of_animals*number_of_foraging+1:number_of_animals*number_of_foraging+number_of_plants);
        effort_dynamic = dB(:,number_of_animals*number_of_foraging+number_of_plants+1:end);
        
        % densites trop faibles (<seuil) = 0 
        animal_density = animal_density(:,:).*(animal_density(:,:)>seuil_abondance);
        plant_density = plant_density(:,:).*(plant_density(:,:)>seuil_abondance);

        Output_a_f = [Output_a_f,animal_density];
        Output_p_f = [Output_p_f,plant_density];
        Effort_f = [Effort_f,effort_dynamic];

        % SIMU TONDEUSE 
        D = 0;
        e = ones(number_of_foraging,1);
        A_foraging = spdiags([e,-2*e,e],-1:1,number_of_foraging,number_of_foraging);
        A_foraging(1,1) = -1; A_foraging(end,end) = -1;
        A_foraging = D*A_foraging/(dx^2);

        animal_density = zeros(number_of_animals,number_of_foraging);
        plant_density = zeros(1,number_of_plants);
        plant_density((number_of_plants+1)/2) = 1.1;
        animal_density((number_of_animals+1)/2,1) = 1.1;

        animal_density = reshape(animal_density,number_of_animals*number_of_foraging,1);
        B = [animal_density ; plant_density' ; effort_ij0];
        tmax = 1e3;
        tspan = 0:10:tmax;

        [t,dB]=ode45(@(t,B) demographic_system_evol_foraging_bis(t, B, number_of_animals, delta_ij, extraction_coeff, conversion_coeff,...
            animal_intrinsic_growth, animal_intraspe_compet, plant_intrinsic_growth, seuil_effort, A_foraging,...
            handling_time, number_of_plants, A_animal, A_plant, dx, dy, dz, Compet, K, seuil_abondance, effort_speed_of_change,...
            number_of_foraging, Foraging_trait, plant_intraspe_compet),tspan, B, options);

        animal_density = dB(:,1:number_of_animals*number_of_foraging);
        plant_density = dB(:,number_of_animals*number_of_foraging+1:number_of_animals*number_of_foraging+number_of_plants);
        effort_dynamic = dB(:,number_of_animals*number_of_foraging+number_of_plants+1:end);

        animal_density = animal_density(:,:).*(animal_density(:,:)>seuil_abondance);
        plant_density = plant_density(:,:).*(plant_density(:,:)>seuil_abondance);

        Output_a_t = [Output_a_t,animal_density];
        Output_p_t = [Output_p_t,plant_density];
        Effort_t = [Effort_t,effort_dynamic];

        Param = [Param,[SIGMA;sigmaK;hmax;animal_intrinsic_growth]];
    end
    save(['community_comparison_b','.mat'],'Output_a_f','Output_p_f','Effort_f'...
    ,'Output_a_t','Output_p_t','Effort_t','Param');
end
save(['community_comparison_b','.mat'],'Output_a_f','Output_p_f','Effort_f'...
    ,'Output_a_t','Output_p_t','Effort_t','Param');

