%% NICHE-JUMP PERTURBATION AVEC Z FIXE
% load('prep.mat') 

%% INITIALISATION %%
number_of_animals = 77; 
number_of_plants  = 77;
% number_of_foraging = 11;
number_of_foraging = 101;
foraging_trait = linspace(0,1,number_of_foraging);
dz = 1;
Foraging_trait = reshape(foraging_trait,1,1,number_of_foraging);
% FOR EACH FORAGING TRAIT
y0_dyn_all = cell(length(foraging_trait),1);
parfor i = 1:length(foraging_trait)

    % FONCTION K(y)
    K0 = 50; % maximal carrying capacity
    y0 = -25; % optimal trait pour deplacemnt niche %
    sigmaK = 2.5;

    % TRAITS %
    % Animal traits
    xmin = -35; %pour deplacement niche
    xmax = 35;
    xx = linspace(xmin,xmax,number_of_animals);
    dx = 1;
    traits_of_animals = xx';

    % Plants traits
    ymin = -35; % pour deplacement niche
    ymax = 35;
    yy = linspace(ymin,ymax,number_of_plants);
    dy = 1;
    traits_of_plants = yy';
    % FONCTION K(y)
    Kf = @(trait_plant) K0*exp(-(trait_plant - y0).^2./(2*sigmaK^2));
    % Carrying capacity %
    K = Kf(traits_of_plants);
    % FONCTION C(y-y0)
    Beta = 0; % Beta = 0 : symetrical competition
    sigmaC = sigmaK-1;
    C = @(trait_plant_y,trait_plant_yi) exp(sigmaC^2*Beta^2/2).*exp(-(trait_plant_y-trait_plant_yi + sigmaC^2*Beta^2).^2./(2*sigmaC^2));
    % Competition %
    Compet = C(traits_of_plants,traits_of_plants');

    [YY,XX] = meshgrid(yy,xx);
    SIGMA = .9;
    delta_ij = complementary_traits(SIGMA,XX,YY);
    delta_ij = delta_ij';

    effort_speed_of_change = .5;

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

%     DENSITIES
    animal_density = zeros(number_of_animals,number_of_foraging);
    plant_density = zeros(1,number_of_plants);
%     pour deplacement niche : on initialise en bordure de la grande niche
%     permettant de faire des sauts, on reste tout de même assez eloigne de
%     la bordure pour laisser la place a la communaute d'emerger
%     pleinement, a chaque iteration on fixe un Z different (fixe)
    animal_density(10,i) = .1; 
    plant_density(10) = .1;
    
%     EFFORT
    effort_ij0 = delta_ij;
    effort_ij0 = effort_ij0./sum(effort_ij0,2);
    effort_ij0 = reshape(effort_ij0,[],1);
    effort_ij0 = repmat(effort_ij0,number_of_foraging,1);

%     DIFFUSION MATRIX
    D = 1e-3; % Diffusion coefficient
    e = ones(number_of_animals,1);
    A_animal = spdiags([e,-2*e,e],-1:1,number_of_animals,number_of_animals);
    A_animal(1,1) = -1; A_animal(end,end) = -1;
    A_animal = D*A_animal/(dx^2);

    D = 1e-3; % Diffusion coefficient
    e = ones(number_of_plants,1);
    A_plant = spdiags([e,-2*e,e],-1:1,number_of_plants,number_of_plants);
    A_plant(1,1) = -1; A_plant(end,end) = -1;
    A_plant = D*A_plant/(dx^2);
 
    Df = 0; % Z fixe
    e = ones(number_of_foraging,1);
    A_foraging = spdiags([e,-2*e,e],-1:1,number_of_foraging,number_of_foraging);
    A_foraging(1,1) = -1; A_foraging(end,end) = -1;
    A_foraging = Df*A_foraging/(dx^2) ;

    %%%%%%% ODE %%%%%%
    animal_density = reshape(animal_density,number_of_animals*number_of_foraging,1);
    B = [animal_density ; plant_density' ; effort_ij0];
    tmax = 5000 ; % 5000
    tspan = 0 : 50 : tmax ;
%     tspan = 0 : 1 : tmax ;
    options = odeset('RelTol',1e-3,'AbsTol',1e-4);

    [t,dB]=ode45(@(t,B) demographic_system_evol_foraging_bis(t, B, number_of_animals, delta_ij, extraction_coeff, conversion_coeff,...
    animal_intrinsic_growth, animal_intraspe_compet, plant_intrinsic_growth, seuil_effort, A_foraging,...
    handling_time, number_of_plants, A_animal, A_plant, dx, dy, dz, Compet, K, seuil_abondance, effort_speed_of_change,...
    number_of_foraging, Foraging_trait, plant_intraspe_compet),tspan, B, options);

    animal_density = dB(:,1:number_of_animals*number_of_foraging);
    plant_density = dB(:,number_of_animals*number_of_foraging+1:number_of_animals*number_of_foraging+number_of_plants);
    effort_dynamic = dB(:,number_of_animals*number_of_foraging+number_of_plants+1:end);

    animal_density = animal_density.*(animal_density>seuil_abondance);
    plant_density = plant_density.*(plant_density>seuil_abondance);

    t_pre = t;

%     ISOLATE LAST TIME LAPS to initialize perturbation 
    animal_density_last_step = animal_density(end,:) ;
    animal_density_last_step = animal_density_last_step' ;
    plant_density_last_step = plant_density(end,:) ;
    plant_density_last_step = plant_density_last_step' ;
    effort_ij0_last_step = effort_dynamic(end,:);
    effort_ij0_last_step = effort_ij0_last_step';

%     POST-PERTURB
    test = 1 ;
    y0_perturb = y0 ;
    perturb_increment = 1;

    animal_density_perturb = animal_density_last_step;
    plant_density_perturb = plant_density_last_step;
    effort_ij0_perturb = effort_ij0_last_step;

    animal_total_perturb = [];
    plant_total_perturb = [];
    effort_dyn_tot_perturb = [];
    y0_dyn = [];
    while test==1    
        animal_density_perturb = animal_density_last_step;
        plant_density_perturb = plant_density_last_step;
%         Increasing 'amplitude' of y
        y0_perturb = y0_perturb + perturb_increment ; % optimal trait 
        Kf = @(trait_plant) K0*exp(-(trait_plant - y0_perturb).^2./(2*sigmaK^2));
%         Carrying capacity
        K = Kf(traits_of_plants);

        %%%%%%%%%%%%%%%% ODE %%%%%%%%%%%%%%%%   
        B = [animal_density_perturb ; plant_density_perturb ; effort_ij0_perturb]; % plant density retourné avant
        tmax = 3000; % 3000
        tspan = 0:30:tmax ;
%         tspan = 0:1:tmax ;
    
        [t,dB]=ode45(@(t,B) demographic_system_evol_foraging_bis(t, B, number_of_animals, delta_ij, extraction_coeff, conversion_coeff,...
        animal_intrinsic_growth, animal_intraspe_compet, plant_intrinsic_growth, seuil_effort, A_foraging,...
        handling_time, number_of_plants, A_animal, A_plant, dx, dy, dz, Compet, K, seuil_abondance, effort_speed_of_change,...
        number_of_foraging, Foraging_trait, plant_intraspe_compet),tspan, B, options);
    
        animal_density_perturb = dB(:,1:number_of_animals*number_of_foraging);
        plant_density_perturb = dB(:,number_of_animals*number_of_foraging+1:number_of_animals*number_of_foraging+number_of_plants);
        effort_dynamic_perturb = dB(:,number_of_animals*number_of_foraging+number_of_plants+1:end);
    
        animal_density_perturb = animal_density_perturb.*(animal_density_perturb>seuil_abondance);
        plant_density_perturb = plant_density_perturb.*(plant_density_perturb>seuil_abondance);
    
        animal_total_perturb = [animal_total_perturb ; animal_density_perturb] ;
        plant_total_perturb = [plant_total_perturb ; plant_density_perturb] ;
        effort_dyn_tot_perturb = [effort_dyn_tot_perturb ; effort_dynamic_perturb] ;
    
        animal_density_perturb = animal_density_perturb(end,:)' ;
        plant_density_perturb = plant_density_perturb(end,:)' ;
    
        y0_dyn = [y0_dyn ; y0_perturb] ;

%         determine if extinction of the whole animal community occured
        test = sum(animal_density_perturb,'all')>0;    
    end

y0_dyn_all{i} = y0_dyn;
end
 
% save(['jump_with_fixed_z','.mat'],'y0_dyn_all');
save(['jump_with_fixed_z_101','.mat'],'y0_dyn_all');

%%