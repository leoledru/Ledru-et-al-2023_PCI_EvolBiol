%% Increasing mortality to extinction with fixed Z %%%

count_all = [];
parfor i = 1:101
    %%% INITIALISATION %%%
    number_of_animals = 11 ;
    number_of_plants  = 11 ;
    number_of_foraging = 101;

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
    foraging_trait = linspace(0,1,number_of_foraging);
    dz = 1;
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

    % DENSITIES
    animal_density = zeros(number_of_animals,number_of_foraging);
    plant_density = zeros(1,number_of_plants);
    animal_density((number_of_animals+1)/2,i) = .1; % one on niche's center
    plant_density((number_of_plants+1)/2) = .1;

    % EFFORT
    effort_ij0 = delta_ij;
    effort_ij0 = effort_ij0./sum(effort_ij0,2);
    effort_ij0 = reshape(effort_ij0,[],1);
    effort_ij0 = repmat(effort_ij0,number_of_foraging,1);
    
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
  
    Df = 0;    
    e = ones(number_of_foraging,1);
    A_foraging = spdiags([e,-2*e,e],-1:1,number_of_foraging,number_of_foraging);
    A_foraging(1,1) = -1; A_foraging(end,end) = -1;
    %         A_foraging(1,2) = -1; A_foraging(end,end-1) = -1;
    A_foraging = Df*A_foraging/(dx^2) ;


    %%%%%%%%%%%%%%%%%%%%%%%% ODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    tmax = 5000;
    tspan = 0:50:tmax;
    options = odeset('RelTol',1e-4,'AbsTol',1e-5);
    animal_density = reshape(animal_density,number_of_animals*number_of_foraging,1);
    B = [animal_density ; plant_density' ; effort_ij0];

    [t,dB]=ode45(@(t,B) demographic_system_evol_foraging_bis(t, B, number_of_animals, delta_ij, extraction_coeff, conversion_coeff,...
    animal_intrinsic_growth, animal_intraspe_compet, plant_intrinsic_growth, seuil_effort, A_foraging,...
    handling_time, number_of_plants, A_animal, A_plant, dx, dy, dz, Compet, K, seuil_abondance, effort_speed_of_change,...
    number_of_foraging, Foraging_trait, plant_intraspe_compet),tspan, B, options);

    animal_density = dB(:,1:number_of_animals*number_of_foraging);
    plant_density = dB(:,number_of_animals*number_of_foraging+1:number_of_animals*number_of_foraging+number_of_plants);
    effort_dynamic = dB(:,number_of_animals*number_of_foraging+number_of_plants+1:end);

    animal_density = animal_density.*(animal_density>seuil_abondance);
    plant_density = plant_density.*(plant_density>seuil_abondance);
    %ISOLATE LAST TIME LAPS to initialize perturbation 
    animal_density_last_step = animal_density(end,:) ;
    plant_density_last_step = plant_density(end,:) ;
    effort_ij0_last_step = effort_dynamic(end,:);

    
    % PERTURBATION
    tmax = 3000 ;
    tspan = 0:20:tmax ;
    options = odeset('RelTol',1e-6,'AbsTol',1e-8); %ATT plus petit

    % loop initialization
    animal_intrinsic_growth = .1 ; 
    plant_intrinsic_growth = .8 ;
    test_a = 1 ;
 
    animal_density_perturb_loop = animal_density_last_step ;
    plant_density_perturb_loop = plant_density_last_step ;
    effort_ij0_perturb_loop = effort_ij0_last_step ;

    animal_density_perturb = [] ;
    plant_density_perturb = [] ;
    effort_dynamic_perturb = [] ;
    animal_mortality_dynamic = [] ;
    t_post = 0 ;
    
    count = 0;
    perturb_increment_a = 0.1*animal_intrinsic_growth; % fixed increment
    while test_a==1
%         perturb_increment_a = 0.1*animal_intrinsic_growth; % accelerating increment

        % increasing animal mortality to global extinction in the community 
         animal_intrinsic_growth = animal_intrinsic_growth + perturb_increment_a ; %animal_intrinsic_growth = mortality           

         %%%%%%%%%%%%%%%% ODE %%%%%%%%%%%%%%%%%% 
        animal_density_perturb_loop = reshape(animal_density_perturb_loop,[],1);
        plant_density_perturb_loop = reshape(plant_density_perturb_loop,[],1);
        effort_ij0_perturb_loop = reshape(effort_ij0_perturb_loop,[],1);
        B = [animal_density_perturb_loop ; plant_density_perturb_loop ; effort_ij0_perturb_loop];

        [t,dB]=ode45(@(t,B) demographic_system_evol_foraging_bis(t, B, number_of_animals, delta_ij, extraction_coeff, conversion_coeff,...
        animal_intrinsic_growth, animal_intraspe_compet, plant_intrinsic_growth, seuil_effort, A_foraging,...
        handling_time, number_of_plants, A_animal, A_plant, dx, dy, dz, Compet, K, seuil_abondance, effort_speed_of_change,...
        number_of_foraging, Foraging_trait, plant_intraspe_compet),tspan, B, options);

        animal_density_perturb_loop = dB(:,1:number_of_animals*number_of_foraging);
        plant_density_perturb_loop = dB(:,number_of_animals*number_of_foraging+1:number_of_animals*number_of_foraging+number_of_plants);
        effort_dynamic_perturb_loop = dB(:,number_of_animals*number_of_foraging+number_of_plants+1:end);

        animal_density_perturb_loop = animal_density_perturb_loop.*(animal_density_perturb_loop>seuil_abondance);
        plant_density_perturb_loop = plant_density_perturb_loop.*(plant_density_perturb_loop>seuil_abondance);

        %stocker les moyennes des densites perturb
%         animal_density_perturb = [animal_density_perturb ; mean(animal_density_perturb_loop,1)];
%         plant_density_perturb = [plant_density_perturb ; mean(plant_density_perturb_loop,1)];
%         effort_dynamic_perturb = [effort_dynamic_perturb ; mean(effort_dynamic_perturb_loop,1)] ;
%         animal_mortality_dynamic = [animal_mortality_dynamic; animal_intrinsic_growth] ;

        % Redefinir conditions initiales
        animal_density_perturb_loop = animal_density_perturb_loop(end,:)';
        plant_density_perturb_loop = plant_density_perturb_loop(end,:)';
        effort_ij0_perturb_loop = effort_dynamic_perturb_loop(end,:)' ;

        %determine if global extinction occured
        test_a = any(animal_density_perturb_loop);
        count = count + 1; % compteur du nb d'itération avant crash
    end
count_all(i) = count;
end
save(['increase_morta_z_fixed_101_b','.mat'],'count_all')


%% POSTPROCESS
% load('increase_morta_z_fixed_101.mat')
% 
% animal_intrinsic_growth = .1;
% animal_intrinsic_growth_all = [];

% accelerating increment
% for i = 1:max(count_all)
%     perturb_increment_a = 0.1*animal_intrinsic_growth;
%     animal_intrinsic_growth = animal_intrinsic_growth + perturb_increment_a;
%     animal_intrinsic_growth_all = [animal_intrinsic_growth_all,animal_intrinsic_growth];
% end

% fixed increment
% perturb_increment_a = 0.2*animal_intrinsic_growth;
% for i = 1:max(count_all)
%     animal_intrinsic_growth = animal_intrinsic_growth + perturb_increment_a;
%     animal_intrinsic_growth_all = [animal_intrinsic_growth_all,animal_intrinsic_growth];
% end
% 
% 
% max_morta = [];
% for i = 1:length(count_all)
%     max_morta(i) = animal_intrinsic_growth_all(count_all(i));
% end
% foraging_trait = linspace(0,1,101);
% scatter(foraging_trait,max_morta,'pk','filled')
% ax = gca;
% ax.FontSize = 16;
% xlabel('foraging trait','interpreter','latex','FontSize',25)
% ylabel({'maximum sustainable','consumer mortality'},'interpreter','latex','FontSize',25)





