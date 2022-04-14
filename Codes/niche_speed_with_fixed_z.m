%% NICHE-SPEED PERTURBATION AVEC Z FIXE

%% INITIALISATION
number_of_animals = 77;
number_of_plants  = 77;
% number_of_foraging = 11;
number_of_foraging = 101;

% TRAITS %
% Animal traits
xmin = -35;
xmax = 35;
xx = linspace(xmin,xmax,number_of_animals);
dx = 1;
traits_of_animals = xx';

foraging_trait = linspace(0,1,number_of_foraging);
dz = 1;
Foraging_trait = reshape(foraging_trait,1,1,number_of_foraging);

% Plants traits
ymin = -35;
ymax = 35;
yy = linspace(ymin,ymax,number_of_plants);
dy = 1;
traits_of_plants = yy';
traits_of_plants_b = linspace(-70,70,77*2-1)';

[YY,XX] = meshgrid(yy,xx);
SIGMA = .5;
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

% Set the niche center and shift K according to the new center
% y0 = 35;
% K = circshift(K,y0);

y0_dyn_all = cell(length(foraging_trait),1);
parfor i = 1:length(foraging_trait)   
    % FONCTION K(y)
    K0 = 50; % maximal carrying capacity
    y0 = 0; % optimal trait
    sigmaK = 2.5;
    Kf = @(trait_plant) K0*exp(-(trait_plant - y0).^2./(2*sigmaK^2));
    % FONCTION C(y-y0)
    Beta = 0; % Beta = 0 : symetrical competition
    sigmaC = sigmaK-1;
    C = @(trait_plant_y,trait_plant_yi) exp(sigmaC^2*Beta^2/2).*exp(-(trait_plant_y-trait_plant_yi + sigmaC^2*Beta^2).^2./(2*sigmaC^2));
    % HANDLING TIME TRADE-OFF
    % alpha_h = 5; % expo
    alpha_h = 1; % lineaire
    % alpha_h = .5; % racine
    hmin = .1;
    hmax = 0.55;
    h = @(z) hmin + (hmax - hmin).*z.^alpha_h;
    handling_time = h(foraging_trait);
    handling_time = reshape(handling_time,1,1,number_of_foraging);
    % COMPETITION %
    Compet = C(traits_of_plants,traits_of_plants');
    % CARRYING CAPACITY %
    K = Kf(traits_of_plants);
    
    delta_ij = complementary_traits(SIGMA,XX,YY);
    delta_ij = delta_ij';
    % Wrapping delta_ij
    delta_ij(delta_ij < 1e-10) = 0;
    for ii = 1:size(delta_ij,1)
        if delta_ij(ii,1)>0
            index = find(round(delta_ij(ii,:),12)==round(delta_ij(ii,1),12));
            values = delta_ij(ii,index(end)+1:end);
            values(values==0) = [];
            values = flip(values);
            delta_ij(ii,end-length(values)+1:end) = values;
        end
        if delta_ij(ii,end)>0
            index = find(round(delta_ij(ii,:),12)==round(delta_ij(ii,end),12));
            values = delta_ij(ii,1:index(1)-1);
            values(values==0) = [];
            values = flip(values);
            delta_ij(ii,1:length(values)) = values;
        end
    end
    % Wrapping Compet
    Compet(Compet < 1e-10) = 0;
    for ii = 1:size(Compet,1)
        if Compet(ii,1)>0
            index = find(round(Compet(ii,:),12)==round(Compet(ii,1),12));
            values = Compet(ii,index(end)+1:end);
            values(values==0) = [];
            values = flip(values);
            Compet(ii,end-length(values)+1:end) = values;
        end
        if Compet(ii,end)>0
            index = find(round(Compet(ii,:),12)==round(Compet(ii,end),12));
            values = Compet(ii,1:index(1)-1);
            values(values==0) = [];
            values = flip(values);
            Compet(ii,1:length(values)) = values;
        end
    end
    
    % DENSITIES
    animal_density = zeros(number_of_animals,number_of_foraging);
    plant_density = zeros(1,number_of_plants);
    animal_density((number_of_animals+1)/2,i) = .1; % one on niche's center, FIXED Z
    plant_density((number_of_plants+1)/2) = .1;
    % EFFORT
    effort_ij0 = delta_ij;
    effort_ij0 = effort_ij0./sum(effort_ij0,2);
    effort_ij0 = reshape(effort_ij0,[],1);
    effort_ij0 = repmat(effort_ij0,number_of_foraging,1);
    % DIFFUSION MATRIX
    D = 1e-3; % Diffusion coefficient
    e = ones(number_of_animals,1);
    A_animal = spdiags([e,-2*e,e],-1:1,number_of_animals,number_of_animals);
    % A_animal(1,1) = -1; A_animal(end,end) = -1; % no wrapping
    A_animal(1,end) = 1; A_animal(end,1) = 1; % wrapping niche mutation
    A_animal = D*A_animal/(dx^2);
    D = 1e-3; % Diffusion coefficient
    e = ones(number_of_plants,1);
    A_plant = spdiags([e,-2*e,e],-1:1,number_of_plants,number_of_plants);
    % A_plant(1,1) = -1; A_plant(end,end) = -1; % no wrapping
    A_plant(1,end) = 1; A_plant(end,1) = 1; % wrapping niche mutation
    A_plant = D*A_plant/(dx^2);    
    Df = 0; % FIXED Z 
    e = ones(number_of_foraging,1);
    A_foraging = spdiags([e,-2*e,e],-1:1,number_of_foraging,number_of_foraging);
    A_foraging(1,1) = -1; A_foraging(end,end) = -1;
    A_foraging = Df*A_foraging/(dx^2) ;

%     ODE          
    increment = 0; % the niche don't move for the initialisation
    animal_density = reshape(animal_density,number_of_animals*number_of_foraging,1);
    B = [animal_density ; plant_density' ; effort_ij0; y0];
    tmax = 5000;
    tspan = 0:50:tmax;
    % options = odeset('RelTol',1e-3,'AbsTol',1e-5,'Events',@myEvent);
    options = odeset('RelTol',1e-4,'AbsTol',1e-5);

    [t,dB] = ode45(@(t,B) demographic_system_evol_foraging_bis_perturb(t, B, number_of_animals, delta_ij, extraction_coeff, conversion_coeff,...
    animal_intrinsic_growth, animal_intraspe_compet, plant_intrinsic_growth, A_foraging,...
    handling_time, number_of_plants, A_animal, A_plant, dx, dy, dz, Compet, seuil_abondance,...
    effort_speed_of_change, number_of_foraging, Foraging_trait, traits_of_plants_b, traits_of_plants, K0, sigmaK, increment), tspan, B, options) ;

    animal_density = dB(:,1:number_of_animals*number_of_foraging);
    plant_density = dB(:,number_of_animals*number_of_foraging+1:number_of_animals*number_of_foraging+number_of_plants);
    effort_dynamic = dB(:,number_of_animals*number_of_foraging+number_of_plants+1:end-1);
    y0_dyn = dB(:,end);
    animal_density = animal_density.*(animal_density>seuil_abondance);
    plant_density = plant_density.*(plant_density>seuil_abondance);

%     Isolate last time laps to initialize perturbation 
    % DENSITIES
    animal_density_last_step = animal_density(end,:) ;
    animal_density_last_step = animal_density_last_step' ;
    plant_density_last_step = plant_density(end,:) ;
    plant_density_last_step = plant_density_last_step' ;
%     EFFORT
    effort_ij0_last_step = effort_dynamic(end,:);
    effort_ij0_last_step = effort_ij0_last_step';
    
    %%%%%%%%%%%%%
    % PERURBATION
    test = 1;
    increment_loop = linspace(.05,.2,20);
    animal_total_perturb = [];
    plant_total_perturb = [];
    effort_dynamic_perturb_tot = [];
    t_post_tot = 0 ;
    ji = 0;
    y0_dyn_tot = [];
    while test == 1
        ji = ji+1;
        increment = increment_loop(ji);
        animal_density_perturb = animal_density_last_step;
        plant_density_perturb = plant_density_last_step;
        effort_ij0_perturb = effort_ij0_last_step;

        %%%%%%%%%%%%%%%% ODE %%%%%%%%%%%%%%%%%% 
        B = [animal_density_perturb ; plant_density_perturb ; effort_ij0_perturb; y0];       
        tmax = 2000;
        tspan = 0:20:tmax;
        [t, dB] = ode45(@(t,B) demographic_system_evol_foraging_bis_perturb(t, B, number_of_animals, delta_ij, extraction_coeff, conversion_coeff,...
        animal_intrinsic_growth, animal_intraspe_compet, plant_intrinsic_growth, A_foraging,...
        handling_time, number_of_plants, A_animal, A_plant, dx, dy, dz, Compet, seuil_abondance,...
        effort_speed_of_change, number_of_foraging, Foraging_trait, traits_of_plants_b, traits_of_plants, K0, sigmaK, increment), tspan, B, options) ;

        animal_density_perturb = dB(:,1:number_of_animals*number_of_foraging);
        plant_density_perturb = dB(:,number_of_animals*number_of_foraging+1:number_of_animals*number_of_foraging+number_of_plants);
        effort_dynamic_perturb = dB(:,number_of_animals*number_of_foraging+number_of_plants+1:end-1);
        y0_dyn = dB(:,end);

        animal_density_perturb = animal_density_perturb.*(animal_density_perturb>seuil_abondance);
        plant_density_perturb = plant_density_perturb.*(plant_density_perturb>seuil_abondance);

%         animal_total_perturb = [animal_total_perturb ; animal_density_perturb] ;
%         plant_total_perturb = [plant_total_perturb ; plant_density_perturb] ;
%         effort_dynamic_perturb_tot = [effort_dynamic_perturb_tot ; effort_dynamic_perturb] ;
        y0_dyn_tot = [y0_dyn_tot ; y0_dyn] ; 
%         t_post_tot = [t_post_tot ; t + t_post_tot(end)] ;
        % determine if extinction occured (animal)
        test = sum(animal_density_perturb(end,:),'all')>0; % if crash, test=0 and the loop break
    end
y0_dyn_all{i} = y0_dyn_tot; 
end 

% save(['niche_speed_with_fixed_z','.mat'],'y0_dyn_all');
save(['niche_speed_with_fixed_z_101','.mat'],'y0_dyn_all');
