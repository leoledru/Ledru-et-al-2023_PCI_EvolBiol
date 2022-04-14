%% Oceane's code for niche-jump perturbation
%March 2020

%ATT en preperturb, 1 ligne = 1 pas de temps, mais en post perturb = 1 loop
%(à chaque augmentation de deplacement de niche)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('prep.mat') 

%% % INITIALISATION %%%

number_of_animals = 77 ; %
number_of_plants  = 77 ;%
number_of_foraging = 11 ;%


% FONCTION K(y)
K0 = 50; % maximal carrying capacity
y0 = -25; % optimal trait pour deplacemnt niche %
sigmaK = 2.5;
Kf = @(trait_plant) K0*exp(-(trait_plant - y0).^2./(2*sigmaK^2));

% FONCTION C(y-y0)
Beta = 0; % Beta = 0 : symetrical competition
sigmaC = sigmaK-1;
C = @(trait_plant_y,trait_plant_yi) exp(sigmaC^2*Beta^2/2).*exp(-(trait_plant_y-trait_plant_yi + sigmaC^2*Beta^2).^2./(2*sigmaC^2));


% TRAITS %
% Animal traits

xmin = -35; %pour deplacement niche
xmax = 35;
xx = linspace(xmin,xmax,number_of_animals);
%dx = xx(2)-xx(1);
dx = 1;
traits_of_animals = xx';

foraging_trait = linspace(0,1,number_of_foraging);
%dz = foraging_trait(2) - foraging_trait(1);
dz = 1;
Foraging_trait = reshape(foraging_trait,1,1,number_of_foraging);

% Plants traits
ymin = -35; % pour deplacement niche
ymax = 35;
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
hmax = 0.55;
h = @(z) hmin + (hmax - hmin).*z.^alpha_h;
handling_time = h(foraging_trait);
handling_time = reshape(handling_time,1,1,number_of_foraging);



%% parametres pour lancer commu tondeuse et af

N = [0,1e-3] ; 

for aj = 1:length(prep)

% FONCTION K(y)
Kf = @(trait_plant) K0*exp(-(trait_plant - y0).^2./(2*sigmaK^2));
% Carrying capacity %
K = Kf(traits_of_plants);
    

% DENSITIES

% animal_density = (xx<3*dx).*(xx>-3*dx);
% plant_density  = (yy<3*dy).*(yy>-3*dy);
animal_density = zeros(number_of_animals,number_of_foraging);
plant_density = zeros(1,number_of_plants);
% animal_density((number_of_animals+1)/2,1) = .1; %one on niche's center
% plant_density((number_of_plants+1)/2) = .1;

animal_density(10,1) = .1; %pour dep niche
plant_density(10) = .1;

% animal_density(number_of_animals/2,1) = .1; %one on niche's center pour impair
% plant_density(number_of_plants/2) = .1;

%   animal_density(1,1) = .1; % one "tondeuse"

% animal_density(:,:) = 1; % animal everywhere
% plant_density(:) = 1; % plant everywhere
%animal_density(15,:) = 1;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
 
%% 

  prep{aj}.Df = N(aj) ; %commu foraging ou tondeuse
        
e = ones(number_of_foraging,1);
A_foraging = spdiags([e,-2*e,e],-1:1,number_of_foraging,number_of_foraging);
A_foraging(1,1) = -1; A_foraging(end,end) = -1;
%         A_foraging(1,2) = -1; A_foraging(end,end-1) = -1;
A_foraging = prep{aj}.Df*A_foraging/(dx^2) ;

%

%%%%%%%%%%%%%%%%%%%%%%%% ODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Effort_ij0 = reshape(effort_ij,[],1);

animal_density = reshape(animal_density,number_of_animals*number_of_foraging,1);

B = [animal_density ; plant_density' ; effort_ij0];

tmax = 5000 ;
tspan = 0 : 50 : tmax ;
options = odeset('RelTol',1e-3,'AbsTol',1e-4);

[t,dB]=ode45(@(t,B) demographic_system_evol_foraging_bis(t, B, number_of_animals, delta_ij, extraction_coeff, conversion_coeff,...
animal_intrinsic_growth, animal_intraspe_compet, plant_intrinsic_growth, seuil_effort, A_foraging,...
handling_time, number_of_plants, A_animal, A_plant, dx, dy, dz, Compet, K, seuil_abondance, effort_speed_of_change,...
number_of_foraging, Foraging_trait, plant_intraspe_compet),tspan, B, options);


animal_density = dB(:,1:number_of_animals*number_of_foraging);
plant_density = dB(:,number_of_animals*number_of_foraging+1:number_of_animals*number_of_foraging+number_of_plants);
effort_dynamic = dB(:,number_of_animals*number_of_foraging+number_of_plants+1:end);

prep{aj}.niche_displacement_range.effort_dynamic = effort_dynamic ; 

animal_density = animal_density.*(animal_density>seuil_abondance);
prep{aj}.niche_displacement_range.animal_density = animal_density ;

plant_density = plant_density.*(plant_density>seuil_abondance);
prep{aj}.niche_displacement_range.plant_density = plant_density ;


t_pre = t;
prep{aj}.niche_displacement_range.t_pre = t ;

%ISOLATE LAST TIME LAPS to initialize perturbation 
 animal_density_last_step = animal_density(end,:) ;
 animal_density_last_step = animal_density_last_step' ;
 plant_density_last_step = plant_density(end,:) ;
 plant_density_last_step = plant_density_last_step' ;

 effort_ij0_last_step = effort_dynamic(end,:);
 effort_ij0_last_step = effort_ij0_last_step';



 %% POST-PERTURB
 
test = 1 ;
y0_perturb = y0 ;
perturb_increment = 1;

animal_density_perturb = animal_density_last_step;
plant_density_perturb = plant_density_last_step;
effort_ij0_perturb = effort_ij0_last_step;


animal_total_perturb = [];
plant_total_perturb = [];
effort_dyn_tot_perturb = [] ;
y0_dyn = [] ;
T_post = 0 ;


    while test==1 
    
        animal_density_perturb = animal_density_last_step;
        plant_density_perturb = plant_density_last_step;

        %Increasing 'amplitude' of y
        y0_perturb = y0_perturb + perturb_increment ; % optimal trait 
        Kf = @(trait_plant) K0*exp(-(trait_plant - y0_perturb).^2./(2*sigmaK^2));
        % Carrying capacity %
        K = Kf(traits_of_plants);

        %%%%%%%%%%%%%%%% ODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        B = [animal_density_perturb ; plant_density_perturb ; effort_ij0_perturb]; %plant density retourné avant
        
        tmax = 3000 ;
        tspan = 0:30:tmax ;
    
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
        T_post = [T_post; t + T_post(end)] ;

        %determine if extinction of the whole animal community occured
        test = sum(animal_density_perturb,'all')>0; 
    
    end

T_post(1) = [] ;

%save results
 prep{aj}.niche_displacement_range.animal_density_perturb =  animal_total_perturb ;  %attention ne correspond pas à des pas de temps mais à chaque baisse de mortalité
 prep{aj}.niche_displacement_range.plant_density_perturb = plant_total_perturb ;
 prep{aj}.niche_displacement_range.effort_dynamic_perturb = effort_dyn_tot_perturb ; 
 prep{aj}.niche_displacement_range.y0_dyn = y0_dyn ;
 prep{aj}.niche_displacement_range.t_post = t ;
 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PROXIS AND METRICS %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

name_prox = {'biomass_post_p', 'biomass_post_a', 'biomass_post_tot', 'biomass_pre_p', 'biomass_pre_a', 'biomass_pre_tot' ...
              'productivity_post', 'productivity_pre', ...
              'FDis_post_p','FDis_post_a', 'FDis_pre_p', 'FDis_pre_a', ...
              'FRO_pre_p', 'FRO_pre_a', 'FRO_post_p', 'FRO_post_a'} ;
    
    
          % loading var %utile car chgmnt de nom des var en haut (animal total & animal perturb)
        animal_density = prep{aj}.niche_displacement_range.animal_density ;
        plant_density = prep{aj}.niche_displacement_range.plant_density ;
        effort_dynamic = prep{aj}.niche_displacement_range.effort_dynamic ;
        animal_density_perturb = prep{aj}.niche_displacement_range.animal_density_perturb ;
        plant_density_perturb = prep{aj}.niche_displacement_range.plant_density_perturb ;
        effort_dynamic_perturb = prep{aj}.niche_displacement_range.effort_dynamic_perturb ;
        t_pre = prep{aj}.niche_displacement_range.t_pre ;        
        t_post = prep{aj}.niche_displacement_range.t_post ;
        
        % PROXIS %
        proxis_process
        
        % save results of proxis_process
        prep{aj}.niche_displacement_range.proxis = proxis ;
      
        for ji = 1:length(name_prox)
            
            % mean and variance of proxis
            prep{aj}.niche_displacement_range.proxis.([sprintf('%s', name_prox{ji}), '_', 'mean']) = mean(prep{aj}.niche_displacement_range.proxis.(name_prox{ji})) ;
            prep{aj}.niche_displacement_range.proxis.([sprintf('%s', name_prox{ji}), '_', 'var']) = var(prep{aj}.niche_displacement_range.proxis.(name_prox{ji})) ;
        end

        % METRICS  
        Tolerance = tolerance(biomass_pre_p,biomass_post_p, biomass_pre_a, biomass_post_a, biomass_pre_tot, ...
            biomass_post_tot, productivity_pre, productivity_post, productivity_pond_pre, ...
            productivity_pond_post, FDis_pre_p, FDis_post_p, FDis_pre_a, FDis_post_a, ...
            FRO_pre_p, FRO_post_p, FRO_pre_a, FRO_post_a, t_pre, t_post);    

        %save results
        prep{aj}.niche_displacement_range.tolerance = Tolerance ;    
        
end

save niche.mat prep