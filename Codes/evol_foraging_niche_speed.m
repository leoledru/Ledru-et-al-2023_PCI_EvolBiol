%% CODE DE LA PERTURBATION NICHE-SPEED
%% L'OBJET PREP.MAT EST APPELE CAR C'EST LA STRUCTURE DANS LAQUELLE TOUS
%% LES RESULTATS SONT ENREGISTRES

load('prep.mat')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%INITIALISATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%

number_of_animals = 77;
number_of_plants  = 77;
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
xmin = -35;
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
ymin = -35;
ymax = 35;
yy = linspace(ymin,ymax,number_of_plants);
%dy =  yy(2)-yy(1);
dy = 1;
traits_of_plants = yy';
traits_of_plants_b = linspace(-70,70,77*2-1)';

[YY,XX] = meshgrid(yy,xx);
SIGMA = .5;
delta_ij = complementary_traits(SIGMA,XX,YY);
delta_ij = delta_ij';

% effort_speed_of_change = .5;
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

% Set the niche center and shift K according to the new center
% y0 = 35;
% K = circshift(K,y0);


%% parametres pour lancer commu tondeuse et af

N = [0,1e-3] ; 

for jk = 1:length(prep)


    % DENSITIES

    % animal_density = (xx<3*dx).*(xx>-3*dx);
    % plant_density  = (yy<3*dy).*(yy>-3*dy);
    animal_density = zeros(number_of_animals,number_of_foraging);
    plant_density = zeros(1,number_of_plants);
    animal_density((number_of_animals+1)/2,1) = .1; %one on niche's center
    plant_density((number_of_plants+1)/2) = .1;
    % animal_density(1,1) = .1;
    % plant_density(1) = .1;

    % y0 = traits_of_plants(5);
    % Kf = @(trait_plant) K0*exp(-(trait_plant - y0).^2./(2*sigmaK^2));
    % K = Kf(traits_of_plants);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % EFFORT

    effort_ij0 = delta_ij;
    effort_ij0 = effort_ij0./sum(effort_ij0,2);
    effort_ij0 = reshape(effort_ij0,[],1);
    effort_ij0 = repmat(effort_ij0,number_of_foraging,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DIFFUSION MATRIX

    D = 1e-3; % Diffusion coefficient
    %         D = 0;
    e = ones(number_of_animals,1);
    A_animal = spdiags([e,-2*e,e],-1:1,number_of_animals,number_of_animals);
    % A_animal(1,1) = -1; A_animal(end,end) = -1; % no wrapping
    A_animal(1,end) = 1; A_animal(end,1) = 1; % wrapping niche mutation
    A_animal = D*A_animal/(dx^2);

    D = 1e-3; % Diffusion coefficient
    %         D = 0;
    e = ones(number_of_plants,1);
    A_plant = spdiags([e,-2*e,e],-1:1,number_of_plants,number_of_plants);
    % A_plant(1,1) = -1; A_plant(end,end) = -1; % no wrapping
    A_plant(1,end) = 1; A_plant(end,1) = 1; % wrapping niche mutation
    A_plant = D*A_plant/(dx^2);

    % Soit commu foraging soit tondeuse
    prep{jk}.Df = N(jk) ;
        
    e = ones(number_of_foraging,1);
    A_foraging = spdiags([e,-2*e,e],-1:1,number_of_foraging,number_of_foraging);
    A_foraging(1,1) = -1; A_foraging(end,end) = -1;
    %         A_foraging(1,2) = -1; A_foraging(end,end-1) = -1;
    A_foraging = prep{jk}.Df*A_foraging/(dx^2) ;

%%%%%%%%%%%%%%%%%%%%%%%% ODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
                    
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
prep{jk}.niche_speed.effort_dynamic = effort_dynamic ;

y0_dyn = dB(:,end);

animal_density = animal_density.*(animal_density>seuil_abondance);
prep{jk}.niche_speed.animal_density = animal_density ;
plant_density = plant_density.*(plant_density>seuil_abondance);
prep{jk}.niche_speed.plant_density = plant_density ;

t_pre = t ;
prep{jk}.niche_speed.t_pre = t ;

%Isolate last time laps to initialize perturbation 
%%
% DENSITIES
  animal_density_last_step = animal_density(end,:) ;
  animal_density_last_step = animal_density_last_step' ;
  plant_density_last_step = plant_density(end,:) ;
  plant_density_last_step = plant_density_last_step' ;

% EFFORT
  effort_ij0_last_step = effort_dynamic(end,:);
  effort_ij0_last_step = effort_ij0_last_step';

%%

%loop initialization
test = 1;
increment_loop = linspace(.05,.2,20);

animal_total_perturb = [];
plant_total_perturb = [];
y0_dyn_tot = [];
effort_dynamic_perturb_tot = [];
t_post_tot = 0 ;

ji = 0;

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

      %  options = odeset('RelTol',1e-3,'AbsTol',1e-5,'Events',@myEvent);
        
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

        animal_total_perturb = [animal_total_perturb ; animal_density_perturb] ;
        plant_total_perturb = [plant_total_perturb ; plant_density_perturb] ;
        effort_dynamic_perturb_tot = [effort_dynamic_perturb_tot ; effort_dynamic_perturb] ;
        y0_dyn_tot = [y0_dyn_tot ; y0_dyn] ; 
        t_post_tot = [t_post_tot ; t + t_post_tot(end)] ;

        % determine if extinction occured (animal)
        test = sum(animal_density_perturb(end,:),'all')>0; % if crash, test=0 and the loop break

    end

    t_post_tot(1) =[] ; %enlever 1e ligne (0)
  %  t_post = t ;
 
 %save results
 prep{jk}.niche_speed.animal_density_perturb = animal_total_perturb ;
 prep{jk}.niche_speed.plant_density_perturb = plant_total_perturb ;
 prep{jk}.niche_speed.effort_dynamic_perturb = effort_dynamic_perturb_tot ; 
 prep{jk}.niche_speed.y0_dyn = y0_dyn_tot ;
 prep{jk}.niche_speed.t_post = t_post_tot ;
 
 prep{jk}.niche_speed.number_of_animals = number_of_animals ;
 prep{jk}.niche_speed.number_of_plants = number_of_plants;
 prep{jk}.niche_speed.number_of_foraging = number_of_foraging ;
  
end 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% PROXIS AND METRICS %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
name_prox = {'biomass_post_p', 'biomass_post_a', 'biomass_post_tot', 'biomass_pre_p', 'biomass_pre_a', 'biomass_pre_tot' ...
              'productivity_post', 'productivity_pre', ...
              'FDis_post_p','FDis_post_a', 'FDis_pre_p', 'FDis_pre_a', ...
              'FRO_pre_p', 'FRO_pre_a', 'FRO_post_p', 'FRO_post_a' } ;
    
for jk = 1:length(prep)
    
    % PROXIS %
      % loading var %utile car chgmnt de nom des var en haut (animal total & animal perturb
        animal_density = prep{jk}.niche_speed.animal_density ;
        plant_density = prep{jk}.niche_speed.plant_density ;
        effort_dynamic = prep{jk}.niche_speed.effort_dynamic ;
        animal_density_perturb = prep{jk}.niche_speed.animal_density_perturb ;
        plant_density_perturb = prep{jk}.niche_speed.plant_density_perturb ;
        effort_dynamic_perturb = prep{jk}.niche_speed.effort_dynamic_perturb ;
        t_pre = prep{jk}.niche_speed.t_pre ;        
        t_post = prep{jk}.niche_speed.t_post ;
    
    % compute proxis
    proxis_process
        
    % save results of proxis_process
    prep{jk}.niche_speed.proxis = proxis ;        
       
    
   for a = 1:length(name_prox)
        
        % mean and variance of proxis
        prep{jk}.niche_speed.proxis.([sprintf('%s', name_prox{a}), '_', 'mean']) = mean(prep{jk}.niche_speed.proxis.(name_prox{a})) ;
        prep{jk}.niche_speed.proxis.([sprintf('%s', name_prox{a}), '_', 'var']) = var(prep{jk}.niche_speed.proxis.(name_prox{a})) ;
   end

    % METRICS 
    %tolerance
      Tolerance = tolerance(biomass_pre_p,biomass_post_p, biomass_pre_a, biomass_post_a, biomass_pre_tot, ...
    biomass_post_tot, productivity_pre, productivity_post, productivity_pond_pre, ...
    productivity_pond_post, FDis_pre_p, FDis_post_p, FDis_pre_a, FDis_post_a, ...
    FRO_pre_p, FRO_post_p, FRO_pre_a, FRO_post_a, t_pre , t_post) ;
    
    %save results
    prep{jk}.niche_speed.tolerance = Tolerance ; 
     
end    
    
save prep_speed.mat prep
