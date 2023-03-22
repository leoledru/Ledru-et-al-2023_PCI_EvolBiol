%% CODE DE LA PERTURBATION NICHE-SPEED
%% L'OBJET PREP.MAT EST APPELE CAR C'EST LA STRUCTURE DANS LAQUELLE TOUS
%% LES RESULTATS SONT ENREGISTRES
clear all;close all

options = odeset('RelTol',1e-3,'AbsTol',1e-5);
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

%% Initial conditions
load('../Data/prep.mat')
tmax = 2e3;
tspan = 0:20:tmax;
nt = length(tspan);
% DENSITIES
animal_density_pre_perturb_f = prep{1,2}.niche_speed.animal_density(end,:);
plant_density_pre_perturb_f  = prep{1,2}.niche_speed.plant_density(end,:);
animal_density_pre_perturb_t = prep{1,1}.niche_speed.animal_density(end,:);
plant_density_pre_perturb_t = prep{1,1}.niche_speed.plant_density(end,:);

% EFFORT
effort_ij0_pre_perturb_f = prep{1,2}.niche_speed.effort_dynamic(end,:);
effort_ij0_pre_perturb_t = prep{1,1}.niche_speed.effort_dynamic(end,:);

clear prep
%% Speed of change
test = 1;
number_c = 1;%0;
increment_loop = linspace(0,.15,number_c+1);
prep = {};

for jk = 1:2
    animal_density_perturb = zeros(number_c*length(tspan),number_of_animals*number_of_foraging);
    plant_density_perturb =  zeros(number_c*length(tspan),number_of_plants);
    effort_dynamic_perturb = zeros(number_c*length(tspan),number_of_plants*number_of_animals*number_of_foraging);
    t_post = zeros(20*length(tspan),1) ;
    
    ji = 1;
    
    while (test == 1)&&(ji<number_c+1)
        
        increment = increment_loop(ji+1);
        
        %%%%%%%%%%%%%%%% ODE %%%%%%%%%%%%%%%%%%
        %%% Initial conditions
        if (jk==2)
            animal_density = animal_density_pre_perturb_f;
            plant_density = plant_density_pre_perturb_f;
            effort_ij0 = effort_ij0_pre_perturb_f;
            %%% Diffusion A_foraging
            D = 1e-3; % Diffusion coefficient
            %         D = 0;
            e = ones(number_of_foraging,1);
            A_foraging = spdiags([e,-2*e,e],-1:1,number_of_foraging,number_of_foraging);
            A_foraging(1,1) = -1; A_foraging(end,end) = -1;
            A_foraging = D*A_foraging/(dx^2);
        elseif(jk==1)
            animal_density = animal_density_pre_perturb_t;
            plant_density = plant_density_pre_perturb_t;
            effort_ij0 = effort_ij0_pre_perturb_t;
            %%% Diffusion A_foraging
            % D = 1e-3; % Diffusion coefficient
            D = 0;
            e = ones(number_of_foraging,1);
            A_foraging = spdiags([e,-2*e,e],-1:1,number_of_foraging,number_of_foraging);
            A_foraging(1,1) = -1; A_foraging(end,end) = -1;
            A_foraging = D*A_foraging/(dx^2);
        end
        B0 = [animal_density' ; plant_density' ; effort_ij0'; y0'];
        
        %%% Solving ODE
        [t, dB] = ode45(@(t,B) demographic_system_evol_foraging_bis_perturb(t, B, number_of_animals, delta_ij, extraction_coeff, conversion_coeff,...
            animal_intrinsic_growth, animal_intraspe_compet, plant_intrinsic_growth, A_foraging,...
            handling_time, number_of_plants, A_animal, A_plant, dx, dy, dz, Compet, seuil_abondance,...
            effort_speed_of_change, number_of_foraging, Foraging_trait, traits_of_plants_b, traits_of_plants, K0, sigmaK, increment), tspan, B0, options) ;
        
        animal_density = dB(:,1:number_of_animals*number_of_foraging);
        plant_density = dB(:,number_of_animals*number_of_foraging+1:number_of_animals*number_of_foraging+number_of_plants);
        effort_dynamic = dB(:,number_of_animals*number_of_foraging+number_of_plants+1:end-1);
        y0_dyn = dB(:,end);
        
        animal_density = animal_density.*(animal_density>seuil_abondance);
        plant_density = plant_density.*(plant_density>seuil_abondance);
        
        animal_density_perturb(nt*(ji-1)+1:nt*ji,:) = animal_density;
        plant_density_perturb(nt*(ji-1)+1:nt*ji,:)  = plant_density;
        effort_dynamic_perturb(nt*(ji-1)+1:nt*ji,:) = effort_dynamic;
        
        t_post(nt*(ji-1)+1:nt*ji,1) = t ;
        
        % determine if extinction occured (animal)
        test = sum(animal_density(end,:),'all')>0; % if crash, test=0 and the loop break
        ji = ji+1;
    end
    
    prep{jk}.niche_speed.animal_density_perturb = animal_density_perturb ;
    prep{jk}.niche_speed.plant_density_perturb = plant_density_perturb ;
    prep{jk}.niche_speed.effort_dynamic_perturb = effort_dynamic_perturb ;
    prep{jk}.niche_speed.t_post = t_post ;
    
    prep{jk}.niche_speed.number_of_animals = number_of_animals ;
    prep{jk}.niche_speed.number_of_plants = number_of_plants;
    prep{jk}.niche_speed.number_of_foraging = number_of_foraging ;
end

