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
% Kf = @(trait_plant) K0*exp(-(trait_plant - y0).^2./(2*sigmaK^2));
Kf = @(trait_plant) K0*exp(-(trait_plant - y0).^4./(12*sigmaK^4));

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
%% PARFOR
tmax = 1e3;
tspan = 0:10:tmax;
Nt = 10;

% Output_a_f = zeros(loop,Nt,number_of_animals,number_of_foraging);
% Output_p_f = zeros(loop,Nt,number_of_plants);
% Output_a_t = zeros(loop,Nt,number_of_animals,number_of_foraging);
% Output_p_t = zeros(loop,Nt,number_of_plants);
% Effort_f = zeros(loop,Nt,number_of_animals,number_of_plants,number_of_foraging);
% Effort_t = zeros(loop,Nt,number_of_animals,number_of_plants,number_of_foraging);

sigma_loop  = parametres(:,1);
sigmaK_loop = parametres(:,2);
hmax_loop   = parametres(:,3);
animal_intrinsic_growth_loop = parametres(:,4);

Iloop = 20;
for j = 1:loop/Iloop
    idx = (j-1)*20 + (1:20);
    sigma_parfor = sigma_loop(idx);
    sigmaK_parfor = sigmaK_loop(idx);
    hmax_parfor = hmax_loop(idx);
    animal_intrinsic_growth_parfor = animal_intrinsic_growth_loop(idx);

    output_a_f = zeros(Iloop,Nt+1,number_of_animals,number_of_foraging);
    output_p_f = zeros(Iloop,Nt+1,number_of_plants);
    output_a_t = zeros(Iloop,Nt+1,number_of_animals,number_of_foraging);
    output_p_t = zeros(Iloop,Nt+1,number_of_plants);
    effort_f = zeros(Iloop,Nt+1,number_of_animals,number_of_plants,number_of_foraging);
    effort_t = zeros(Iloop,Nt+1,number_of_animals,number_of_plants,number_of_foraging);
    parfor i = 1:Iloop

        %% SIGMA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         SIGMA = sigma_parfor(i);
         delta_ij = complementary_traits(SIGMA,XX,YY);
         delta_ij = delta_ij';

        %% SIGMA_K %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         sigmaK = sigmaK_parfor(i);
         Kf = @(trait_plant) K0*exp(-(trait_plant - y0).^2./(2*sigmaK^2));
         sigmaC = sigmaK-1;
         C = @(trait_plant_y,trait_plant_yi) exp(sigmaC^2*Beta^2/2).*exp(-(trait_plant_y-trait_plant_yi + sigmaC^2*Beta^2).^2./(2*sigmaC^2));
         Compet = C(traits_of_plants,traits_of_plants');

        %% h_max & alpha_h %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hmax = hmax_parfor(i);
        h = @(z) hmin + (hmax - hmin).*z.^alpha_h;
        handling_time = h(foraging_trait);
        handling_time = reshape(handling_time,1,1,number_of_foraging);

        % animal_growth
        animal_intrinsic_growth = animal_intrinsic_growth_parfor(i);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% EFFORT
        effort_ij0 = delta_ij;
        effort_ij0 = effort_ij0./sum(effort_ij0,2);
        effort_ij0 = reshape(effort_ij0,[],1);
        effort_ij0 = repmat(effort_ij0,number_of_foraging,1);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% DIFFUSION MATRIX
        D = 1e-3;
        e = ones(number_of_animals,1);
        A_animal = spdiags([e,-2*e,e],-1:1,number_of_animals,number_of_animals);
        A_animal(1,1) = -1; A_animal(end,end) = -1;
        A_animal = D*A_animal/(dx^2);

        D = 1e-3;
        e = ones(number_of_plants,1);
        A_plant = spdiags([e,-2*e,e],-1:1,number_of_plants,number_of_plants);
        A_plant(1,1) = -1; A_plant(end,end) = -1;
        A_plant = D*A_plant/(dx^2);

        D = 1e-3;
        e = ones(number_of_foraging,1);
        A_foraging = spdiags([e,-2*e,e],-1:1,number_of_foraging,number_of_foraging);
        A_foraging(1,1) = -1; A_foraging(end,end) = -1;
        A_foraging = D*A_foraging/(dx^2);

        %% SIMU AF
        animal_density = zeros(number_of_animals,number_of_foraging);
        plant_density  = zeros(1,number_of_plants);
        plant_density((number_of_plants+1)/2)     = .1;
        animal_density((number_of_animals+1)/2,:) = .1;

        options = odeset('RelTol',1e-3,'AbsTol',1e-5);
        animal_density = reshape(animal_density,number_of_animals*number_of_foraging,1);
        B = [animal_density ; plant_density' ; effort_ij0];


        [t,dB]=ode45(@(t,B) demographic_system_evol_foraging_bis(t, B, number_of_animals, delta_ij, extraction_coeff, conversion_coeff,...
            animal_intrinsic_growth, animal_intraspe_compet, plant_intrinsic_growth, seuil_effort, A_foraging,...
            handling_time, number_of_plants, A_animal, A_plant, dx, dy, dz, Compet, K, seuil_abondance, effort_speed_of_change,...
            number_of_foraging, Foraging_trait, plant_intraspe_compet),tspan, B, options);

        animal_density = dB(end-Nt:end,1:number_of_animals*number_of_foraging);
        plant_density  = dB(end-Nt:end,number_of_animals*number_of_foraging+1:number_of_animals*number_of_foraging+number_of_plants);
        effort_dynamic = dB(end-Nt:end,number_of_animals*number_of_foraging+number_of_plants+1:end);
        
        % densites trop faibles (<seuil) = 0 
        animal_density = animal_density(:,:).*(animal_density(:,:)>seuil_abondance);
        plant_density  = plant_density(:,:).*(plant_density(:,:)>seuil_abondance);

        output_a_f(i,:,:,:) = reshape(animal_density,[Nt+1,number_of_animals,number_of_foraging]);
        output_p_f(i,:,:)   = reshape(plant_density,[Nt+1,number_of_plants]);
        effort_f(i,:,:,:,:) = reshape(effort_dynamic,[Nt+1,number_of_animals,number_of_plants,number_of_foraging]);

        %% SIMU TONDEUSE 
        D = 0;
        e = ones(number_of_foraging,1);
        A_foraging = spdiags([e,-2*e,e],-1:1,number_of_foraging,number_of_foraging);
        A_foraging(1,1) = -1; A_foraging(end,end) = -1;
        A_foraging = D*A_foraging/(dx^2);

        animal_density = zeros(number_of_animals,number_of_foraging);
        plant_density = zeros(1,number_of_plants);
        plant_density((number_of_plants+1)/2)     = .1;
        animal_density((number_of_animals+1)/2,1) = .1;

        animal_density = reshape(animal_density,number_of_animals*number_of_foraging,1);
        B = [animal_density ; plant_density' ; effort_ij0];


        [t,dB]=ode45(@(t,B) demographic_system_evol_foraging_bis(t, B, number_of_animals, delta_ij, extraction_coeff, conversion_coeff,...
            animal_intrinsic_growth, animal_intraspe_compet, plant_intrinsic_growth, seuil_effort, A_foraging,...
            handling_time, number_of_plants, A_animal, A_plant, dx, dy, dz, Compet, K, seuil_abondance, effort_speed_of_change,...
            number_of_foraging, Foraging_trait, plant_intraspe_compet),tspan, B, options);

        animal_density = dB(end-Nt:end,1:number_of_animals*number_of_foraging);
        plant_density  = dB(end-Nt:end,number_of_animals*number_of_foraging+1:number_of_animals*number_of_foraging+number_of_plants);
        effort_dynamic = dB(end-Nt:end,number_of_animals*number_of_foraging+number_of_plants+1:end);

        animal_density = animal_density(:,:).*(animal_density(:,:)>seuil_abondance);
        plant_density  = plant_density(:,:).*(plant_density(:,:)>seuil_abondance);

        output_a_t(i,:,:,:) = reshape(animal_density,[Nt+1,number_of_animals,number_of_foraging]);
        output_p_t(i,:,:)   = reshape(plant_density,[Nt+1,number_of_plants]);
        effort_t(i,:,:,:,:) = reshape(effort_dynamic,[Nt+1,number_of_animals,number_of_plants,number_of_foraging]);

    end
%     Output_a_f(Iloop*(j-1)+1:Iloop*j,:,:,:) = output_a_f;
%     Output_p_f(Iloop*(j-1)+1:Iloop*j,:,:,:) = output_p_f;
%     Output_a_t(Iloop*(j-1)+1:Iloop*j,:,:,:) = output_a_t;
%     Output_p_t(Iloop*(j-1)+1:Iloop*j,:,:,:) = output_p_t;
    Param = parametres(1:Iloop*j,:);
%     save(['community_comparison_b','.mat'],'Output_a_f','Output_p_f','Effort_f'...
%     ,'Output_a_t','Output_p_t','Effort_t','Param');
    save(['community_comparison_K4_',num2str(j),'.mat'],'output_a_f','output_p_f','effort_f'...
    ,'output_a_t','output_p_t','effort_t','Param');

end
% Param = parametres';
% % save(['community_comparison_b','.mat'],'Output_a_f','Output_p_f','Effort_f'...
% %     ,'Output_a_t','Output_p_t','Effort_t','Param');
% save(['community_comparison_K4','.mat'],'Output_a_f','Output_p_f','Effort_f'...
%     ,'Output_a_t','Output_p_t','Effort_t','Param');

