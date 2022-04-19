%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Océane's code of pulse perturbation on biomass %%

% RQ: on proportionnalise par 10 fois la biomasse %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plant_density_last_step = prep{1,2}.plant_density_last_step ;
animal_density_last_step = prep{1,2}.animal_density_last_step ;



% Biomass of the last step %
% biomasse_totale = sum(animal_density_last_step) + sum(plant_density_last_step) ;
biomass_last_p = sum(plant_density_last_step) ;
biomass_last_a = sum(animal_density_last_step) ;

% tirer vecteur de taile des densités %
% v_pert_tot = rand(1 , length(animal_density_last_step) + length(plant_density_last_step)) ;
v_pert_p = rand(1 , length(plant_density_last_step)) ;
v_pert_a = rand(1 , length(animal_density_last_step)) ;


% normalize on 1 %
% v_pert_tot = v_pert_tot./sum(v_pert_tot) ;
v_pert_p = v_pert_p./sum(v_pert_p) ;
v_pert_a = v_pert_a./sum(v_pert_a) ;

% proportionnaliser par la biomasse % ATT Ici par 3* la biomasse 
% v_pert_tot = biomasse_totale.*v_pert_tot ; %sum = biomasse tot
v_pert_p = 3*biomass_last_p.*v_pert_p ;
v_pert_a = 3*biomass_last_a.*v_pert_a ;

% tirer vecteur de meme taille avec 1 ou -1 %
    %vec_tot = zeros(1,length(v_pert_tot)) ;
    %vec_tot([1:length(vec_tot)/2]) = 1 ; vec_tot([length(vec_tot)/2:end]) = -1 ;
    %idx = randperm(length(vec_tot)) ;
    %vec_tot = vec_tot((idx)) ;

vec_p = zeros(1,length(v_pert_p)) ;
vec_p(1:round((length(vec_p))/2)) = 1 ; vec_p(round((length(vec_p))/2):end) = -1 ; %nb d'ind impair, kékonfé? un -1 de plus
vec_p = vec_p(randperm(length(vec_p)));


vec_a = zeros(1,length(v_pert_a)) ; 
vec_a(1:round(length(vec_a)/2)) = 1 ; vec_a((round(length(vec_a)+1)/2):end) = -1 ;
vec_a = vec_a(randperm(length(vec_a)));

%  vecteur densite .* ce vecteur %
    %v_pert_tot = v_pert_tot.*vec_tot ; 
v_pert_p = v_pert_p.*vec_p ;
v_pert_a = v_pert_a.*vec_a ; 

%chercher les densités nulles
null_density_p = plant_density_last_step == 0 ;
null_density_a = animal_density_last_step == 0 ;

%%
% somme vecteur des densités + vec1 %
    % total_density_perturb = v_pert_tot + plant_density_last_step ; 
plant_density_last_step = v_pert_p' + plant_density_last_step ; 
%animal_density_last_step = v_pert_a' + animal_density_last_step ;

% enlever les perturb sur les densités nulles (pour éviter la création d'espèces)
plant_density_last_step(null_density_p) = 0 ;
%animal_density_last_step(null_density_a) = 0 ;

% ramener les valeurs neg à 0
plant_density_last_step(plant_density_last_step<0)=0 ;
animal_density_last_step(animal_density_last_step<0)=0 ;


%% %%%%%%%%%%%%%%%%%%%% ODE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B = [animal_density_last_step ; plant_density_last_step ; effort_ij0_last_step];
    
options = odeset('RelTol',1e-3,'AbsTol',1e-4);
% tmax = 10000 ;
%tspan = [0,tmax] ;  
        
[t,dB]=ode45(@(t,B) demographic_system_evol_foraging_bis(t, B, number_of_animals, delta_ij, extraction_coeff, conversion_coeff,...
animal_intrinsic_growth, animal_intraspe_compet, plant_intrinsic_growth, seuil_effort, A_foraging,...
handling_time, number_of_plants, A_animal, A_plant, dx, dy, dz, Compet, K, seuil_abondance, effort_speed_of_change,...
number_of_foraging, Foraging_trait, plant_intraspe_compet),tspan, B, options);
      
animal_density_perturb = dB(:,1:number_of_animals*number_of_foraging);
plant_density_perturb = dB(:,number_of_animals*number_of_foraging+1:number_of_animals*number_of_foraging+number_of_plants);
effort_dynamic_perturb = dB(:,number_of_animals*number_of_foraging+number_of_plants+1:end);
    
animal_density_perturb = animal_density_perturb.*(animal_density_perturb>seuil_abondance);
plant_density_perturb = plant_density_perturb.*(plant_density_perturb>seuil_abondance);

t_post = t ;
