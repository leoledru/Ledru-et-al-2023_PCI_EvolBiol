%% Oceane's code for response 'distance to threshold' with increasing
%% mortality until global extinction 

% ATT: animal_density_perturb et plant_density_perturb donnent la moyenne
% des pas de temps de l'ode à chaque itération de la boucle (à chaque fois
% que la perturb est augmentée) et non pas à chaque pas de temps d'un même
% ODE comme dans les autres.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Increasing mortality to first extinction %%%%%%%

%loop initialization
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

while test_a==1
    
    perturb_increment_a = 0.1*animal_intrinsic_growth ;
    
    % increasing animal mortality to global extinction in the community 
     animal_intrinsic_growth = animal_intrinsic_growth + perturb_increment_a ; %animal_intrinsic_growth = mortality           

     %%%%%%%%%%%%%%%% ODE %%%%%%%%%%%%%%%%%% 
    
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
    animal_density_perturb = [animal_density_perturb ; mean(animal_density_perturb_loop,1)];
    plant_density_perturb = [plant_density_perturb ; mean(plant_density_perturb_loop,1)];
    effort_dynamic_perturb = [effort_dynamic_perturb ; mean(effort_dynamic_perturb_loop,1)] ;
    animal_mortality_dynamic = [animal_mortality_dynamic; animal_intrinsic_growth] ;
    
    % Redefinir conditions initiales
    animal_density_perturb_loop = animal_density_perturb_loop(end,:)';
    plant_density_perturb_loop = plant_density_perturb_loop(end,:)';
    effort_ij0_perturb_loop = effort_dynamic_perturb_loop(end,:)' ;

    %determine if global extinction occured
    test_a = any(animal_density_perturb_loop); 
    
    t_post = [t_post; t+t_post(end)] ;

end  
   
t_post(1) = [] ;
