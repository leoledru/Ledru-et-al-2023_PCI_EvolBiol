%% CODE CONTENANT TOUTES LES PERTURBATIONS EN COMMENÇANT PAR APPELER 
%% LE RESULTAT DE L'ÉMERGENCE (initialisation.mat).
%% SAVE TOUS LES RÉSULTATS DANS LA STRUCTURE PREP.MAT

% %%%%%%%%%%%%%%%%%%%%%%%%%
%% INITIALISATION %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Paramètres communs pour pre
% tmax = 10000 ; % +haut que post_p
% tspan = 0:75:tmax ;
% options = odeset('RelTol',1e-6,'AbsTol',1e-8);
% 
% % tspan = [0,tmax] ;
% 
% %%%
% initialisation_pre_perturb_process


%% %%%%%%%%%%%%%%%%%%%%%%%%%
%%% PERTURBATIONS %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% %%%%%%%%%%%%%%%%%%%%%%%% PULSE PERTURBATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
load initialisation.mat
load prep.mat

% Parametres pour ODE post
tmax = 3000 ;
tspan = 0 : 10 : tmax ;


%%% Pulse perturb on plants, animals or all
name = {'pulse_on_plant' , 'pulse_on_animal' , 'pulse_on_all'} ;     
name_metric = {'resilience_p_B', 'resilience_a_B', 'resilience_tot_B',...
                'resilience_P', 'resilience_p_FD', 'resilience_a_FD',...
                'resilience_p_FRO', 'resilience_a_FRO'} ;


            
            
%% Pour plusieurs simus

% Pour initialiser les 10 simus
for cc = 1:length(name_metric)

    prep{1}.pulse_on_plant.resilience_all_simus.(name_metric{cc}) = [] ; 
    prep{2}.pulse_on_plant.resilience_all_simus.(name_metric{cc}) = [] ;

    prep{1}.pulse_on_animal.resilience_all_simus.(name_metric{cc}) = [] ;
    prep{2}.pulse_on_animal.resilience_all_simus.(name_metric{cc}) = [] ;

    prep{1}.pulse_on_all.resilience_all_simus.(name_metric{cc}) = [] ; 
    prep{2}.pulse_on_all.resilience_all_simus.(name_metric{cc}) = [] ;

end


for aa = 1:1000 % nb simus
        
        %code perturb    
        pulse_perturb_process
        
    for bc = 1 %type de pulse  
        for cc = 1:length(name_metric) 

            %stocker
             prep{1}.(name{bc}).resilience_all_simus.(name_metric{cc}) = [prep{1}.(name{bc}).resilience_all_simus.(name_metric{cc}); prep{1}.(name{bc}).resilience.(name_metric{cc})] ; %tondeuse
             prep{2}.(name{bc}).resilience_all_simus.(name_metric{cc}) = [prep{2}.(name{bc}).resilience_all_simus.(name_metric{cc}); prep{2}.(name{bc}).resilience.(name_metric{cc})] ; %af
        end            
    end
    
    pulse_on_plant = {prep{1}.pulse_on_plant, prep{2}.pulse_on_plant} ;
%     save pulse_en_cours pulse_on_plant    
end


% save results_100_pulse_on_plants.mat pulse_on_plant -v7.3


%% Compute mean and variance 
% remplacer les never leave equilibrum & not resilient
tmax = 3000 ; % a verifier si nouveau code

for j = 1:2 %tondeuse ou af
    
    resilience_all_simus_num = struct2cell(prep{j}.pulse_on_plant.resilience_all_simus);
    names = fieldnames(prep{j}.pulse_on_plant.resilience_all_simus) ;

    for i = 1:length(resilience_all_simus_num) %mm taille que af

%         plot_p_t = prep{1}.pulse_on_plant.resilience_all_simus.resilience_p_B ;     
%         plot_p_af = prep{2}.pulse_on_plant.resilience_all_simus.resilience_p_B ;     
% 
%         plot_a_t = prep{1}.pulse_on_plant.resilience_all_simus.resilience_a_B ;     
%         plot_a_af = prep{2}.pulse_on_plant.resilience_all_simus.resilience_a_B ;     


        if isstring(resilience_all_simus_num{i}) ==1
            
        %   Remplacer never et not resilient par 0 et tmax        
%             resilience_all_simus_num{i}(resilience_all_simus_num{i} == "never leave equilibrum") = 0 ;            
%             resilience_all_simus_num{i}(resilience_all_simus_num{i} == "not resilient") = tmax ;        
%             resilience_all_simus_num{i} = str2double(resilience_all_simus_num{i}) ;
    
        % ENLEVER never et not resilient du calcul des moyennes
            resilience_all_simus_num{i} = str2double(resilience_all_simus_num{i}) ;
            resilience_all_simus_num{i} = rmmissing(resilience_all_simus_num{i}) ;

        end
      
    end 

    
    prep{j}.pulse_on_plant.resiliences_nostring = cell2struct(resilience_all_simus_num, names) ; % donne resiliences sans str
    prep{j}.pulse_on_plant.resiliences_mean = structfun(@mean,prep{j}.pulse_on_plant.resiliences_nostring,'UniformOutput',false ) ; 
    prep{j}.pulse_on_plant.resiliences_var = structfun(@var,prep{j}.pulse_on_plant.resiliences_nostring,'UniformOutput',false ) ; 
    prep{j}.pulse_on_plant.resiliences_mode = structfun(@mode, prep{j}.pulse_on_plant.resiliences_nostring,'UniformOutput',false) ;
    %+ mode??? 
   
end



 % %%%%%%%%%%%%%%%%%%%% PRESS PERTURBATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Press perturbation on sigma_K %%%%%%%%%%%%%%%%

load initialisation.mat
load prep.mat
tmax = 3000 ;
tspan = 0 : 20 : tmax ;

for vv = 1:length(prep) %Tondeuse et AF


% load('pre_perturb.mat') ; % reecrire les parametres ecrasés par les autres perturbs

% enregistrer var nécessaires pour les perturb
animal_density_last_step = prep{vv}.pre_perturb.animal_density_last_step ;
plant_density_last_step = prep{vv}.pre_perturb.plant_density_last_step ;
effort_ij0_last_step = prep{vv}.pre_perturb.effort_ij0_last_step ; 
A_foraging = prep{vv}.pre_perturb.A_foraging ;

% Perturb
sigmaK_perturb_code 

% Save results
prep{vv}.press_sigmaK_decrease = press_sigmaK{1,1} ;
prep{vv}.press_sigmaK_increase = press_sigmaK{1,2} ; 

end

save prep.mat prep -v7.3

%

%% PRESS PERTURBATION on mortality %%%%%%%
% Increasing mortality one time at 10% of intrinsic mortality 


%enregistrer var nécessaires pour les perturb
load initialisation.mat
load prep.mat
tmax = 3000 ;
tspan = 0 : 20 : tmax ;

%load('pre_perturb.mat')


for vv = 1:length(prep) %Tondeuse et AF
    
animal_density_last_step = prep{vv}.pre_perturb.animal_density_last_step ;
plant_density_last_step = prep{vv}.pre_perturb.plant_density_last_step ;
effort_ij0_last_step = prep{vv}.pre_perturb.effort_ij0_last_step ; 
A_foraging = prep{vv}.pre_perturb.A_foraging ;

%perturb
increased_mortality_10percent 

%save results
prep{vv}.press_mortality_decrease = press_mortality{1} ;
prep{vv}.press_mortality_increase = press_mortality{2} ; 

end

save prep.mat prep -v7.3

%%  PRESS PERTURBATIONS "SEUIL" %%%%%%%%%%%%%%%%%%%
% Increasing mortality to extinction %%%

load initialisation.mat
% load prep.mat
% load('pre_perturb.mat')
tmax = 3000 ;
tspan = 0:20:tmax ;
options = odeset('RelTol',1e-6,'AbsTol',1e-8); %ATT plus petit


for vv = 1:length(prep)
    
%enregistrer var nécessaires pour les perturb

animal_density_last_step = prep{vv}.pre_perturb.animal_density_last_step ;
plant_density_last_step = prep{vv}.pre_perturb.plant_density_last_step ;
effort_ij0_last_step = prep{vv}.pre_perturb.effort_ij0_last_step ;  
A_foraging = prep{vv}.pre_perturb.A_foraging ;

%Perturb
increasing_mortality_threshold

%Save results
prep{vv}.press_mortality_threshold.animal_density_perturb = animal_density_perturb ;
prep{vv}.press_mortality_threshold.plant_density_perturb = plant_density_perturb ;
prep{vv}.press_mortality_threshold.effort_dynamic_perturb = effort_dynamic_perturb ;
prep{vv}.press_mortality_threshold.t_post = t_post ;
prep{vv}.press_mortality_threshold.animal_mortality_dynamic = animal_mortality_dynamic ;

end

save prep.mat prep -v7.3

%% %%%%%%%%%%%% ENVIRONNEMENTAL STOCHASTICITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load initialisation.mat
load prep.mat 
%load('pre_perturb.mat')
tmax = 3000 ;
tspan = 0 : 20 : tmax ;
options = odeset('RelTol',1e-3,'AbsTol',1e-4);


for vi = 1:100 %nb de simus
    
%for vv = 1:length(prep)
for vv = 1:2
    
%enregistrer var nécessaires pour les perturb

animal_density_last_step = prep{vv}.pre_perturb.animal_density_last_step ;
plant_density_last_step = prep{vv}.pre_perturb.plant_density_last_step ;
effort_ij0_last_step = prep{vv}.pre_perturb.effort_ij0_last_step ; 
A_foraging = prep{vv}.pre_perturb.A_foraging ;
 
%perturb
stochasticity

prep{vv}.stochasticity.animal_density_perturb = animal_density_perturb ;
prep{vv}.stochasticity.plant_density_perturb = plant_density_perturb ;
prep{vv}.stochasticity.effort_dynamic_perturb = effort_dynamic_perturb ;
prep{vv}.stochasticity.y0_dyn = y0_dyn ; 
prep{vv}.stochasticity.t_post = t_post ;

stoch{vv,vi} = prep{vv}.stochasticity ;

end

save stoch.mat stoch -v7.3


end



%% Displacement range of the niche optimum %%%%%%%%


%Code (tout est dedans)

options = odeset('RelTol',1e-3,'AbsTol',1e-4);

niche_displacement_range

%attention sortie animal_density_perturb ne correpsond pas à des pas de
%temps mais à chaque augmentation de mortalité

%% Vitesse de deplacement de l'optimum de niche %%%

options = odeset('RelTol',1e-4,'AbsTol',1e-5);

%Code (tout est dedans)
evol_foraging_niche_speed

%tout est dans evol_foraging_speed

save prep.mat prep -v7.3
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PROXIS and METRICS %%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    name = {'press_sigmaK_decrease' , 'press_sigmaK_increase', 'press_mortality_decrease', ...
                'press_mortality_increase', 'press_mortality_threshold', 'stochasticity' } ;  

        name_prox = {'biomass_post_p', 'biomass_post_a', 'biomass_post_tot', 'biomass_pre_p', 'biomass_pre_a', 'biomass_pre_tot' ...
                      'productivity_post', 'productivity_pre', ...
                      'FDis_post_p','FDis_post_a', 'FDis_pre_p', 'FDis_pre_a', ...
                      'FRO_pre_p', 'FRO_pre_a', 'FRO_post_p', 'FRO_post_a'} ;

for ii = 1:length(prep)
%     for ji = 1:length(name)
    for ji = 6 

%         % PROXIS %
%         % loading var
        animal_density = prep{ii}.pre_perturb.animal_density ;
        plant_density = prep{ii}.pre_perturb.plant_density ;
        effort_dynamic = prep{ii}.pre_perturb.effort_dynamic ;
        animal_density_perturb = prep{ii}.(name{ji}).animal_density_perturb ;
        plant_density_perturb = prep{ii}.(name{ji}).plant_density_perturb ;
        effort_dynamic_perturb = prep{ii}.(name{ji}).effort_dynamic_perturb ;
        t_pre = prep{ii}.pre_perturb.t_pre ;        
        t_post = prep{ii}.(name{ji}).t_post ;

        % compute proxis
        proxis_process
        
        % save results of proxis_process
        prep{ii}.(name{ji}).proxis = proxis ; 
        
%           for a = 1:length(name_prox) 
          for a = 7:8  
              % mean and variance of proxis
              prep{ii}.(name{ji}).proxis.([sprintf('%s', name_prox{a}), '_', 'mean']) = mean(prep{ii}.(name{ji}).proxis.(name_prox{a})) ;
              prep{ii}.(name{ji}).proxis.([sprintf('%s', name_prox{a}), '_', 'var']) = var(prep{ii}.(name{ji}).proxis.(name_prox{a})) ;
          end

        % METRICS 
        %%% ATT necessite que name soit dans le bon ordre! %%%
        if ji <= 4
 
           Resistance = resistance(biomass_pre_p,biomass_post_p, biomass_pre_a, biomass_post_a, biomass_pre_tot, ...
        biomass_post_tot, productivity_pre, productivity_post, productivity_pond_pre, ...
        productivity_pond_post, FDis_pre_p, FDis_post_p, FDis_pre_a, FDis_post_a, ...
        FRO_pre_p, FRO_post_p, FRO_pre_a, FRO_post_a, t_pre , t_post) ;
        
        %save results
        prep{ii}.(name{ji}).resistance = Resistance ;
    
       elseif ji==5
           Tolerance = tolerance(biomass_pre_p,biomass_post_p, biomass_pre_a, biomass_post_a, biomass_pre_tot, ...
         biomass_post_tot, productivity_pre, productivity_post, productivity_pond_pre, ...
         productivity_pond_post, FDis_pre_p, FDis_post_p, FDis_pre_a, FDis_post_a, ...
         FRO_pre_p, FRO_post_p, FRO_pre_a, FRO_post_a, t_pre , t_post) ;

         %save results
         prep{ii}.(name{ji}).tolerance = Tolerance ; 
     
        elseif ji==6   
           Stochasticity = stoch_metric(biomass_pre_p,biomass_post_p, biomass_pre_a, biomass_post_a, biomass_pre_tot, ...
         biomass_post_tot, productivity_pre, productivity_post, productivity_pond_pre, ...
         productivity_pond_post, FDis_pre_p, FDis_post_p, FDis_pre_a, FDis_post_a, ...
         FRO_pre_p, FRO_post_p, FRO_pre_a, FRO_post_a, t_pre , t_post) ;
      
         %save results
         prep{ii}.(name{ji}).stochasticity_metric = Stochasticity ;
     
        end        
    
    end
end

stoch_metric{vv,ii} = prep{vv}.stochasticity.stochasticity_metric ; % ça va pas ça écrase les densités 

save stoch_metric.mat stoch_metric
%%
%tondeuse = prep{1} ;
%af = prep{2} ;

save prep.mat prep -v7.3