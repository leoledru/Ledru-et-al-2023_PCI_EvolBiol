%% CODE CONTENANT : EMERGENCE DE COMMU PRE-PERTURB, PUIS CHAQUE PERTURB
%% EN APPELANT A CHAQUE FOIS LE CODE CORRESPONDANT A LA PERTURB.

clear all; clc;
%%%
% Version Oceane avec l'introduction de pertubations basée sur la version 
% avec la méthode de calcul des efforts de Jimmy (matrice de
% transition Q). Ici pour les efforts et delta, les plantes sont en lignes
% et les animaux en colonnes
%%%%


%% %%%%%%%%%%%%%%%%%%%%%% INITIALISATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Paramètres communs pour ODE pre & post
tmax = 1000 ;
tspan = 0:20:tmax ;
% tspan = [0,tmax] ;
% Df = 0 

%%%

initialisation_pre_perturb


% save(['simu_post_perturb_increased_threshold','.mat']) ;

%% RESULTATS SIMUS PRE ENREGISTRES

load('simu_initiale_preperturb_t') ; %tmax = 1000
load('simu_initiale_preperturb_af') ;
load('simu_10000_t') ;


load('simu_ref_preperturb_af') ; %att simu avec y0 et taille niche modifiés pour amplitude 
load('simu_ref_preperturb_t') ; %a reenregistrer

load('simu_post_perturb_10percent') ;
load('simu_post_perturb_pulse') ;
load('simu_post_perturb_increased_threshold') ;

%% %%%%%%%%%%%%%%%%%%%%
%%% PERTURBATIONS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%% PULSE PERTURBATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% !!!!!! meme pulse pour tondeuse et af !!!!!!!!

pulse_perturb %code


%% %%%%%%%%%%%%%%%%%%%% PRESS PERTURBATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Press perturbation on sigma_K
sigmaK_perturb_code

%% Press perturbation on mortality (+10%)
increased_mortality_10percent %code

%% Increasing mortality to global extinction ("seuil")
increasing_mortality_threshold

%% Displacement range of the niche optimum ("seuil")
niche_displacement_range

%% Vitesse de deplacement de l'optimum de niche ("seuil")

% ATT pas le même tmax !!!

evol_foraging_niche_speed

t_post = t ;


%% %%%%%%%%%%%% ENVIRONNEMENTAL STOCHASTICITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stochasticity

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%% CALCUL DES PROXIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
proxis 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%% METRIQUES %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%%% PULSE PERTURBATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%REACTIVITY (neubert, cited in Kefi)

%% RESILIENCE R
clearvars Resilience ; %necessite que le str n'existe pas pour la créer

Resilience = resilience(biomass_pre_p,biomass_post_p, biomass_pre_a, biomass_post_a, biomass_pre_tot, ...
    biomass_post_tot, productivity_pre, productivity_post, productivity_pond_pre, ...
    productivity_pond_post, FDis_pre_p, FDis_post_p, FDis_pre_a, FDis_post_a, ...
    FRO_pre_p, FRO_post_p, FRO_pre_a, FRO_post_a, t_pre , t_post) ;

%% %%%%%%%%%%%%%%%%%%%% PRESS PERTURBATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% RESISTANCE (à l'augmentation unique de la mortalité RM ou diminution largeur de niche)
clearvars Resistance ;

Resistance = resistance(biomass_pre_p,biomass_post_p, biomass_pre_a, biomass_post_a, biomass_pre_tot, ...
    biomass_post_tot, productivity_pre, productivity_post, productivity_pond_pre, ...
    productivity_pond_post, FDis_pre_p, FDis_post_p, FDis_pre_a, FDis_post_a, ...
    FRO_pre_p, FRO_post_p, FRO_pre_a, FRO_post_a) ; 

%Resistance.RG


%% TOLERANCE to increased mortality TM (augmentation ou amplitude) / to y0 deplacement (amplitude ou vitesse)

clearvars Tolerance ; 

Tolerance = tolerance(biomass_pre_p,biomass_post_p, biomass_pre_a, biomass_post_a, biomass_pre_tot, ...
    biomass_post_tot, productivity_pre, productivity_post, productivity_pond_pre, ...
    productivity_pond_post, FDis_pre_p, FDis_post_p, FDis_pre_a, FDis_post_a, ...
    FRO_pre_p, FRO_post_p, FRO_pre_a, FRO_post_a, t_pre, t_post) ;

%% STOCHASTICITY METRICS

clearvars stoch_metric ;

stoch_metric = stoch_metric(biomass_pre_p,biomass_post_p, biomass_pre_a, biomass_post_a, biomass_pre_tot, ...
    biomass_post_tot, productivity_pre, productivity_post, productivity_pond_pre, ...
    productivity_pond_post, FDis_pre_p, FDis_post_p, FDis_pre_a, FDis_post_a, ...
    FRO_pre_p, FRO_post_p, FRO_pre_a, FRO_post_a, t_pre, t_post) ;


 %% %%%%%%%%%%%%%%%%%%%%%%%%% FIGURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure ()

% Pre-perturb
    for i=1:size(plant_density,1)
       if mod(i,100)==0
%             pause
          drawnow
        
            subplot(2,2,1)
            plot(plant_density(i,:),'g') ;
            title('Dynamique des ressources pré-perturbation') ;
            
            subplot (2,2,2)
            plot(animal_density(i,:),'g') ; %Tondeuse
%             imagesc(reshape(animal_density(i,:),number_of_animals,number_of_foraging)); %AF
     
            colorbar
            caxis([0 1])  %bloquer l'echelle de couleur
            title('Dynamique des consommateurs pré-perturbation ') ;
 
       end
    end

% Post-perturb
    for i=1:size(plant_density_perturb,1)
%         if mod(i,100)==0
            drawnow
%             pause 
            subplot(2,2,3)
            plot(plant_density_perturb(i,:),'g');
            title('Dynamique des ressources post-perturbation') ;
        
%         end
            
            subplot(2,2,4)
            imagesc(reshape(animal_density_perturb(i,:),number_of_animals,number_of_foraging));
            colorbar
            caxis([0 1]) 
            title('Dynamique des consomateurs post-perturbation') ;
    end

%%
figure() %ATT choose between tondeuse and AF

T_post = t_post + (t_pre(end)+1) ;


%Plant
subplot(121)  

plot(t_pre, biomass_pre_p, 'k ') ; hold on ; plot(T_post, biomass_post_p, 'r') ; %Tondeuse
plot(t_pre(end),biomass_pre_p(end),'-s', 'MarkerSize', 7, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b') ; 
plot(T_post(1), biomass_post_p(1), '-s', 'MarkerSize', 7, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b') ;

%plot(t_pre, movmean_pre_p) ; hold on ; plot(t_post, movmean_post_p) ; %AF

title('Biomasse des ressources') ;

%Animal
subplot(122)

plot(t_pre, biomass_pre_a, 'k') ; hold on ; plot(T_post, biomass_post_a, 'r') ; %Tondeuse
plot(t_pre(end),biomass_pre_a(end),'-s', 'MarkerSize', 7, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b') ; 
plot(T_post(1), biomass_post_a(1), '-s', 'MarkerSize', 7, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b') ;

%plot(t_pre, movmean_pore_a) ; hold on ; plot(t_post, movmean_post_a) ; %AF

title('Biomasse des consommateurs') ;

% %All plotted
% %Plant
% plot(t_pre, biomass_pre_p, 'k', 'LineStyle', ':') ; hold on ; plot(t_post, biomass_post_p, 'r', 'LineStyle', ':') ; %Tondeuse
% %Animals
% plot(t_pre, biomass_pre_a, 'k') ; hold on ; plot(t_post, biomass_post_a, 'r') ; %AF


%% Productivity
figure() 

plot(t_pre, productivity_pre, 'k ') ; hold on ; plot(T_post, productivity_post, 'r') ; %Tondeuse
plot(t_pre(end),productivity_pre(end),'-s', 'MarkerSize', 7, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b') ; 
plot(T_post(1), productivity_post(1), '-s', 'MarkerSize', 7, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b') ;
