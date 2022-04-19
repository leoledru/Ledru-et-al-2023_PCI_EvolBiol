%% Post_process: figures et table des metriques
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('prep.mat')
figure ()

% Pre-perturb
    for i=1:size(prep{1}.pre_perturb.plant_density,1)
       if mod(i,10)==0
         %  pause
          drawnow
        
            subplot(2,2,1)
            plot(prep{1}.pre_perturb.plant_density(i,:),'g') ;
            title('Dynamique des ressources pré-perturbation') ;
            
            subplot (2,2,2)
            plot(prep{1}.pre_perturb.animal_density(i,:),'g') ; %Tondeuse
            imagesc(reshape(prep{1}.pre_perturb.animal_density(i,:),number_of_animals,number_of_foraging)); %AF
     
            colorbar
            caxis([0 1])  %bloquer l'echelle de couleur
            title('Dynamique des consommateurs pré-perturbation ') ;
 
       end
    end

%% Post-perturb

    for i=1:size(prep{1, 1}.niche_displacement_range.animal_density_perturb,1)
        
        
        if mod(i,1)==0
            drawnow
%             pause 
            subplot(2,2,3)
            plot(prep{1, 1}.niche_displacement_range.plant_density_perturb(i,:),'g');
            title('Dynamique des ressources post-perturbation') ;
        
%         end
            
            subplot(2,2,4)
            imagesc(reshape(prep{1, 1}.niche_displacement_range.animal_density_perturb(i,:),number_of_animals,number_of_foraging));
            colorbar
            caxis([0 1]) 
            title('Dynamique des consomateurs post-perturbation') ;
    
        end
    end
    
%% Generate plots of proxis

%%Pour une meilleure visibilite, on ne garde que les 85% finaux du pre_perturb

name = {'pulse_on_plant' , 'pulse_on_animal' , 'pulse_on_all', 'press_sigmaK_decrease' , ...
        'press_sigmaK_increase', 'press_mortality_decrease', 'press_mortality_increase', ... 
         'press_mortality_threshold', 'stochasticity', 'niche_speed', 'niche_displacement_range' } ;  

        
name_title = {'Pulse on plant' , 'Pulse on animal' , 'Pulse on all', 'Press sigmaK decrease' , ...
              'Press sigmaK increase', 'Press mortality decrease', 'Press mortality increase', ...
              'Press mortality threshold', 'Stochasticity', 'Niche speed', 'Niche displacement range'} ; 

            
name_prox = { 'biomass_pre_p','FDis_pre_p', 'FRO_pre_p' ; ...
            'biomass_pre_a', 'FDis_pre_a', 'FRO_pre_a' ; ...
            'biomass_post_p','FDis_post_p', 'FRO_post_p' ;  ...
            'biomass_post_a','FDis_post_a', 'FRO_post_a'
            'biomass',        'FDis',       'FRO'     } ;  %last row for title
            

          
%reste biomass tot non plotée

nom = {'mower community', 'af community'} ;


for j = 1:2 %tondeuse ou af
    for i = 1: size(name,2)
        for k = 1 : size(name_prox, 2)
            

        if (strcmp(name{i},'press_mortality_threshold'))==1   %ATT necessite que name soit dans l'ordre 
            idx = find(prep{j}.pre_perturb.t_pre> 0.25*max(prep{j}.pre_perturb.t_pre)) ;
            t_pre = 1:size(prep{j}.pre_perturb.animal_density(idx)) ; %pour ne pas exprimer en fonctiond du temps
            T_post = (1:size(prep{j}.(name{i}).animal_density_perturb,1)) +  t_pre(end) ;

        elseif (strcmp(name{i},'niche_speed'))==1 || (strcmp(name{i},'niche_displacement_range'))==1 
            idx = find(prep{j}.(name{i}).t_pre> 0.25*max(prep{j}.(name{i}).t_pre)) ;
            t_pre = 1:size(prep{j}.pre_perturb.animal_density(idx)) ; % pas d'emergence
            T_post = (1:size(prep{j}.(name{i}).animal_density_perturb,1)) + t_pre(end) ;
            
        else 
            idx = find(prep{j}.pre_perturb.t_pre> 0.25*max(prep{j}.pre_perturb.t_pre)) ;
            t_pre = prep{j}.pre_perturb.t_pre(idx) ; % pas d'emergence
            T_post = (prep{j}.(name{i}).t_post) + t_pre(end) ;                

        end
        
       
        % Plants  

       figure('visible' , 'off')
%        figure()
        ji(1) = subplot(121)  ;
        
        sgtitle([sprintf('%s', name_title{1,i}), ' ', 'in', ' ', sprintf('%s', nom{j}), ' ', sprintf('%s', name_prox{5,k})], 'FontSize', 10) ;

        fig = plot(t_pre, prep{j}.(name{i}).proxis.(name_prox{1,k}), 'k ') ; hold on ; plot(T_post, prep{j}.(name{i}).proxis.(name_prox{3,k}), 'r') ;
        plot(t_pre(end),prep{j}.(name{i}).proxis.(name_prox{1,k})(end),'-s', 'MarkerSize', 3, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b') ; 
        plot(T_post(1), prep{j}.(name{i}).proxis.(name_prox{3,k})(1), '-s', 'MarkerSize', 3, 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'm') ;

        title('Ressources') 
        

        % Animals

        ji(2)= subplot(122);

        fig = plot(t_pre, prep{j}.(name{i}).proxis.(name_prox{2,k}), 'k ') ; hold on ; plot(T_post, prep{j}.(name{i}).proxis.(name_prox{4,k}), 'r') ;
        
        plot(t_pre, prep{j}.(name{i}).proxis.(name_prox{2,k}), 'k ') ; hold on ; plot(T_post, prep{j}.(name{i}).proxis.(name_prox{4,k}), 'r') ;
        plot(t_pre(end),prep{j}.(name{i}).proxis.(name_prox{2,k})(end),'-s', 'MarkerSize', 3, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b') ; 
        plot(T_post(1), prep{j}.(name{i}).proxis.(name_prox{4,k})(1), '-s', 'MarkerSize', 3, 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'm') ;
 
        title('Consumers')

        %linkaxes(ji) ;
        filename = [sprintf('%s', nom{j}), '_', sprintf('%s', name_title{1,i}), '_', sprintf('%s', name_prox{5,k})] ;
        fpath = '~/Bureau/Etudes/Stages/Stage_LECA/Eco-evo_foraging/Ecosystem/Figures/perturb' ;
        saveas(fig, fullfile(fpath, filename), 'png')

        end
        
        
        %Productivty %
        
        figure('visible', 'off')
        
        fig = plot(t_pre, prep{j}.(name{i}).proxis.productivity_pre, 'k ') ; hold on ; plot(T_post, prep{j}.(name{i}).proxis.productivity_post, 'r') ;
        plot(t_pre(end),prep{j}.(name{i}).proxis.productivity_pre(end),'-s', 'MarkerSize', 3, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b') ; 
        plot(T_post(1), prep{j}.(name{i}).proxis.productivity_post(1), '-s', 'MarkerSize', 3, 'MarkerFaceColor', 'm', 'MarkerEdgeColor', 'm') ;

        title([sprintf('%s', name_title{1,i}), ' ', 'in a', ' ', sprintf('%s', nom{j}), ' ', 'productivity'] , 'FontSize', 10) ;
  
        filename = [sprintf('%s', nom{j}),' ', sprintf('%s', name_title{1,i}), ' ', sprintf('productivity')] ;
        fpath = '~/Bureau/Etudes/Stages/Stage_LECA/Eco-evo_foraging/Ecosystem/Figures' ;
        saveas(fig, fullfile(fpath, filename), 'png')

        
        %prod + ponderee %%%%%%%%%%
             %%%%%%%%%%%%%%%%%%%%%%%%%
        
        
    end
end

           


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% METRICS %%%%
%%%%%%%%%%%%%%%

%% RESILIENCE FOR PULSE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

type_p = {'pulse_on_plant', 'pulse_on_animal', 'pulse_on_all'
          'pulse on plant', 'pulse on animal', 'pulse on all' } ;

for i= 1:length(type_p)
    
    %remplacer les never leave equilibrum par -1
    
    plot_p_t = prep{1}.(type_p{i}).resilience_all_simus.resilience_p_B ;     
    plot_p_af = prep{2}.(type_p{i}).resilience_all_simus.resilience_p_B ;     
    
    plot_a_t = prep{1}.(type_p{i}).resilience_all_simus.resilience_a_B ;     
    plot_a_af = prep{2}.(type_p{i}).resilience_all_simus.resilience_a_B ;     

    
    if isstring(plot_p_t) ==1
        plot_p_t(plot_p_t == "never leave equilibrum") = -1 ;
        plot_p_t = str2double(plot_p_t) ;
        
    elseif isstring(plot_p_af) ==1    
        plot_p_af(plot_p_af == "never leave equilibrum") = -1 ;
        plot_p_af = str2double(plot_p_af) ;
        
    elseif isstring(plot_a_af) ==1    
        plot_a_af(plot_a_af == "never leave equilibrum") = -1 ;
        plot_a_af = str2double(plot_a_af) ;
        
    elseif isstring(plot_a_t) ==1     
         plot_a_t(plot_a_t == "never leave equilibrum") = -1 ;
        plot_a_t = str2double(plot_a_t) ;
        
    end
    
    
    %Biomass
    figure
    sgtitle(['Resilience of communities after a ', sprintf('%s', type_p{2,i})], 'FontSize', 10)
    
    subplot(121)
    boxplot(plot_p_t , plot_p_af, 'Labels', {'tondeuse', 'af'}) 
    title('Ressources')

    subplot(122)
    boxplot(plot_a_t , plot_a_af )
    title('Consumers')

    %Productivité
    
    %FRO
    
    %FDIs

end


%% RESISTANCE %%


%% TOLERANCE

name = {'press_mortality_threshold', 'niche_displacement_range' , 'niche_speed'} ;  
name_prox = { 'biomass', 'FDis', 'FRO'}; ...


for j = 1: length(name)
    for k = 1:length(name_prox)
            
        figure()
        sgtitle(['tolerance of ', ' ', sprintf('%s', name_prox{k}), ' to ',  sprintf('%s', name{j})], 'FontSize', 10) ;
        
        subplot(121);
        fig = plot(prep{1}.(name{j}).tolerance.(name_prox{k}).tolerance_p) ; hold on ; plot(prep{2}.(name{j}).tolerance.(name_prox{k}).tolerance_p) %tondeuse
       
        text(length(prep{1}.(name{j}).tolerance.(name_prox{k}).tolerance_p), prep{1}.(name{j}).tolerance.(name_prox{k}).tolerance_p(end), sprintf('%d', length(prep{1}.(name{j}).tolerance.(name_prox{k}).tolerance_p)), 'Color', 'blue')
        text(length(prep{2}.(name{j}).tolerance.(name_prox{k}).tolerance_p), prep{2}.(name{j}).tolerance.(name_prox{k}).tolerance_p(end), sprintf('%d', length(prep{2}.(name{j}).tolerance.(name_prox{k}).tolerance_p)), 'Color', 'red')

        title('Ressources')
        legend('tondeuse', 'af', 'Location', 'best', 'FontSize', 10 )       
        
        
        subplot(122)
        fig = plot(prep{1}.(name{j}).tolerance.(name_prox{k}).tolerance_a) ; hold on ; plot(prep{2}.(name{j}).tolerance.(name_prox{k}).tolerance_a)
        
        text(length(prep{1}.(name{j}).tolerance.(name_prox{k}).tolerance_a), prep{1}.(name{j}).tolerance.(name_prox{k}).tolerance_a(end), sprintf('%d', length(prep{1}.(name{j}).tolerance.(name_prox{k}).tolerance_a)), 'Color', 'blue')
        text(length(prep{2}.(name{j}).tolerance.(name_prox{k}).tolerance_a), prep{2}.(name{j}).tolerance.(name_prox{k}).tolerance_a(end), sprintf('%d', length(prep{2}.(name{j}).tolerance.(name_prox{k}).tolerance_a)), 'Color', 'red')

        title('Consummers')
        legend('tondeuse', 'af', 'Location', 'best')
        

%         filename = ['tolerance', '_', sprintf('%s', name{j}), '_', sprintf('%s', name_prox{k})] ;
%         fpath = '~/Bureau/Etudes/Stages/Stage_LECA/Eco-evo_foraging/Ecosystem/Figures/perturb' ;
%         saveas(fig, fullfile(fpath, filename), 'png')
        
        % Productivity
        
        figure()
        sgtitle(['tolerance of ', ' ', sprintf('%s', name_prox{k}), ' to ',  sprintf('%s', name{j})], 'FontSize', 10) ;
               
        fig = plot(prep{1}.(name{j}).tolerance.productivity.tolerance_a) ; hold on ; plot(prep{2}.(name{j}).tolerance.productivity.tolerance_a)
        
        text(length(prep{1}.(name{j}).tolerance.productivity.tolerance_a), prep{1}.(name{j}).tolerance.productivity.tolerance_a(end), sprintf('%d', length(prep{1}.(name{j}).tolerance.productivity.tolerance_a)), 'Color', 'blue')
        text(length(prep{2}.(name{j}).tolerance.productivity.tolerance_a), prep{2}.(name{j}).tolerance.productivity.tolerance_a(end), sprintf('%d', length(prep{2}.(name{j}).tolerance.productivity.tolerance_a)), 'Color', 'red')

        title('Consummers')
        legend('tondeuse', 'af', 'Location', 'best', 'FontSize', 10)

    end
end
        
  
       


