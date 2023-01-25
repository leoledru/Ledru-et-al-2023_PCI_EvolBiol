clear all;clc
%%
number_of_animals = 11;
number_of_plants  = 11;
number_of_foraging = 11;

xmin = -5;
xmax = 5;
xx = linspace(xmin,xmax,number_of_animals);
dx = xx(2)-xx(1);
dx = 1;
traits_of_animals = xx';
foraging_trait = linspace(0,1,number_of_foraging);

% load('community_comparison6.mat');
% load('community_comparison_CK.mat');
load('community_comparison.mat');

% each column is a simulation, row 1 = animal simu AF; row 2 = plant simu AF;
% row 3 = Effort simu AF; row 4 = animal simu tondeuse; row 5 = plant simu
% tondeuse; row 6 = Effort tondeuse; row 7 = params

% postprocess = cell(7,size(Param,2));
% load('postprocess_densities')
%%
for ii = 1:size(Param,2)
% for ii = 1:1
  postprocess{1,ii} = Output_a_f(:,(ii-1)*number_of_animals*number_of_foraging+1:ii*number_of_animals*number_of_foraging);
  postprocess{2,ii} = Output_p_f(:,(ii-1)*number_of_plants+1:ii*number_of_plants);
%   postprocess{3,ii} = Effort_f(:,(ii-1)*number_of_animals*number_of_foraging+1:ii*number_of_animals*number_of_foraging*number_of_plants);
  postprocess{4,ii} = Output_a_t(:,(ii-1)*number_of_animals*number_of_foraging+1:ii*number_of_animals*number_of_foraging);
  postprocess{5,ii} = Output_p_t(:,(ii-1)*number_of_plants+1:ii*number_of_plants);
%   postprocess{6,ii} = Effort_t(:,(ii-1)*number_of_animals*number_of_foraging+1:ii*number_of_animals*number_of_foraging*number_of_plants);
  postprocess{7,ii} = Param(:,ii); 
end
%%
% Check for global collapse ?
sum_simu = cell(4,size(postprocess,2));
for ii = 1:size(postprocess,2)
    sum_simu{1,ii} = sum(postprocess{1,ii},'all'); 
    sum_simu{2,ii} = sum(postprocess{2,ii},'all'); 
    sum_simu{3,ii} = sum(postprocess{4,ii},'all'); 
    sum_simu{4,ii} = sum(postprocess{5,ii},'all'); 
end
sum_simu = cell2mat(sum_simu);

check_for_collapse = sum(sum_simu==0,'all');
if check_for_collapse == 0
    message2 = sprintf('No simulation collapse')
%     msgbox(message2);
else
    message1 = sprintf('The simulation %d collapse\n',find(sum(sum_simu==0)))
%     msgbox(message1);
end

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comparaison of final states with/without evolution of AF

% count = 6; % number of simulations to look
% for ii = 1:count
%     final_animal_f = postprocess{1,ii};
%     final_plant_f = postprocess{2,ii};
%     final_animal_t = postprocess{4,ii};
%     final_plant_t = postprocess{5,ii};
% 
%     final_animal_f = reshape(final_animal_f(end,:),number_of_animals,number_of_foraging);
%     final_animal_t = reshape(final_animal_t(end,:),number_of_animals,number_of_foraging);
%     final_plant_f = final_plant_f(end,:);
%     final_plant_t = final_plant_t(end,:);
% 
%     subplot(count,3,ii*3-2)
%     plot(xx,sum(final_animal_f,1),'b',xx,sum(final_animal_t,1),'r');
%     title('foraging animal')
%     subplot(count,3,ii*3-1)
%     plot(foraging_trait,sum(final_animal_f,2),'b',foraging_trait,sum(final_animal_t,2),'r');
%     title('niche animal')
%     subplot(count,3,ii*3)
%     plot(xx,final_plant_f,xx,final_plant_t);
%     title('niche plant')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% All final states superposed (hold on)

for ii = 1:size(postprocess,2)
    final_animal_f = postprocess{1,ii};
    final_plant_f = postprocess{2,ii};
    final_animal_t = postprocess{4,ii};
    final_plant_t = postprocess{5,ii};

    final_animal_f = reshape(final_animal_f(end,:),number_of_animals,number_of_foraging);
    final_animal_t = reshape(final_animal_t(end,:),number_of_animals,number_of_foraging);
    final_plant_f = final_plant_f(end,:);
    final_plant_t = final_plant_t(end,:);

    % DENSITIES
    subplot(3,2,1)
    plot(foraging_trait,sum(final_animal_f,1),'b')
    hold on
    title('Consumers foraging stratgies','interpreter','latex')
    ylabel('Abundances','interpreter','latex')
    xlabel('Foraging trait','interpreter','latex')
    subplot(3,2,2)
    plot(foraging_trait,sum(final_animal_t,1),'r');
    subplot(3,2,3)
    plot(foraging_trait,sum(final_animal_f,2),'b');
    hold on
    title('consumers niche position','interpreter','latex')
    ylabel('Abundances','interpreter','latex')
    xlabel('Niche trait','interpreter','latex')
    subplot(3,2,4)
    plot(foraging_trait,sum(final_animal_t,2),'r');
    hold on
    subplot(3,2,5)
    plot(xx,final_plant_f,'b');
    hold on
    title('resources niche position','interpreter','latex')
    subplot(3,2,6)
    plot(xx,final_plant_t,'r');
    hold on
%     
    % ONLY PEAKS
%     peaks_animal_f = islocalmax(sum(final_animal_f,2));
%     peaks_animal_t = islocalmax(sum(final_animal_t,2));
%     peaks_plant_f = islocalmax(final_plant_f);
%     peaks_plant_t = islocalmax(final_plant_t);
%     peaks_foraging_f = islocalmax([0,sum(final_animal_f,1),0]);
%     peaks_foraging_t = islocalmax([0,sum(final_animal_t,1),0]);
%     peaks_foraging_f = peaks_foraging_f(2:end-1);
%     peaks_foraging_t = peaks_foraging_t(2:end-1);
%     
%     subplot(3,2,1)
%     plot(xx,sum(final_animal_f,2).*peaks_animal_f,'*b')
%     xlabel('niche trait','interpreter','latex')
%     ylabel('densities','interpreter','latex')
%     hold on
%     title('Animals','interpreter','latex')
%     subplot(3,2,2)
%     plot(xx,sum(final_animal_t,2).*peaks_animal_t,'*r')
%     hold on    
%     subplot(3,2,3)
%     plot(xx,final_plant_f.*peaks_plant_f,'*b')
%     xlabel('niche trait','interpreter','latex')
%     ylabel('densities','interpreter','latex')
%     hold on
%     title('Plants','interpreter','latex')
%     subplot(3,2,4)
%     plot(xx,final_plant_t.*peaks_plant_t,'*r')
%     hold on   
%     subplot(3,2,5)
%     plot(foraging_trait,sum(final_animal_f,1).*peaks_foraging_f,'*b')
%     xlabel('foraging trait','interpreter','latex')
%     ylabel('densities','interpreter','latex')
%     hold on
%     title('Foraging strategies','interpreter','latex')
%     subplot(3,2,6)
%     plot(foraging_trait,sum(final_animal_t,1).*peaks_foraging_t,'*r')
%     hold on
end

sgtitle('Aggregation of all final states ($\sigma_K > \sigma_C)$','interpreter','latex')
%%

% Compute the number of species averaged on all simulations : peaks in trait-space

% Reshape animal density (postprocess rows 1 and 4)
Out_f = cell(1,size(postprocess,2));
Out_t = cell(1,size(postprocess,2));
Out_F = cell(1,size(postprocess,2));
Out_T = cell(1,size(postprocess,2));
for ii = 1:size(postprocess,2)
    Output_a_f = [];
    Output_a_t = [];
    output_a_f = postprocess{1,ii};
    output_a_t = postprocess{4,ii};
    for jj = 1:size(output_a_f,1)
        Output_a_f = [Output_a_f;reshape(output_a_f(jj,:),number_of_animals,number_of_foraging)];
        Output_a_t = [Output_a_t;reshape(output_a_t(jj,:),number_of_animals,number_of_foraging)];    
    end
    Out_F{1,ii} = Output_a_f;
    Out_T{1,ii} = Output_a_t;
end
%
% Sum animal densities on foraging trait
animal_niche_f = [];
for ii = 1:size(postprocess,2)
    animal_niche_f = [];
    animal_f = Out_F{:,ii};
    animal_niche_t = [];
    animal_t = Out_T{:,ii};
    for jj = 1:101
    animal_niche_f = [animal_niche_f;sum(animal_f((jj-1)*number_of_foraging+1:jj*number_of_foraging,:),2)'];
    animal_niche_t = [animal_niche_t;sum(animal_t((jj-1)*number_of_foraging+1:jj*number_of_foraging,:),2)'];
    end
    Out_f{:,ii} = animal_niche_f;
    Out_t{:,ii} = animal_niche_t;
end

% % Compute peaks for each simulation
% peak_animal_f = cell(1,size(Param,2));
% peak_animal_t = peak_animal_f;
% peak_plant_f = peak_animal_f;
% peak_plant_t = peak_animal_f;
% 
% for ii = 1:size(Param,2)
%     peak_animal_f{:,ii} = islocalmax(Out_f{:,ii},2);
%     peak_animal_t{:,ii} = islocalmax(Out_t{:,ii},2);
%     peak_plant_f{:,ii} = islocalmax(postprocess{2,ii},2);
%     peak_plant_t{:,ii} = islocalmax(postprocess{5,ii},2);
% end
% 
% % Plot of peaks dynamic (for each simulation)
% clf
% time = 1:101;
% 
% for ii = 1:size(Param,2)
% % for ii = 5
%     subplot(1,2,1)
%     plot(time,sum(peak_animal_f{:,ii},2),'*b',time,sum(peak_animal_t{:,ii},2),'or');
%     hold on
%     title('animal peaks')
%     subplot(1,2,2)
%     plot(time,sum(peak_plant_f{:,ii},2),'*b',time,sum(peak_plant_t{:,ii},2),'or');
%     hold on
%     title('plant peaks')
% end
%%
%%%%%%%%%%%%%%%%%%%%%
% SPECIFIC RICHNESS %
%%%%%%%%%%%%%%%%%%%%%

% Sum of peaks of all simulations
% peak_animal_f = zeros(101,number_of_animals);
% peak_animal_t = peak_animal_f;
% peak_plant_f = peak_animal_f;
% peak_plant_t = peak_animal_f;
% 
% for ii = 1:size(Param,2)
%     peak_animal_f = peak_animal_f + islocalmax(Out_f{:,ii},2);
%     peak_animal_t = peak_animal_t + islocalmax(Out_t{:,ii},2);
%     peak_plant_f = peak_plant_f + islocalmax(postprocess{2,ii},2);
%     peak_plant_t = peak_plant_t + islocalmax(postprocess{5,ii},2);
% end

% Sum of all traits>0 of all simulations
richness_animal_f = zeros(101,number_of_animals);
richness_animal_t = richness_animal_f;
richness_plant_f = richness_animal_f;
richness_plant_t = richness_animal_f;

for ii = 1:size(Param,2)
    richness_animal_f = richness_animal_f + (Out_f{:,ii}>0);
    richness_animal_t = richness_animal_t + (Out_t{:,ii}>0);
    richness_plant_f = richness_plant_f + (postprocess{2,ii}>0);
    richness_plant_t = richness_plant_t + (postprocess{5,ii}>0);
end

% Plot of peaks dynamic (averaged on all simulations)
% time = 1:101;
% subplot(2,2,1)
% plot(time,sum(peak_animal_f(:,:),2)/size(postprocess,2),'b',time,sum(peak_animal_t(:,:),2)/size(postprocess,2),'r');
% xlim([0 101])
% xlabel('Time','interpreter','latex')
% ylabel('Number of species','interpreter','latex')
% title('Specific richness of consumers peaks','interpreter','latex')
% legend('with AF evolution','without','interpreter','latex');
% subplot(2,2,2)
% plot(time,sum(peak_plant_f(:,:),2)/size(postprocess,2),'b',time,sum(peak_plant_t(:,:),2)/size(postprocess,2),'r');
% xlim([0 101])
% title('Specific richness of resources peaks','interpreter','latex')
% 
% subplot(2,2,3)
% plot(time,sum(richness_animal_f(:,:),2)/size(postprocess,2),'*b',time,sum(richness_animal_t(:,:),2)/size(postprocess,2),'r');
% xlim([0 101])
% xlabel('Time','interpreter','latex')
% ylabel('Number of species','interpreter','latex')
% title('Specific richness of all consumers','interpreter','latex')
% subplot(2,2,4)
% plot(time,sum(richness_plant_f(:,:),2)/size(postprocess,2),'*b',time,sum(richness_plant_t(:,:),2)/size(postprocess,2),'r');
% xlim([0 101])
% title('Specific richness of all resources','interpreter','latex')

% sgtitle('Specific Richness','interpreter','latex')
%%
%%%%%%%%%%%
% BIOMASS %
%%%%%%%%%%%

% Compute the total biomass of all simulations
biomass_animal_f = 0;
biomass_animal_t = 0;
biomass_plant_f = 0;
biomass_plant_t = 0;
for ii = 1:size(Param,2)
%     biomass_animal_f = biomass_animal_f + sum(postprocess{1,ii},2);
%     biomass_animal_t = biomass_animal_t + sum(postprocess{4,ii},2);
%     biomass_plant_f = biomass_plant_f + sum(postprocess{2,ii},2);
%     biomass_plant_t = biomass_plant_t + sum(postprocess{5,ii},2);

    peak_animal_f = islocalmax(Out_f{:,ii},2);
    peak_animal_t = islocalmax(Out_t{:,ii},2);
    peak_plant_f = islocalmax(postprocess{2,ii},2);
    peak_plant_t = islocalmax(postprocess{5,ii},2);
    biomass_animal_f = biomass_animal_f + sum(Out_f{:,ii}.*peak_animal_f,2);
    biomass_animal_t = biomass_animal_t + sum(Out_t{:,ii}.*peak_animal_t,2);
    biomass_plant_f = biomass_plant_f + sum(postprocess{2,ii}.*peak_plant_f,2);
    biomass_plant_t = biomass_plant_t + sum(postprocess{5,ii}.*peak_plant_t,2);
end

subplot(2,2,1)
plot(time,biomass_animal_f/size(postprocess,2),'b',time,biomass_animal_t/size(postprocess,2),'r');
xlim([0 101])
title('Biomass of consumers peaks','interpreter','latex')
legend('with AF evolution','without','interpreter','latex');
xlabel('Time','interpreter','latex');ylabel('Biomass','interpreter','latex');
subplot(2,2,2)
plot(time,biomass_plant_f/size(postprocess,2),'b',time,biomass_plant_t/size(postprocess,2),'r');
xlim([0 101])
title('Biomass of resources peaks','interpreter','latex')

biomass_animal_f = 0;
biomass_animal_t = 0;
biomass_plant_f = 0;
biomass_plant_t = 0;
for ii = 1:size(Param,2)
    biomass_animal_f = biomass_animal_f + sum(postprocess{1,ii},2);
    biomass_animal_t = biomass_animal_t + sum(postprocess{4,ii},2);
    biomass_plant_f = biomass_plant_f + sum(postprocess{2,ii},2);
    biomass_plant_t = biomass_plant_t + sum(postprocess{5,ii},2);
end

subplot(2,2,3)
plot(time,biomass_animal_f/size(postprocess,2),'b',time,biomass_animal_t/size(postprocess,2),'r');
xlim([0 101])
xlabel('Time','interpreter','latex');ylabel('Biomass','interpreter','latex');
title('Biomass of all consumers','interpreter','latex')
subplot(2,2,4)
plot(time,biomass_plant_f/size(postprocess,2),'b',time,biomass_plant_t/size(postprocess,2),'r');
xlim([0 101])
title('Biomass of all resources','interpreter','latex')

% sgtitle('Biomass summed on all traits','interpreter','latex')
% sgtitle('Biomass summed only on peaks traits','interpreter','latex
% sgtitle('Biomass','interpreter','latex')
%%
% BOXPLOT (average on the last 100 time steps)

biomass_animal_f = [];
biomass_animal_t = [];
biomass_plant_f = [];
biomass_plant_t = [];
for ii = 1:size(Param,2)
    biomass_animal_f = [biomass_animal_f;sum(postprocess{1,ii},2)'];
    biomass_animal_t = [biomass_animal_t;sum(postprocess{4,ii},2)'];
    biomass_plant_f = [biomass_plant_f;sum(postprocess{2,ii},2)'];
    biomass_plant_t = [biomass_plant_t;sum(postprocess{5,ii},2)'];
end

Biomass_plant_f = mean(biomass_plant_f(:,end-10:end),2);
Biomass_plant_t = mean(biomass_plant_t(:,end-10:end),2);
Biomass_plant_boxplot = [Biomass_plant_f,Biomass_plant_t];

Biomass_animal_f = mean(biomass_animal_f(:,end-10:end),2);
Biomass_animal_t = mean(biomass_animal_t(:,end-10:end),2);
Biomass_animal_boxplot = [Biomass_animal_f,Biomass_animal_t];

subplot(1,2,1)
boxplot(Biomass_animal_boxplot)
title('Total consumers biomass','interpreter','latex','FontSize',18)
xticks([1 2 3 4])
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTickLabel', {'with foraging','without'});
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)

subplot(1,2,2)
boxplot(Biomass_plant_boxplot)
title('Total resources biomass','interpreter','latex','FontSize',18)
xticks([1 2 3 4])
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTickLabel', {'with foraging','without'});
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)

%%
% coefficient of variation (average on the last 100 time steps)

Biomass_plant_f_cv = std(biomass_plant_f(:,end-10:end),0,2)./mean(biomass_plant_f(:,end-10:end),2);
Biomass_plant_t_cv = std(biomass_plant_t(:,end-10:end),0,2)./mean(biomass_plant_t(:,end-10:end),2);
Biomass_plant_boxplot_cv = [Biomass_plant_f_cv,Biomass_plant_t_cv];

Biomass_animal_f_cv = std(biomass_animal_f(:,end-10:end),0,2)./mean(biomass_animal_f(:,end-10:end),2);
Biomass_animal_t_cv = std(biomass_animal_t(:,end-10:end),0,2)./mean(biomass_animal_t(:,end-10:end),2);
Biomass_animal_boxplot_cv = [Biomass_animal_f_cv,Biomass_animal_t_cv];

subplot(1,2,1)
boxplot(Biomass_animal_boxplot_cv)
title('Coefficient of variation of consumers biomass','interpreter','latex','FontSize',18)
xticks([1 2 3 4])
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTickLabel', {'with foraging','without'});
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)

subplot(1,2,2)
boxplot(Biomass_plant_boxplot_cv)
title('Coefficient of variation of resources biomass','interpreter','latex','FontSize',18)
xticks([1 2 3 4])
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTickLabel', {'with foraging','without'});
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)


%%
%%%%%%%%%%%%%%%%%%%%%
% SHANNON / SIMPSON %
%%%%%%%%%%%%%%%%%%%%%

% Compute Shannon index on densities peaks (species) 

Shannon_animal_f = [];
Shannon_plant_f = [];
Shannon_animal_t = [];
Shannon_plant_t = [];

Simpson_animal_f = [];
Simpson_animal_t = [];
Simpson_plant_f = [];
Simpson_plant_t = [];

for ii = 1:size(postprocess,2)
    
    peak_animal_f = islocalmax(Out_f{:,ii},2);
    peak_animal_t = islocalmax(Out_t{:,ii},2);
    peak_plant_f = islocalmax(postprocess{2,ii},2);
    peak_plant_t = islocalmax(postprocess{5,ii},2);

    % Only peaks
%     animal_density_f = Out_f{:,ii}.*peak_animal_f;
%     plant_density_f = postprocess{2,ii}.*peak_plant_f;
%     animal_density_t = Out_t{:,ii}.*peak_animal_t;
%     plant_density_t = postprocess{5,ii}.*peak_plant_t;
    
    % All traits
    animal_density_f = Out_f{:,ii};
    plant_density_f = postprocess{2,ii};
    animal_density_t = Out_t{:,ii};
    plant_density_t = postprocess{5,ii};
    
    %%%%%%%%%%
    shannon_animal_f = zeros(101,number_of_animals);
    shannon_plant_f = zeros(101,number_of_plants);
    shannon_animal_t = zeros(101,number_of_animals);
    shannon_plant_t = zeros(101,number_of_plants);
    
    simpson_animal_f = zeros(101,number_of_animals);
    simpson_animal_t = zeros(101,number_of_animals);
    simpson_plant_f = zeros(101,number_of_animals);
    simpson_plant_t = zeros(101,number_of_animals);

    for jj = 1:number_of_animals
        shannon_animal_f(:,jj) = animal_density_f(:,jj)./sum(animal_density_f(:,:),2)...
                                    .*log2(animal_density_f(:,jj)./sum(animal_density_f(:,:),2));
        shannon_animal_f(isnan(shannon_animal_f)) = 0;

        shannon_plant_f(:,jj) = plant_density_f(:,jj)./sum(plant_density_f(:,:),2)...
                                    .*log2(plant_density_f(:,jj)./sum(plant_density_f(:,:),2));
        shannon_plant_f(isnan(shannon_plant_f)) = 0; 
        
        shannon_animal_t(:,jj) = animal_density_t(:,jj)./sum(animal_density_t(:,:),2)...
                                    .*log2(animal_density_t(:,jj)./sum(animal_density_t(:,:),2));
        shannon_animal_t(isnan(shannon_animal_t)) = 0;

        shannon_plant_t(:,jj) = plant_density_t(:,jj)./sum(plant_density_t(:,:),2)...
                                    .*log2(plant_density_t(:,jj)./sum(plant_density_t(:,:),2));
        shannon_plant_t(isnan(shannon_plant_t)) = 0;
        
        simpson_animal_f(:,jj) = animal_density_f(:,jj).*(animal_density_f(:,jj)-1)./(sum(animal_density_f(:,:),2).*(sum(animal_density_f(:,:),2)-1));
        simpson_animal_t(:,jj) = animal_density_t(:,jj).*(animal_density_t(:,jj)-1)./(sum(animal_density_t(:,:),2).*(sum(animal_density_t(:,:),2)-1));
        simpson_plant_f(:,jj) = plant_density_f(:,jj).*(plant_density_f(:,jj)-1)./(sum(plant_density_f(:,:),2).*(sum(plant_density_f(:,:),2)-1));
        simpson_plant_t(:,jj) = plant_density_t(:,jj).*(plant_density_t(:,jj)-1)./(sum(plant_density_t(:,:),2).*(sum(plant_density_t(:,:),2)-1));
        
    end
    Shannon_animal_f = [Shannon_animal_f -sum(shannon_animal_f(:,:),2)];
    Shannon_plant_f = [Shannon_plant_f -sum(shannon_plant_f(:,:),2)];
    Shannon_animal_t = [Shannon_animal_t -sum(shannon_animal_t(:,:),2)];
    Shannon_plant_t = [Shannon_plant_t -sum(shannon_plant_t(:,:),2)];
    
    Simpson_animal_f = [Simpson_animal_f 1-sum(simpson_animal_f(:,:),2)];
    Simpson_animal_t = [Simpson_animal_t 1-sum(simpson_animal_t(:,:),2)];
    Simpson_plant_f = [Simpson_plant_f 1-sum(simpson_plant_f(:,:),2)];
    Simpson_plant_t = [Simpson_plant_t 1-sum(simpson_plant_t(:,:),2)];
end

time = 1:101;
subplot(1,2,1)
plot(time,sum(Shannon_animal_f,2)/size(postprocess,2),'b',time,sum(Shannon_animal_t,2)/size(postprocess,2),'r');
hold on
plot(time,sum(Simpson_animal_f,2)/size(postprocess,2),'--b',time,sum(Simpson_animal_t,2)/size(postprocess,2),'--r');
xlim([0 101])
legend('Shannon with AF','Shannon without AF','Simpson with AF','Simpson without AF','interpreter','latex')
title('Animals','interpreter','latex')
subplot(1,2,2)
plot(time,sum(Shannon_plant_f,2)/size(postprocess,2),'b',time,sum(Shannon_plant_t,2)/size(postprocess,2),'r');
hold on
plot(time,sum(Simpson_plant_f,2)/size(postprocess,2),'--b',time,sum(Simpson_plant_t,2)/size(postprocess,2),'--r');
xlim([0 101])
title('Plants','interpreter','latex')

sgtitle('Biodiversity Indices','interpreter','latex')

%%

%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONAL DIVERSITY %
% Laliberté & Legendre 2010 %
%%%%%%%%%%%%%%%%%%%%%%%%

FDis_animal_f = [];
FDis_animal_t = [];
FDis_plant_f = [];
FDis_plant_t = [];
for ii = 1:size(postprocess,2)  
    
    % Peaks
    peak_animal_f = islocalmax(Out_f{:,ii},2);
    peak_animal_t = islocalmax(Out_t{:,ii},2);
    peak_plant_f = islocalmax(postprocess{2,ii},2);
    peak_plant_t = islocalmax(postprocess{5,ii},2);
    
    % Densities
    animal_density_f = Out_f{:,ii};
    plant_density_f = postprocess{2,ii};
    animal_density_t = Out_t{:,ii};
    plant_density_t = postprocess{5,ii};
    
%     animal_density_f = Out_f{:,ii}.*peak_animal_f;
%     plant_density_f = postprocess{2,ii}.*peak_plant_f;
%     animal_density_t = Out_t{:,ii}.*peak_animal_t;
%     plant_density_t = postprocess{5,ii}.*peak_plant_t;

    % Centroid
    C_animal_f = sum(animal_density_f.*xx,2)./sum(animal_density_f,2);
    C_animal_t = sum(animal_density_t.*xx,2)./sum(animal_density_t,2);
    C_plant_f = sum(plant_density_f.*xx,2)./sum(plant_density_f,2);
    C_plant_t = sum(plant_density_t.*xx,2)./sum(plant_density_t,2);
    
    % Functional dispersion
    z_animal_f = abs(xx-C_animal_f);
    FDis_animal_f = [FDis_animal_f,sum(animal_density_f.*z_animal_f,2)./sum(animal_density_f,2)]; 
    z_animal_t = abs(xx-C_animal_t);
    FDis_animal_t = [FDis_animal_t,sum(animal_density_t.*z_animal_t,2)./sum(animal_density_t,2)];
    z_plant_f = abs(xx-C_plant_f);
    FDis_plant_f = [FDis_plant_f,sum(plant_density_f.*z_plant_f,2)./sum(plant_density_f,2)]; 
    z_plant_t = abs(xx-C_plant_t);
    FDis_plant_t = [FDis_plant_t,sum(plant_density_t.*z_plant_t,2)./sum(plant_density_t,2)]; 
end

% FDis_animal_f = mean(FDis_animal_f,2);
% FDis_animal_t = mean(FDis_animal_t,2);
% FDis_plant_f = mean(FDis_plant_f,2);
% FDis_plant_t = mean(FDis_plant_t,2);

% FDis_animal_f = var(FDis_animal_f,1,2);
% FDis_animal_t = var(FDis_animal_t,1,2);
% FDis_plant_f = var(FDis_plant_f,1,2);
% FDis_plant_t = var(FDis_plant_t,1,2);

% subplot(2,2,1)
% plot(1:101,FDis_animal_f,1:101,FDis_animal_t)
% legend('with AF evolution','without','interpreter','latex');
% xlabel('Time','interpreter','latex');
% ylabel('FDis','interpreter','latex');
% xlim([0 101])
% title('Functional dispersion for consumers peaks')
% subplot(2,2,2)
% plot(1:101,FDis_plant_f,1:101,FDis_plant_t)
% xlim([0 101])
% title('Functional dispersion for resources peaks')

% subplot(2,2,3)
% plot(1:101,FDis_animal_f,1:101,FDis_animal_t)
% legend('with AF evolution','without','interpreter','latex');
% xlabel('Time','interpreter','latex');
% ylabel('FDis','interpreter','latex');
% xlim([0 101])
% title('Functional dispersion for all consumers')
% subplot(2,2,4)
% plot(1:101,FDis_plant_f,1:101,FDis_plant_t)
% xlim([0 101])
% title('Functional dispersion for all resources')

%%
%%%%%%%%%%%
% Boxplot %

% FDis_animal_end_KC = [FDis_animal_f(end,:);FDis_animal_t(end,:)];
% FDis_plant_end_KC = [FDis_plant_f(end,:);FDis_plant_t(end,:)];
% FDis_animal_end_CK = [FDis_animal_f(end,:);FDis_animal_t(end,:)];
% FDis_plant_end_CK = [FDis_plant_f(end,:);FDis_plant_t(end,:)];
% FDis_animal_boxplot = [FDis_animal_end_KC;FDis_animal_end_CK];
% FDis_plant_boxplot = [FDis_plant_end_KC;FDis_plant_end_CK];

% average on the last 100 time steps
FDis_animal_boxplot = [mean(FDis_animal_f(end-10:end,:),1);mean(FDis_animal_t(end-10:end,:),1)];
FDis_plant_boxplot = [mean(FDis_plant_f(end-10:end,:),1);mean(FDis_plant_t(end-10:end,:),1)];

subplot(1,2,1)
boxplot(FDis_animal_boxplot')
title('Functional dispersion of consumers','interpreter','latex','FontSize',18)
xticks([1 2 3 4])
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTickLabel', {'with foraging','without'});
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)

subplot(1,2,2)
boxplot(FDis_plant_boxplot')
title('Functional dispersion of resources','interpreter','latex','FontSize',18)
xticks([1 2 3 4])
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTickLabel', {'with foraging','without'});
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)

%%

%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTIONAL DIVERSITY %
% Villéger et al 2008 %
%%%%%%%%%%%%%%%%%%%%%%%%

FDiv_animal_f = [];
FDiv_animal_t = [];
FDiv_plant_f = [];
FDiv_plant_t = [];
for ii = 1:size(postprocess,2)  
    % Peaks
    peak_animal_f = islocalmax(Out_f{:,ii},2);
    peak_animal_t = islocalmax(Out_t{:,ii},2);
    peak_plant_f = islocalmax(postprocess{2,ii},2);
    peak_plant_t = islocalmax(postprocess{5,ii},2);
    
%     Densities
    animal_density_f = Out_f{:,ii};
    plant_density_f = postprocess{2,ii};
    animal_density_t = Out_t{:,ii};
    plant_density_t = postprocess{5,ii};
    
%     animal_density_f = Out_f{:,ii}.*peak_animal_f;
%     plant_density_f = postprocess{2,ii}.*peak_plant_f;
%     animal_density_t = Out_t{:,ii}.*peak_animal_t;
%     plant_density_t = postprocess{5,ii}.*peak_plant_t;

    % Gravity center
    g_animal_f = 0;
    g_animal_t = 0;    
    g_plant_f = 0;
    g_plant_t = 0;
    
    % Animal_f
    % Euclidean distance to g center for each species
    dGi = sqrt((xx-0).^2).*(animal_density_f>0);
    % Mean distance of species to center
    S = sum(animal_density_f(:,:)>0,2);
    dG = 1./S.*sum(dGi,2);
    % Sum of abundance-weighted deviances
    weighted_abundances = animal_density_f./sum(animal_density_f,2);
    delta_d = sum(weighted_abundances.*(dGi-dG.*(weighted_abundances>0)),2);
    % Sum of absolute abundance weighted deviances
    delta_D = sum(weighted_abundances.*abs(dGi-dG.*(weighted_abundances>0)),2);

    FDiv_animal_f = [FDiv_animal_f,(delta_d+dG)./(delta_D+dG)];
    FDiv_animal_f(isnan(FDiv_animal_f))=0;
    
    % Animal_t
    dGi = sqrt((xx-0).^2).*(animal_density_t>0);
    S = sum(animal_density_t(:,:)>0,2);
    dG = 1./S.*sum(dGi,2);
    weighted_abundances = animal_density_t./sum(animal_density_t,2);
    delta_d = sum(weighted_abundances.*(dGi-dG.*(weighted_abundances>0)),2);
    delta_D = sum(weighted_abundances.*abs(dGi-dG.*(weighted_abundances>0)),2);

    FDiv_animal_t = [FDiv_animal_t,(delta_d+dG)./(delta_D+dG)];
    FDiv_animal_t(isnan(FDiv_animal_t))=0;
    
    % Plant_f
    dGi = sqrt((xx-0).^2).*(plant_density_f>0);
    S = sum(plant_density_f(:,:)>0,2);
    dG = 1./S.*sum(dGi,2);
    weighted_abundances = plant_density_f./sum(plant_density_f,2);
    delta_d = sum(weighted_abundances.*(dGi-dG.*(weighted_abundances>0)),2);
    delta_D = sum(weighted_abundances.*abs(dGi-dG.*(weighted_abundances>0)),2);

    FDiv_plant_f = [FDiv_plant_f,(delta_d+dG)./(delta_D+dG)];
    FDiv_plant_f(isnan(FDiv_plant_f))=0;
    
    % Plant_t
    dGi = sqrt((xx-0).^2).*(plant_density_t>0);
    S = sum(plant_density_t(:,:)>0,2);
    dG = 1./S.*sum(dGi,2);
    weighted_abundances = plant_density_t./sum(plant_density_t,2);
    delta_d = sum(weighted_abundances.*(dGi-dG.*(weighted_abundances>0)),2);
    delta_D = sum(weighted_abundances.*abs(dGi-dG.*(weighted_abundances>0)),2);

    FDiv_plant_t = [FDiv_plant_t,(delta_d+dG)./(delta_D+dG)];
    FDiv_plant_t(isnan(FDiv_plant_t))=0;
end

% FDiv_animal_f = mean(FDiv_animal_f,2);
% FDiv_animal_t = mean(FDiv_animal_t,2);
% FDiv_plant_f = mean(FDiv_plant_f,2);
% FDiv_plant_t = mean(FDiv_plant_t,2);

% FDiv_animal_f = median(FDiv_animal_f,2);
% FDiv_animal_t = median(FDiv_animal_t,2);
% FDiv_plant_f = median(FDiv_plant_f,2);
% FDiv_plant_t = median(FDiv_plant_t,2);

% FDiv_animal_f = var(FDiv_animal_f,1,2);
% FDiv_animal_t = var(FDiv_animal_t,1,2);
% FDiv_plant_f = var(FDiv_plant_f,1,2);
% FDiv_plant_t = var(FDiv_plant_t,1,2);

FDiv_animal_f(1) = 0;
FDiv_animal_t(1) = 0;
FDiv_plant_f(1) = 0;
FDiv_plant_t(1) = 0;

% subplot(2,2,1)
% plot(1:101,FDiv_animal_f,1:101,FDiv_animal_t)
% legend('with AF evolution','without','interpreter','latex');
% xlabel('Time','interpreter','latex');
% ylabel('FDiv','interpreter','latex');
% xlim([0 101]);
% title('Functional diversity of consumers peaks','interpreter','latex')
% subplot(2,2,2)
% plot(1:101,FDiv_plant_f,1:101,FDiv_plant_t)
% xlim([0 101]);
% title('Functional diversity of resources peaks','interpreter','latex')

% subplot(2,2,3)
% plot(1:101,FDiv_animal_f,1:101,FDiv_animal_t)
% legend('with AF evolution','without','interpreter','latex');
% xlabel('Time','interpreter','latex');
% ylabel('FDiv','interpreter','latex');
% xlim([0 101]);
% title('Functional diversity of all consumers','interpreter','latex')
% subplot(2,2,4)
% plot(1:101,FDiv_plant_f,1:101,FDiv_plant_t)
% xlim([0 101]);
% title('Functional diversity of all resources','interpreter','latex')

%%
%%%%%%%%%%%
% Boxplot %

% FDiv_animal_end_KC = [FDiv_animal_f(end,:);FDiv_animal_t(end,:)];
% FDiv_plant_end_KC = [FDiv_plant_f(end,:);FDiv_plant_t(end,:)];
% FDiv_animal_end_CK = [FDiv_animal_f(end,:);FDiv_animal_t(end,:)];
% FDiv_plant_end_CK = [FDiv_plant_f(end,:);FDiv_plant_t(end,:)];
FDiv_animal_boxplot = [FDiv_animal_end_KC;FDiv_animal_end_CK];
FDiv_plant_boxplot = [FDiv_plant_end_KC;FDiv_plant_end_CK];
 
subplot(1,2,1)
boxplot(FDiv_animal_boxplot')
title('Functional diversity of consumers','interpreter','latex','FontSize',18)
xticks([1 2 3 4])
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTickLabel', {'with foraging','without'});
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)

subplot(1,2,2)
boxplot(FDis_plant_boxplot')
title('Functional diversity of resources','interpreter','latex','FontSize',18)
xticks([1 2 3 4])
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTickLabel', {'with foraging','without'});
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)

%%

%%%%%%%%%%%%%%%%%%
% FRO : evenness %
%%%%%%%%%%%%%%%%%%

animal_density_f_end = [];
animal_density_t_end = [];
plant_density_f_end = [];
plant_density_t_end = [];

for ii = 1:size(postprocess,2)  
%     Densities
    animal_density_f = Out_f{:,ii};
    plant_density_f = postprocess{2,ii};
    animal_density_t = Out_t{:,ii};
    plant_density_t = postprocess{5,ii};
    
    % final state
%     animal_density_f_end = [animal_density_f_end;animal_density_f(end,:)];
%     animal_density_t_end = [animal_density_t_end;animal_density_t(end,:)];
%     plant_density_f_end = [plant_density_f_end;plant_density_f(end,:)];
%     plant_density_t_end = [plant_density_t_end;plant_density_t(end,:)];
    
    % average on the last 100 time step
    animal_density_f_end = [animal_density_f_end;mean(animal_density_f(end-10:end,:),1)];
    animal_density_t_end = [animal_density_t_end;mean(animal_density_t(end-10:end,:),1)];
    plant_density_f_end = [plant_density_f_end;mean(plant_density_f(end-10:end,:),1)];
    plant_density_t_end = [plant_density_t_end;mean(plant_density_t(end-10:end,:),1)];
end

% S_animal_f = sum(animal_density_f_end>0,2);
% S_animal_t = sum(animal_density_t_end>0,2);

EW_animal_f = [];
EW_animal_t = [];
EW_plant_f = [];
EW_plant_t = [];
for jj = 1:(number_of_animals-1)
    EW_animal_f = [EW_animal_f,abs(xx(jj+1) - xx(jj))./(animal_density_f_end(:,jj+1)+animal_density_f_end(:,jj))];
    EW_animal_t = [EW_animal_t,abs(xx(jj+1) - xx(jj))./(animal_density_t_end(:,jj+1)+animal_density_t_end(:,jj))];
    EW_plant_f = [EW_plant_f,abs(xx(jj+1) - xx(jj))./(plant_density_f_end(:,jj+1)+plant_density_f_end(:,jj))];
    EW_plant_t = [EW_plant_t,abs(xx(jj+1) - xx(jj))./(plant_density_t_end(:,jj+1)+plant_density_t_end(:,jj))];
end

PEW_animal_f = EW_animal_f./sum(EW_animal_f,2);
PEW_animal_t = EW_animal_t./sum(EW_animal_t,2);
PEW_plant_f = EW_plant_f./sum(EW_plant_f,2);
PEW_plant_t = EW_plant_t./sum(EW_plant_t,2);

FRO_animal_f = [];
FRO_animal_t = [];
FRO_plant_f = [];
FRO_plant_t = [];
for k = 1:200
FRO_animal_f = [FRO_animal_f;sum(min(PEW_animal_f(k,:),1/(number_of_animals-1)))];
FRO_animal_t = [FRO_animal_t;sum(min(PEW_animal_t(k,:),1/(number_of_animals-1)))];
FRO_plant_f = [FRO_plant_f;sum(min(PEW_plant_f(k,:),1/(number_of_animals-1)))];
FRO_plant_t = [FRO_plant_t;sum(min(PEW_plant_t(k,:),1/(number_of_animals-1)))];
end

% FRO_animal_KC = [FRO_animal_f,FRO_animal_t];
% FRO_plant_KC = [FRO_plant_f,FRO_plant_t];
% FRO_animal_CK = [FRO_animal_f,FRO_animal_t];
% FRO_plant_CK = [FRO_plant_f,FRO_plant_t];

FRO_boxplot_animal = [FRO_animal_f,FRO_animal_t];
FRO_boxplot_plant = [FRO_plant_f,FRO_plant_t];

subplot(1,2,1)
boxplot(FRO_boxplot_animal)
title('Functional regularity of consumers','interpreter','latex','FontSize',18)
xticks([1 2 3 4])
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTickLabel', {'with foraging','without'});
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)

subplot(1,2,2)
boxplot(FRO_boxplot_plant)
title('Functional regularity of resources','interpreter','latex','FontSize',18)
xticks([1 2 3 4])
set(gca,'TickLabelInterpreter', 'latex');
set(gca,'XTickLabel', {'with foraging','without'});
a = get(gca,'XTickLabel');
set(gca,'XTickLabel',a,'fontsize',14)


%%

%%%%%%%%%%%%%%%
% Q-DIVERSITY %
%%%%%%%%%%%%%%%

Z = abs(xx-xx');
Z = exp(-Z);
q = 2;

Qd_animal_f = [];
Qd_animal_t = [];
Qd_plant_f = [];
Qd_plant_t = [];

for ii = 1:size(postprocess,2)  
animal_density_f = Out_f{:,ii};
plant_density_f = postprocess{2,ii};
animal_density_t = Out_t{:,ii};
plant_density_t = postprocess{5,ii};

animal_density_f = animal_density_f./sum(animal_density_f,2);
animal_density_t = animal_density_t./sum(animal_density_t,2);
plant_density_f = plant_density_f./sum(plant_density_f,2);
plant_density_t = plant_density_t./sum(plant_density_t,2);

Zpi_animal_f = [];
Zpi_animal_t = [];
Zpi_plant_f = [];
Zpi_plant_t = [];
for jj = 1:number_of_animals
    Zpi_animal_f = [Zpi_animal_f,sum(Z(jj,:).*animal_density_f,2)];
    Zpi_animal_t = [Zpi_animal_t,sum(Z(jj,:).*animal_density_t,2)];
    Zpi_plant_f = [Zpi_plant_f,sum(Z(jj,:).*plant_density_f,2)];
    Zpi_plant_t = [Zpi_plant_t,sum(Z(jj,:).*plant_density_t,2)];
end
Qd_animal_f = [Qd_animal_f,sum(animal_density_f.*(Zpi_animal_f.^(q-1)),2).^(1/(1-q))];
Qd_animal_t = [Qd_animal_t,sum(animal_density_t.*(Zpi_animal_t.^(q-1)),2).^(1/(1-q))];
Qd_plant_f = [Qd_plant_f,sum(plant_density_f.*(Zpi_plant_f.^(q-1)),2).^(1/(1-q))];
Qd_plant_t = [Qd_plant_t,sum(plant_density_t.*(Zpi_plant_t.^(q-1)),2).^(1/(1-q))];
end

Qd_animal_f = mean(Qd_animal_f,2);
Qd_animal_t = mean(Qd_animal_t,2);
Qd_plant_f = mean(Qd_plant_f,2);
Qd_plant_t = mean(Qd_plant_t,2);

subplot(2,1,1)
plot(1:101,Qd_animal_f,1:101,Qd_animal_t)
subplot(2,1,2)
plot(1:101,Qd_plant_f,1:101,Qd_plant_t)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORAGING STRATEGY DIVERSIFICATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sum animal densities on niche trait
foraging_f = cell(1,size(postprocess,2));
foraging_t = cell(1,size(postprocess,2));
edge = zeros(101,1);
for ii = 1:size(postprocess,2)
    animal_niche_f = [];
    animal_f = Out_F{:,ii};
    animal_niche_t = [];
    animal_t = Out_T{:,ii};
    for jj = 1:101
    animal_niche_f = [animal_niche_f;sum(animal_f((jj-1)*number_of_foraging+1:jj*number_of_foraging,:),1)];
    animal_niche_t = [animal_niche_t;sum(animal_t((jj-1)*number_of_foraging+1:jj*number_of_foraging,:),1)];
    end
    foraging_f{:,ii} = [edge,animal_niche_f,edge];
    foraging_t{:,ii} = [edge,animal_niche_t,edge];
end

% Sum of peaks (foraging trait) of all simulations
peak_animal_f = zeros(101,number_of_animals+2);
peak_animal_t = peak_animal_f;

for ii = 1:size(Param,2)
    peak_animal_f = peak_animal_f + islocalmax(foraging_f{:,ii},2);
    peak_animal_t = peak_animal_t + islocalmax(foraging_t{:,ii},2);
end
peak_animal_f = peak_animal_f(:,2:end-1);
peak_animal_t = peak_animal_t(:,2:end-1);

% Plot of peaks dynamic (averaged on all simulations)
time = 1:101;
plot(time,sum(peak_animal_f(:,:),2)/size(postprocess,2),'b',time,sum(peak_animal_t(:,:),2)/size(postprocess,2),'r');
xlim([0 101])
title('Dynamics of foraging strategies','interpreter','latex')
ylabel('number of foraging trait peaks','interpreter','latex')
legend('with AF evolution','without','interpreter','latex');
% ylim([0 inf]);

%%
% Compute peaks for each simulation
peak_animal_f = cell(1,size(Param,2));
peak_animal_t = peak_animal_f;

for ii = 1:size(Param,2)
    peak_animal_f{:,ii} = islocalmax(foraging_f{:,ii},2);
    peak_animal_t{:,ii} = islocalmax(foraging_t{:,ii},2);
end

% Plot of peaks dynamic (for each simulation)
clf
time = 1:101;

for ii = 1:size(Param,2)
% for ii = 5
    plot(time,sum(peak_animal_f{:,ii},2),'*b',time,sum(peak_animal_t{:,ii},2),'or');
    hold on
    title('animal peaks')
end

%%
% FORAGING STRATEGY DIMORPHISM

% Percentage of dimorphism
peak_sum = [];
for ii = 1:size(Param,2)
peak_end = peak_animal_f{:,ii};
peak_end = peak_end(end,:);
peak_sum = [peak_sum,sum(peak_end,2)];
end
peak_2 = peak_sum.*(peak_sum>1);
peak_3 = peak_sum.*(peak_sum>2);
dimorphism = sum(peak_2./2)/size(postprocess,2)*100
trimorphism = sum(peak_3)/size(postprocess,2)*100
%%
% Visualization of simulations with foraging dimorphism
[I,J] = find(peak_2==2); % J = simulations with dimorphism
sub = 0;
for ii = J
    final_animal_f = postprocess{1,ii};
    final_plant_f = postprocess{2,ii};
    final_animal_t = postprocess{4,ii};
    final_plant_t = postprocess{5,ii};

    final_animal_f = reshape(final_animal_f(end,:),number_of_animals,number_of_foraging);
    final_animal_t = reshape(final_animal_t(end,:),number_of_animals,number_of_foraging);
    final_plant_f = final_plant_f(end,:);
    final_plant_t = final_plant_t(end,:);

    % PEAKS
    peaks_foraging_f = islocalmax([0,sum(final_animal_f,1),0]);
    peaks_foraging_f = peaks_foraging_f(2:end-1);

    % Visualization of animal densities matrix with only the peaks of
    % foraging strategies
    sub = sub+1;
    subplot(4,4,sub)
    imagesc(final_animal_f.*peaks_foraging_f)
    colormap pink
%     plot(foraging_trait,sum(final_animal_f,1),'b')
%     hold on
end
%%
% Visualization of parameters corresponding to simulation with dimorphism
number_simu = 1:size(J,2);
sigma = []; sigmaK = []; hmax = []; alpha_h = [];
for ii = 1:size(postprocess,2)
    param = postprocess{7,ii};
    sigma = [sigma,param(1)];
    sigmaK = [sigmaK,param(2)];
    hmax = [hmax,param(3)];
    alpha_h = [alpha_h,param(4)]; 
end

hold on
plot(0,sigma,'b*',0,sigma(J),'*r')
plot(.1,sigmaK,'b*',.1,sigmaK(J),'*r')
plot(.2,hmax,'b*',.2,hmax(J),'*r')
plot(.3,alpha_h,'b*',.3,alpha_h(J),'*r')
p = plot(0,.2,'*b',0,.2,'*r');
legend([p(1) p(2)],'monomorphism','dimorphism','interpreter','latex')
ylabel('value of parameter','interpreter','latex')

title('Parameters drawn','interpreter','latex')
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% PRODUCTIVITY %%%
% Sum of all consumers functional responses weigthed by the consumers abundances
% = flux from resources to consumers

number_of_animals = 11;
number_of_plants  = 11;
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
xmin = -5;
xmax = 5;
xx = linspace(xmin,xmax,number_of_animals);
dx = xx(2)-xx(1);
traits_of_animals = xx';

foraging_trait = linspace(0,1,number_of_foraging);
dz = foraging_trait(2) - foraging_trait(1);
Foraging_trait = reshape(foraging_trait,1,1,number_of_foraging);

% Plants traits
ymin = -5;
ymax = 5;
yy = linspace(ymin,ymax,number_of_plants);
dy =  yy(2)-yy(1);
traits_of_plants = yy';

[XX,YY] = meshgrid(xx,yy);
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

% HANDLING TIME TRADE-OFF
% alpha_h = 5; % expo
alpha_h = 1; % lineaire
% alpha_h = .5; % racine
hmin = .1;
hmax = 0.55;
h = @(z) hmin + (hmax - hmin).*z.^alpha_h;
handling_time = h(foraging_trait);
handling_time = reshape(handling_time,1,1,number_of_foraging);


Productivity_f = [];
Productivity_t = [];
for jj=1:size(Param,2)
    productivity_f = [];
    productivity_t = []; 
    
  b_animal_f = Output_a_f(:,(jj-1)*number_of_animals*number_of_foraging+1:jj*number_of_animals*number_of_foraging);
  b_plant_f = Output_p_f(:,(jj-1)*number_of_plants+1:jj*number_of_plants);
  effort_ij_f = Effort_f(:,(jj-1)*number_of_animals*number_of_foraging*number_of_plants+1:jj*number_of_animals*number_of_foraging*number_of_plants);
  b_animal_t = Output_a_t(:,(jj-1)*number_of_animals*number_of_foraging+1:jj*number_of_animals*number_of_foraging);
  b_plant_t = Output_p_t(:,(jj-1)*number_of_plants+1:jj*number_of_plants);
  effort_ij_t = Effort_t(:,(jj-1)*number_of_animals*number_of_foraging*number_of_plants+1:jj*number_of_animals*number_of_foraging*number_of_plants);

    for ii=1:101
        % AF
        B_plant_f = b_plant_f(ii,:)'; % B_plant must be a vertical vector
        B_animal_f = reshape(b_animal_f(ii,:),number_of_animals,number_of_foraging);
        % reshape effort AF
        Effort_ij_f = reshape(effort_ij_f(ii,:),number_of_animals,number_of_plants,number_of_foraging);
        % compute effort mower
        Effort_sans_of_f = zeros(number_of_animals,number_of_plants);
        Effort_sans_of_f(:,B_plant_f>seuil_effort) = Effort_sans_of_f(:,B_plant_f>seuil_effort) + B_plant_f(B_plant_f>seuil_effort)';
        Effort_sans_of_f = Effort_sans_of_f./(sum(Effort_sans_of_f,2)+(sum(Effort_sans_of_f,2)==0));
        % compute real effort
        Effort_real_f(:,:,:) = Effort_ij_f(:,:,:).*Foraging_trait + (1-Foraging_trait).*Effort_sans_of_f(:,:); 

        dc_ij_f = (1 + handling_time.*extraction_coeff.*sum((Effort_real_f.*delta_ij.*B_plant_f'),2).*dy);
        c_ij_f = (extraction_coeff.*delta_ij.*B_plant_f')./dc_ij_f;

        functional_response_animal_f = sum(Effort_real_f.*c_ij_f,2).*dy;
        functional_response_animal_f = reshape(functional_response_animal_f,number_of_animals,number_of_foraging); 
        functional_response_animal_f = functional_response_animal_f.*B_animal_f;

        productivity_f = [productivity_f,sum(functional_response_animal_f,'all')*dz*dx];
        
        % MOWER
        B_plant_t = b_plant_t(ii,:)'; % B_plant must be a vertical vector
        B_animal_t = reshape(b_animal_t(ii,:),number_of_animals,number_of_foraging);
        % reshape effort AF
        Effort_ij_t = reshape(effort_ij_t(ii,:),number_of_animals,number_of_plants,number_of_foraging);
        % compute effort mower
        Effort_sans_of_t = zeros(number_of_animals,number_of_plants);
        Effort_sans_of_t(:,B_plant_t>seuil_effort) = Effort_sans_of_t(:,B_plant_t>seuil_effort) + B_plant_t(B_plant_t>seuil_effort)';
        Effort_sans_of_t = Effort_sans_of_t./(sum(Effort_sans_of_t,2)+(sum(Effort_sans_of_t,2)==0));
        % compute real effort
        Effort_real_t(:,:,:) = Effort_ij_t(:,:,:).*Foraging_trait + (1-Foraging_trait).*Effort_sans_of_t(:,:); 

        dc_ij_t = (1 + handling_time.*extraction_coeff.*sum((Effort_real_t.*delta_ij.*B_plant_t'),2).*dy);
        c_ij_t = (extraction_coeff.*delta_ij.*B_plant_t')./dc_ij_t;

        functional_response_animal_t = sum(Effort_real_t.*c_ij_t,2).*dy;
        functional_response_animal_t = reshape(functional_response_animal_t,number_of_animals,number_of_foraging); 
        functional_response_animal_t = functional_response_animal_t.*B_animal_t;

        productivity_t = [productivity_t,sum(functional_response_animal_t,'all')*dz*dx];
    end
    Productivity_t = [Productivity_t;productivity_t];
    Productivity_f = [Productivity_f;productivity_f];
end

%% BOXPLOT 

% average productivity (on the last 100 time steps)
Productivity_boxplot_t = mean(Productivity_t(:,end-10:end),2);
Productivity_boxplot_f = mean(Productivity_f(:,end-10:end),2);
Productivity_boxplot_mean = [Productivity_boxplot_f,Productivity_boxplot_t];
boxplot(Productivity_boxplot_mean);
title('Productivity','interpreter','latex','FontSize',18)
%%
% coefficient of variation of productivity (on the last 100 time steps)
Productivity_boxplot_t = std(Productivity_t(:,end-10:end),0,2)./mean(Productivity_t(:,end-10:end),2);
Productivity_boxplot_f = std(Productivity_f(:,end-10:end),0,2)./mean(Productivity_f(:,end-10:end),2);
Productivity_boxplot_cv = [Productivity_boxplot_f,Productivity_boxplot_t];
boxplot(Productivity_boxplot_cv);

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% IMAGESC OF ONE COMMUNITY 

% simu_animal = Out_t{1,1};
% simu_plant = Output_p_t(1:50,1:11);
simu_animal = Out_f{1,1};
simu_plant = Output_p_f(1:50,1:11);

subplot(1,2,1)
% imAlpha=ones(size(simu_animal(1:50,:)));
% imAlpha(simu_animal(1:50,:)==0)=0;
% imagesc(heatmap,'AlphaData',imAlpha);
cmap = jet(256);
cmap(1,:) = 0;
imagesc(simu_animal(1:50,:));%,'AlphaData',imAlpha);
colormap(cmap);
h = colorbar;
xlabel('Niche trait $x$','interpreter','latex','FontSize',18);
set(gca,'XTick',1:11,'XTickLabel',xx)
ylabel(h,'Abundance','interpreter','latex','FontSize',18);
ylabel('Time','interpreter','latex','FontSize',18)
title('Consumers','interpreter','latex','FontSize',18)
subplot(1,2,2)
imagesc(simu_plant)
xlabel('Niche trait $y$','interpreter','latex','FontSize',18);
set(gca,'XTick',1:11,'XTickLabel',xx)
ylabel(h,'Abundance','interpreter','latex','FontSize',18);
ylabel('Time ','interpreter','latex','FontSize',18)
title('Resources','interpreter','latex','FontSize',18)
h = colorbar;

% sgtitle('Mower community','interpreter','latex','FontSize',18);
sgtitle('Foraging community','interpreter','latex','FontSize',18);

% xlabel('niche trait $x$','interpreter','latex','FontSize',20);
% set(gca,'XTick',1:11,'XTickLabel',xx)
% caxis([0 1]);
% 
% % lsline;
% % colormap(myColorMap)
% h = colorbar;
% ylabel(h,'foraging trait $z$','interpreter','latex','FontSize',20);


%% FIND PARAMETERS LEADING TO DYNAMIC COMMUNITY

biomass_animal = cell(1,length(Out_f));

for ii = 1:length(Out_f)
    biomass_animal{ii} = sum(Out_f{ii},2);
    biomass_eq = biomass_animal{ii};
    biomass_eq = biomass_eq(end-50:end);
    biomass_var(ii) = var(biomass_eq);
end

%% check the stationnary and dynamic communities
index_statio = find(biomass_var<1);
index_dyn = find(biomass_var>=1);

for ii = 1:length(index_statio)
    subplot(2,1,1)
    plot(biomass_animal{index_statio(ii)});
    hold on
end
title('total animal biomass for stationnary communities','interpreter','latex','FontSize',20)
hold off

for ii = 1:length(index_dyn)
    subplot(2,1,2)
    plot(biomass_animal{index_dyn(ii)});
    hold on
end
title('total animal biomass for dynamic communities','interpreter','latex','FontSize',20)
hold off

%% find if AF evolve in the stationnary communities (not so much, AF is stronger in dynamic communities)
close

Z_final = [];
for ii = index_statio
    final = Out_F{ii};
    final = reshape(final(end-10:end,:),number_of_animals,number_of_foraging);
    final = sum(final,1);
    z_final = final./(sum(final));
    Z_final = [Z_final;z_final];
end

Z_final_dyn = [];
for ii = index_dyn
    final = Out_F{ii};
    final = reshape(final(end-10:end,:),number_of_animals,number_of_foraging);
    final = sum(final,1);
    z_final = final./(sum(final));
    Z_final_dyn = [Z_final_dyn;z_final];
end

plot(foraging_trait,sum(Z_final,1),'-*','LineWidth',2)
hold on
plot(foraging_trait,sum(Z_final_dyn,1),'r-*','LineWidth',2)
hold off
legend('stationnary communities','dynamic communities','interpreter','latex','FontSize',15)
xlabel('foraging trait','interpreter','latex','FontSize',20)
ylabel('mean densities over all runs','interpreter','latex','FontSize',20)

set(gcf, 'Position', get(0, 'Screensize'));

%% parameters leading to stationnary or dynamic communities
% Param = [SIGMA,sigmaK,hmax,alpha_h]];
close

subplot(1,4,1)
plot(Param(1,index_statio),'b*')
hold on
plot(Param(1,index_dyn),'r*')
title('$\sigma$','interpreter','latex','FontSize',15)

subplot(1,4,2)
plot(Param(2,index_statio),'b*')
hold on
plot(Param(2,index_dyn),'r*')
title('$\sigma_K$','interpreter','latex','FontSize',15)

subplot(1,4,3)
plot(Param(3,index_statio),'b*')
hold on
plot(Param(3,index_dyn),'r*')
title('$h_{max}$','interpreter','latex','FontSize',15)

subplot(1,4,4)
plot(Param(4,index_statio),'b*')
hold on
plot(Param(4,index_dyn),'r*')
title('$\alpha_h$','interpreter','latex','FontSize',15)
legend('stationnary communities','dynamical communities','interpreter','latex','Location','northwest','FontSize',15)

set(gcf, 'Position', get(0, 'Screensize'));

%% Test of visualization of community-proxy 

% dominant trait
close

simu = Out_F{5};

dominant_x = cell(1,101);
for ii = 1:101
    matrix = simu((ii-1)*number_of_animals+1:ii*number_of_animals,:);
    sum_on_z = sum(matrix,2);
    dominant_x{ii} = find(sum_on_z==max(sum_on_z));
end

for ii = 1:length(dominant_x)
    plot(ii,xx(dominant_x{ii}),'b-*')
    hold on
end








