%% POSTPROCESS POUR LES FIGURES DE L'EFFET SUR Z SUR BIOMASSE, DIV, PRODUCTIVITY

clear all; clc;

% load('effect_of_z.mat')
load('effect_of_z_101.mat')

number_of_animals = 11;
number_of_plants  = 11;
% number_of_foraging = 11;
number_of_foraging = 101;

xmin = -5;
xmax = 5;
xx = linspace(xmin,xmax,number_of_animals);
dx = 1;
traits_of_animals = xx';
foraging_trait = linspace(0,1,number_of_foraging);

for ii = 1:length(foraging_trait)
  postprocess{1,ii} = Output_a(:,(ii-1)*number_of_animals*number_of_foraging+1:ii*number_of_animals*number_of_foraging);
  postprocess{2,ii} = Output_p(:,(ii-1)*number_of_plants+1:ii*number_of_plants);
end

%% Reshape animal density
Out_a = cell(1,size(postprocess,2));
for ii = 1:size(postprocess,2)
    Output_a = [];
    output_a = postprocess{1,ii};
    for jj = 1:size(output_a,1)
        Output_a = [Output_a;reshape(output_a(jj,:),number_of_animals,number_of_foraging)];  
    end
    Out_A{1,ii} = Output_a;
end

% Sum animal densities on foraging trait
animal_niche = [];
for ii = 1:size(postprocess,2)
    animal_niche = [];
    animal = Out_A{:,ii};
    for jj = 1:101
    animal_niche = [animal_niche;sum(animal((jj-1)*number_of_animals+1:jj*number_of_animals,:),2)'];
    end
    Out_a{:,ii} = animal_niche;
end

%% BIOMASS
%% TEMPORAL DYNAMIC
biomass_animal = [];
biomass_plant = [];
for ii = 1:size(postprocess,2)
    biomass_animal = [biomass_animal;sum(postprocess{1,ii},2)'];
    biomass_plant = [biomass_plant;sum(postprocess{2,ii},2)'];
end
subplot(1,2,1)
newDefaultColors = jet(length(foraging_trait));
set(gca, 'ColorOrder', newDefaultColors, 'NextPlot', 'replacechildren');
newColorOrder = get(gca,'ColorOrder');
plot(1:101,biomass_animal,'LineWidth',3)
xlim([0 101])
legend('$z = 0$','$z = 0.1$','$z = 0.2$','$z = 0.3$','$z = 0.4$','$z = 0.5$',...
    '$z = 0.6$','$z = 0.7$','$z = 0.8$','$z = 0.9$','$z = 1$',...
    'interpreter','latex','FontSize',15,'Location','southeast')
xlabel('time','interpreter','latex','FontSize',20);
ylabel('biomass of consumers','interpreter','latex','FontSize',20)
subplot(1,2,2)
newDefaultColors = jet(length(foraging_trait));
set(gca, 'ColorOrder', newDefaultColors, 'NextPlot', 'replacechildren');
newColorOrder = get(gca,'ColorOrder');
plot(1:101,biomass_plant,'LineWidth',3)
xlim([0 101])
legend('$z = 0$','$z = 0.1$','$z = 0.2$','$z = 0.3$','$z = 0.4$','$z = 0.5$',...
    '$z = 0.6$','$z = 0.7$','$z = 0.8$','$z = 0.9$','$z = 1$',...
    'interpreter','latex','FontSize',15,'Location','southeast')
xlabel('time','interpreter','latex','FontSize',20);
ylabel('biomass of resources','interpreter','latex','FontSize',20)
%% FINAL BIOMASS IN FUNCTION OF Z
% subplot(1,2,1)
final_biomass_a = biomass_animal(:,end-10:end); % average over last 100 steps
final_biomass_a = mean(final_biomass_a,2);
% x = foraging_trait';
% y = final_biomass_a;
% p = polyfit(x, y, 3);
% v = polyval(p, x);
% scatter(foraging_trait,final_biomass_a,50,'filled')
% hold on
% % plot(x,v,'k','LineWidth',3)
% hold off
% xlabel('foraging trait','interpreter','latex','FontSize',20)
% ylabel('consumers biomass','interpreter','latex','FontSize',20)
% 
% subplot(1,2,2)
% 
final_biomass_p = biomass_plant(:,end-10:end); % average over last 100 steps
final_biomass_p = mean(final_biomass_p,2);
% x = foraging_trait';
% y = final_biomass_p;
% p = polyfit(x, y, 3);
% v = polyval(p, x);
% scatter(foraging_trait,final_biomass_p,50,'filled')
% hold on
% % plot(x,v,'k','LineWidth',3)
% hold off
% xlabel('foraging trait','interpreter','latex','FontSize',20)
% ylabel('resources biomass','interpreter','latex','FontSize',20)

% superpose curves
yyaxis left
scatter(foraging_trait,final_biomass_a,50,'filled')
ylabel('consumers biomass','interpreter','latex','FontSize',20)
hold on 
yyaxis right
ylabel('resources biomass','interpreter','latex','FontSize',20)
scatter(foraging_trait,final_biomass_p,50,'s','filled','MarkerFaceColor',[.2 .8 .2])
legend('consumers','resources','interpreter','latex','FontSize',15,'Location','southeast')
xlabel('foraging trait','interpreter','latex','FontSize',20)
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = [.2 .8 .2];
hold off

%% FUNCTIONAL DIVERSITY %
% Laliberté & Legendre 2010 %
%%%%%%%%%%%%%%%%%%%%%%%%

FDis_animal = [];
FDis_plant = [];
for ii = 1:size(postprocess,2)      
    % Densities
    animal_density = Out_a{:,ii};
    plant_density = postprocess{2,ii};
    % Centroid
    C_animal = sum(animal_density.*xx,2)./sum(animal_density,2);
    C_plant = sum(plant_density.*xx,2)./sum(plant_density,2);  
    % Functional dispersion
    z_animal = abs(xx-C_animal);
    FDis_animal = [FDis_animal,sum(animal_density.*z_animal,2)./sum(animal_density,2)]; 
    z_plant = abs(xx-C_plant);
    FDis_plant = [FDis_plant,sum(plant_density.*z_plant,2)./sum(plant_density,2)];   
end

% average on the last 100 time steps
FDis_animal_mean = mean(FDis_animal(end-10:end,:),1);
FDis_plant_mean = mean(FDis_plant(end-10:end,:),1);

subplot(1,2,1)
x = foraging_trait;
y = FDis_animal_mean;
p = polyfit(x, y, 2);
v = polyval(p, x);
scatter(foraging_trait,FDis_animal_mean,50,'filled')
hold on
% plot(x,v,'k','LineWidth',3)
hold off
xlabel('foraging trait','interpreter','latex','FontSize',20)
ylabel('consumers FDis','interpreter','latex','FontSize',20)

subplot(1,2,2)
x = foraging_trait;
y = FDis_plant_mean;
p = polyfit(x, y, 2);

v = polyval(p, x);
scatter(foraging_trait,FDis_plant_mean,50,'filled')
hold on
% plot(x,v,'k','LineWidth',3)
hold off
xlabel('foraging trait','interpreter','latex','FontSize',20)
ylabel('resources FDis','interpreter','latex','FontSize',20)

%% superpose curves
yyaxis left
scatter(foraging_trait,FDis_animal_mean,50,'filled')
ylabel('consumers FDis','interpreter','latex','FontSize',20)
hold on 
yyaxis right
ylabel('resources FDis','interpreter','latex','FontSize',20)
scatter(foraging_trait,FDis_plant_mean,50,'s','filled','MarkerFaceColor',[.2 .8 .2])
legend('consumers','resources','interpreter','latex','FontSize',15,'Location','southeast')
xlabel('foraging trait','interpreter','latex','FontSize',20)
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = [.2 .8 .2];
hold off

%% FRO : evenness %
%%%%%%%%%%%%%%%%%%

animal_density_end = [];
plant_density_end = [];

for ii = 1:size(postprocess,2)  
%     Densities
    animal_density = Out_a{:,ii};
    plant_density = postprocess{2,ii};
       
    % average on the last 100 time step
    animal_density_end = [animal_density_end;mean(animal_density(end-10:end,:),1)];
    plant_density_end = [plant_density_end;mean(plant_density(end-10:end,:),1)];
end

EW_animal = [];
EW_plant = [];
for jj = 1:(number_of_animals-1)
    EW_animal = [EW_animal,abs(xx(jj+1) - xx(jj))./(animal_density_end(:,jj+1)+animal_density_end(:,jj))];
    EW_plant = [EW_plant,abs(xx(jj+1) - xx(jj))./(plant_density_end(:,jj+1)+plant_density_end(:,jj))];
end

PEW_animal = EW_animal./sum(EW_animal,2);
PEW_plant = EW_plant./sum(EW_plant,2);

FRO_animal = [];
FRO_plant = [];
for k = 1:size(postprocess,2) 
FRO_animal = [FRO_animal;sum(min(PEW_animal(k,:),1/(number_of_animals-1)))];
FRO_plant = [FRO_plant;sum(min(PEW_plant(k,:),1/(number_of_animals-1)))];
end

%%

% FRO_animal_mean = mean(FRO_animal(end-10:end,:),1);
% FRO_plant_mean = mean(FRO_plant(end-10:end,:),1);

subplot(1,2,1)
x = foraging_trait';
y = FRO_animal;
p = polyfit(x, y, 1);
v = polyval(p, x);
scatter(foraging_trait,FRO_animal,50,'filled')
hold on
% plot(x,v,'k','LineWidth',3)
hold off
xlabel('foraging trait','interpreter','latex','FontSize',20)
ylabel('consumers FRO','interpreter','latex','FontSize',20)

subplot(1,2,2)
x = foraging_trait';
y = FRO_plant;
p = polyfit(x, y, 3);

v = polyval(p, x);
scatter(foraging_trait,FRO_plant,50,'filled')
hold on
% plot(x,v,'k','LineWidth',3)
hold off
xlabel('foraging trait','interpreter','latex','FontSize',20)
ylabel('resources FRO','interpreter','latex','FontSize',20)

%% superpose curves
yyaxis left
scatter(foraging_trait,FRO_animal,50,'filled')
ylabel('consumers FRO','interpreter','latex','FontSize',20)
hold on 
yyaxis right
ylabel('resources FRO','interpreter','latex','FontSize',20)
scatter(foraging_trait,FRO_plant,50,'s','filled','MarkerFaceColor',[.2 .8 .2])
legend('consumers','resources','interpreter','latex','FontSize',15,'Location','southeast')
xlabel('foraging trait','interpreter','latex','FontSize',20)
ax = gca;
ax.YAxis(1).Color = 'b';
ax.YAxis(2).Color = [.2 .8 .2];
hold off

%% PRODUCTIVITY %%%
% Sum of all consumers functional responses weigthed by the consumers abundances
% = flux from resources to consumers

number_of_animals = 11;
number_of_plants  = 11;
number_of_foraging = 101;

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

[YY,XX] = meshgrid(xx,yy);
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
for jj=1:number_of_foraging
    productivity_f = [];
    productivity_t = []; 
    
  b_animal_f = Output_a(:,(jj-1)*number_of_animals*number_of_foraging+1:jj*number_of_animals*number_of_foraging);
  b_plant_f = Output_p(:,(jj-1)*number_of_plants+1:jj*number_of_plants);
  effort_ij_f = Effort(:,(jj-1)*number_of_animals*number_of_foraging*number_of_plants+1:jj*number_of_animals*number_of_foraging*number_of_plants);
%   b_animal_t = Output_a_t(:,(jj-1)*number_of_animals*number_of_foraging+1:jj*number_of_animals*number_of_foraging);
%   b_plant_t = Output_p_t(:,(jj-1)*number_of_plants+1:jj*number_of_plants);
%   effort_ij_t = Effort_t(:,(jj-1)*number_of_animals*number_of_foraging*number_of_plants+1:jj*number_of_animals*number_of_foraging*number_of_plants);

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
        
%         % MOWER
%         B_plant_t = b_plant_t(ii,:)'; % B_plant must be a vertical vector
%         B_animal_t = reshape(b_animal_t(ii,:),number_of_animals,number_of_foraging);
%         % reshape effort AF
%         Effort_ij_t = reshape(effort_ij_t(ii,:),number_of_animals,number_of_plants,number_of_foraging);
%         % compute effort mower
%         Effort_sans_of_t = zeros(number_of_animals,number_of_plants);
%         Effort_sans_of_t(:,B_plant_t>seuil_effort) = Effort_sans_of_t(:,B_plant_t>seuil_effort) + B_plant_t(B_plant_t>seuil_effort)';
%         Effort_sans_of_t = Effort_sans_of_t./(sum(Effort_sans_of_t,2)+(sum(Effort_sans_of_t,2)==0));
%         % compute real effort
%         Effort_real_t(:,:,:) = Effort_ij_t(:,:,:).*Foraging_trait + (1-Foraging_trait).*Effort_sans_of_t(:,:); 
% 
%         dc_ij_t = (1 + handling_time.*extraction_coeff.*sum((Effort_real_t.*delta_ij.*B_plant_t'),2).*dy);
%         c_ij_t = (extraction_coeff.*delta_ij.*B_plant_t')./dc_ij_t;
% 
%         functional_response_animal_t = sum(Effort_real_t.*c_ij_t,2).*dy;
%         functional_response_animal_t = reshape(functional_response_animal_t,number_of_animals,number_of_foraging); 
%         functional_response_animal_t = functional_response_animal_t.*B_animal_t;
% 
%         productivity_t = [productivity_t,sum(functional_response_animal_t,'all')*dz*dx];
    end
%     Productivity_t = [Productivity_t;productivity_t];
    Productivity_f = [Productivity_f;productivity_f];
end

%%
productivity = mean(Productivity_f(:,end-10:end),2);
scatter(foraging_trait,productivity,50,'filled')
hold on
% plot(x,v,'k','LineWidth',3)
hold off
xlabel('foraging trait','interpreter','latex','FontSize',20)
ylabel('productivity','interpreter','latex','FontSize',20)


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

%% PRODUCTIVITY IN FUNCTION OF Z_MEAN
Z_mean = [];
for ii = 1:size(postprocess,2)
    animal = postprocess{1,ii};
    z_mean = [];
    for i = 0:10
        a = animal(end-i,:);
        a_reshape = reshape(a,number_of_animals,number_of_foraging);
        a_sum_and_norm = sum(a_reshape,1)./sum(sum(a_reshape,1));
        z_mean = [z_mean,sum(a_sum_and_norm.*foraging_trait)];
    end
    Z_mean = [Z_mean, mean(z_mean)];
end

% PLOT DE LA DISTANCE ENTRE BIOMASSE_TONDEUSE ET BIOMASSE_AF EN FONCTION DU Z_MOYEN 
dist_productivity_animal = Productivity_boxplot_f - Productivity_boxplot_t; 
scatter(Z_mean,dist_productivity_animal,40,'b','filled')
lm = fitlm(Z_mean',dist_productivity_animal')
coefs = lm.Coefficients.Estimate; % 2x1 [intercept; slope]
hold on
h = refline(coefs(2),coefs(1)); % plot linear regression fit
h.Color = 'b';
h.LineWidth = 2;
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');

% legend('consumers','resources','interpreter','latex','FontSize',15,'location','southeast')
xlabel('mean foraging trait of community with AF evolution','interpreter','latex','FontSize',20)
ylabel('productivity with AF - productivity without AF','interpreter','latex','FontSize',20)
% title('Functional regularity FRO','interpreter','latex','FontSize',20)
text(.1,-2,'R-Squared = 0.037 / p-value = 8.7e-06','interpreter','latex','FontSize',15)
% text(.1,.7,'resources : R-Squared = 0.321 / p-value = 5.55e-44','interpreter','latex','FontSize',15)