load('prep.mat')

number_of_animals = 11;
number_of_plants  = 11;
number_of_foraging = 11;
xmin = -5;
xmax = 5;
xx = linspace(xmin,xmax,number_of_animals);
foraging_trait = linspace(0,1,number_of_foraging);

animal_t = prep{1, 1}.pre_perturb.animal_density;
animal_f = prep{1, 2}.pre_perturb.animal_density;
plant_t = prep{1, 1}.pre_perturb.plant_density;
plant_f = prep{1, 2}.pre_perturb.plant_density;
effort_t = prep{1, 1}.pre_perturb.effort_dynamic;
effort_f = prep{1, 2}.pre_perturb.effort_dynamic;

%% z_mean
z_mean = [];
for i = 0:10
    a = animal_f(end-i,:);
    a_reshape = reshape(a,11,11);
    a_sum_and_norm = sum(a_reshape,1)./sum(sum(a_reshape,1));
    z_mean = [z_mean,sum(a_sum_and_norm.*foraging_trait)];
end
Z_mean = mean(z_mean);


%% biomass
animal_biomass_t = sum(animal_t,2);
animal_biomass_t = mean(animal_biomass_t(end-10:end));
animal_biomass_f = sum(animal_f,2);
animal_biomass_f = mean(animal_biomass_f(end-10:end));
plant_biomass_t = sum(plant_t,2);
plant_biomass_t = mean(plant_biomass_t(end-10:end));
plant_biomass_f = sum(plant_f,2);
plant_biomass_f = mean(plant_biomass_f(end-10:end));

% reference dots
diff_animal = animal_biomass_f - animal_biomass_t;
diff_plant = plant_biomass_f - plant_biomass_t;


%% FDis
animal_t_x = [];
animal_f_x = [];
for i = 1:size(animal_t,1)
    a = animal_t(i,:);
    a_reshape = reshape(a,11,11);
    a_reshape = sum(a_reshape,2)';
    animal_t_x = [animal_t_x;a_reshape];
    a = animal_f(i,:);
    a_reshape = reshape(a,11,11);
    a_reshape = sum(a_reshape,2)';
    animal_f_x = [animal_f_x;a_reshape];
end

FDis_animal_f = [];
FDis_animal_t = [];
FDis_plant_f = [];
FDis_plant_t = [];
% Centroid
C_animal_f = sum(animal_f_x.*xx,2)./sum(animal_f_x,2);
C_animal_t = sum(animal_t_x.*xx,2)./sum(animal_t_x,2);
C_plant_f = sum(plant_f.*xx,2)./sum(plant_f,2);
C_plant_t = sum(plant_t.*xx,2)./sum(plant_t,2);
% Functional dispersion
z_animal_f = abs(xx-C_animal_f);
FDis_animal_f = [FDis_animal_f,sum(animal_f_x.*z_animal_f,2)./sum(animal_f_x,2)]; 
z_animal_t = abs(xx-C_animal_t);
FDis_animal_t = [FDis_animal_t,sum(animal_t_x.*z_animal_t,2)./sum(animal_t_x,2)];
z_plant_f = abs(xx-C_plant_f);
FDis_plant_f = [FDis_plant_f,sum(plant_f.*z_plant_f,2)./sum(plant_f,2)]; 
z_plant_t = abs(xx-C_plant_t);
FDis_plant_t = [FDis_plant_t,sum(plant_t.*z_plant_t,2)./sum(plant_t,2)];   

FDis_animal_f = mean(FDis_animal_f(end-10:end));
FDis_animal_t = mean(FDis_animal_t(end-10:end));
FDis_plant_f = mean(FDis_plant_f(end-10:end));
FDis_plant_t =  mean(FDis_plant_t(end-10:end));

% reference dots
diff_animal = FDis_animal_f - FDis_animal_t;
diff_plant = FDis_plant_f - FDis_plant_t;


%% productivity

% FONCTION K(y)
K0 = 50; % maximal carrying capacity
y0 = 0; % optimal trait 
sigmaK = 2.5;
Kf = @(trait_plant) K0*exp(-(trait_plant - y0).^2./(2*sigmaK^2));
% FONCTION C(y-y0)
sigmaC = sigmaK-1;
Beta = 0;
C = @(trait_plant_y,trait_plant_yi) exp(sigmaC^2*Beta^2/2).*exp(-(trait_plant_y-trait_plant_yi + sigmaC^2*Beta^2).^2./(2*sigmaC^2));
% TRAITS %
xmin = -5;
xmax = 5;
xx = linspace(xmin,xmax,number_of_animals);
dx = xx(2)-xx(1);
traits_of_animals = xx';
foraging_trait = linspace(0,1,number_of_foraging);
dz = foraging_trait(2) - foraging_trait(1);
Foraging_trait = reshape(foraging_trait,1,1,number_of_foraging);
ymin = -5;
ymax = 5;
yy = linspace(ymin,ymax,number_of_plants);
dy =  yy(2)-yy(1);
traits_of_plants = yy';
[XX,YY] = meshgrid(xx,yy);
SIGMA = .5;
delta_ij = complementary_traits(SIGMA,XX,YY);
delta_ij = delta_ij';
effort_speed_of_change = 0.5;
% Competition %
Compet = C(traits_of_plants,traits_of_plants');
% Carrying capacity %
K = Kf(traits_of_plants);
extraction_coeff = .8;   
conversion_coeff = .3;    
animal_intrinsic_growth = .1; 
animal_intraspe_compet  = .01;
plant_intrinsic_growth  = .8;
seuil_abondance = 1e-5;
seuil_effort = seuil_abondance;
% HANDLING TIME TRADE-OFF
alpha_h = 1;
hmin = .1;
hmax = 0.55;
h = @(z) hmin + (hmax - hmin).*z.^alpha_h;
handling_time = h(foraging_trait);
handling_time = reshape(handling_time,1,1,number_of_foraging);

% computation of productivity
for ii=1:size(animal_t,1)
    % AF
    B_plant_f = plant_f(ii,:)'; % B_plant must be a vertical vector
    B_animal_f = reshape(animal_f(ii,:),number_of_animals,number_of_foraging);
    % reshape effort AF
    Effort_ij_f = reshape(effort_f(ii,:),number_of_animals,number_of_plants,number_of_foraging);
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

    productivity_f(i) = sum(functional_response_animal_f,'all')*dz*dx;
    
    % MOWER
    B_plant_t = plant_t(ii,:)'; % B_plant must be a vertical vector
    B_animal_t = reshape(animal_t(ii,:),number_of_animals,number_of_foraging);
    % reshape effort AF
    Effort_ij_t = reshape(effort_t(ii,:),number_of_animals,number_of_plants,number_of_foraging);
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

    productivity_t(i) = sum(functional_response_animal_t,'all')*dz*dx;
end
Productivity_t = mean(productivity_t(end-10:end));
Productivity_f = mean(productivity_f(end-10:end));

% reference dot
diff = Productivity_f - Productivity_t;

