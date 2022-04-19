%% CALCUL DES PROXIS 

%% %%%%%%%%%%%%%%%%%%%%% BIOMASS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
biomass_pre_p = sum(plant_density,2) ;
biomass_post_p = sum(plant_density_perturb,2) ;

biomass_pre_a = sum(animal_density,2) ;
biomass_post_a = sum(animal_density_perturb,2) ;

biomass_pre_tot = sum([biomass_pre_p biomass_pre_a],2) ; %biomasse plante + animaux pour chaque pas de temps
biomass_post_tot = sum([biomass_post_p biomass_post_a],2) ;

%% %%%%%%%%%%%%%%%%%%% PRODUCTIVITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Sum of all consumers functional responses weigthed by the consumers abundances = fluxes from resources to consumers

productivity_pre = [] ;

for i = 1: size(animal_density,1)
      
        %reshape vectors to matrix
        plant_density_prod = plant_density(i,:)'; %  must be a vertical vector
        animal_density_prod = reshape(animal_density(i,:),number_of_animals, number_of_foraging);
        effort_dynamic_prod = reshape(effort_dynamic(i,:),number_of_plants,number_of_animals,number_of_foraging);
        
        % Efforts tondeuse
        Effort_sans_af = plant_density_prod./sum(plant_density_prod);
        Effort_sans_af = repmat(Effort_sans_af,[1 number_of_animals number_of_foraging]);

        % Agregation des efforts AF et tondeuse pondere par valeur de trait
        Effort_real = effort_dynamic_prod.*Foraging_trait + (1-Foraging_trait).*Effort_sans_af; 
      
        % Calcul des c_ij
        dc_ij = (1 + handling_time.*extraction_coeff.*sum(Effort_real.*delta_ij.*plant_density_prod,1).*dy);
        c_ij = (extraction_coeff.*delta_ij.*plant_density_prod)./dc_ij;    
        
        % Reponses fonctionnelles des animaux et des plantes %
        functional_response_animal = sum(Effort_real.*c_ij,1).*dy;                         
         
        productivity_pre = [productivity_pre ; sum(functional_response_animal,'all')];
end

productivity_post = [] ;

for i = 1: size(animal_density_perturb,1)
      
        %reshape vectors to matrix
        plant_density_prod = plant_density_perturb(i,:)'; %  must be a vertical vector
        animal_density_prod = reshape(animal_density_perturb(i,:),number_of_animals, number_of_foraging);
        effort_dynamic_prod = reshape(effort_dynamic_perturb(i,:),number_of_plants,number_of_animals,number_of_foraging);
        
        % Efforts tondeuse
        Effort_sans_af = plant_density_prod./sum(plant_density_prod);
        Effort_sans_af = repmat(Effort_sans_af,[1 number_of_animals number_of_foraging]);

        % Agregation des efforts AF et tondeuse pondere par valeur de trait
        Effort_real = effort_dynamic_prod.*Foraging_trait + (1-Foraging_trait).*Effort_sans_af; 
      
        % Calcul des c_ij
        dc_ij = (1 + handling_time.*extraction_coeff.*sum(Effort_real.*delta_ij.*plant_density_prod,1).*dy);
        c_ij = (extraction_coeff.*delta_ij.*plant_density_prod)./dc_ij;    
        
        % Reponses fonctionnelles des animaux et des plantes %
        functional_response_animal = sum(Effort_real.*c_ij,1).*dy;                         
         
        productivity_post = [productivity_post;sum(functional_response_animal,'all')];
end


% Productivité pondérée par la biomasse d'animaux
productivity_pond_pre = productivity_pre./biomass_pre_a ;
productivity_pond_post = productivity_post./biomass_post_a ;

%% %%%%%%%%%%%%%%%%%%%%%%% FUNCTIONAL DIVERSITY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Laliberté & Legendre 2010 %
%%%%%%%%%%%%%%%%%%%%%%%%

% Pre-perturb

% ecraser foraging trait z
animal_density_without_z = zeros(size(animal_density,1), number_of_animals) ;

for j = 1:size(animal_density,1)
    animal_density_sans_z = reshape(animal_density(j,:),number_of_animals,number_of_foraging) ;
    animal_density_sans_z = sum(animal_density_sans_z,2) ;    
    animal_density_without_z(j,:) = animal_density_sans_z' ;
end
 
% Centroid
C_animal = sum(animal_density_without_z.*xx,2)./sum(animal_density_without_z,2);
C_plant = sum(plant_density.*xx,2)./sum(plant_density,2);
    
% Functional dispersion
z_animal = abs(xx-C_animal);
FDis_pre_a = sum(animal_density_without_z.*z_animal,2)./sum(animal_density_without_z,2); 
z_plant = abs(xx-C_plant);
FDis_pre_p = sum(plant_density.*z_plant,2)./sum(plant_density,2);     


% Post-perturb
animal_density_without_z = zeros(size(animal_density_perturb,1), number_of_animals) ;

for j = 1:size(animal_density_perturb,1)
    animal_density_sans_z = reshape(animal_density_perturb(j,:),number_of_animals,number_of_foraging) ;
    animal_density_sans_z = sum(animal_density_sans_z,2) ;    
    animal_density_without_z(j,:) = animal_density_sans_z' ;
end

% Centroid
C_animal = sum(animal_density_without_z.*xx,2)./sum(animal_density_without_z,2);
C_plant = sum(plant_density.*xx,2)./sum(plant_density,2);
    
% Functional dispersion
z_animal = abs(xx-C_animal);
FDis_post_a = sum(animal_density_without_z.*z_animal,2)./sum(animal_density_without_z,2); 
z_plant = abs(xx-C_plant);
FDis_post_p = sum(plant_density.*z_plant,2)./sum(plant_density,2);     


%% %%%%%%%%%%%%%%%%%%%%%% FRO : evenness %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%à verifier

%% Pre-p
%enlever l'emergence pre-perturb
animal_density_eq = animal_density(t_pre> 0.1*max(t_pre),:) ;
plant_density_eq = plant_density(t_pre> 0.1*max(t_pre),:) ;

% écraser foraging trait z
animal_density_without_z = zeros(size(animal_density_eq,1), number_of_animals) ;

for j = 1:size(animal_density_eq,1)
    animal_density_sans_z = reshape(animal_density_eq(j,:),number_of_animals,number_of_foraging) ;
    animal_density_sans_z = sum(animal_density_sans_z,2) ;    
    animal_density_without_z(j,:) = animal_density_sans_z' ;
end

%
xx_a = xx ;
xx_p = xx ;
FRO_pre_a = [] ;
FRO_pre_p = [] ;



for k = 1:size(animal_density_without_z,1)
    animal_density_without_zero = animal_density_without_z(k,:) ;
    plant_density_without_zero = plant_density_eq(k,:) ;
    
    if (all(animal_density_without_zero))==0 || (all(plant_density_without_zero))==0
       val_a = find(animal_density_without_zero>0) ;
       val_p = find(plant_density_without_zero>0) ;
       
       animal_density_without_zero = animal_density_without_zero(animal_density_without_zero>0) ; 
       xx_a = xx(val_a) ;
       plant_density_without_zero = plant_density_without_zero(plant_density_without_zero>0) ;
       xx_p = xx(val_p) ;
       
    end     
    
    EW_animal = [] ;
    EW_plant = [] ;
    
      
        for jj = 1:(length(xx_a)-1)
            EW_animal = [EW_animal,abs(xx_a(jj+1) - xx_a(jj))./(animal_density_without_zero(:,jj+1)+animal_density_without_zero(:,jj))];
        end     
        
        for jj = 1:(length(xx_p)-1)
            EW_plant = [EW_plant,abs(xx_p(jj+1) - xx_p(jj))./(plant_density_without_zero(:,jj+1)+plant_density_without_zero(:,jj))];
        end
    
        PEW_animal = EW_animal./sum(EW_animal,2);
        PEW_plant = EW_plant./sum(EW_plant,2);

        FRO_pre_a = [FRO_pre_a ; sum(min(PEW_animal,1/(number_of_animals-1)))];
        FRO_pre_p = [FRO_pre_p ; sum(min(PEW_plant,1/(number_of_plants-1)))];
       
end 


%% Post-perturbation

% écraser foraging trait z
animal_density_without_z = zeros(size(animal_density_perturb,1), number_of_animals) ;

for j = 1:size(animal_density_perturb,1)
    animal_density_sans_z = reshape(animal_density_perturb(j,:),number_of_animals,number_of_foraging) ;
    animal_density_sans_z = sum(animal_density_sans_z,2) ;    
    animal_density_without_z(j,:) = animal_density_sans_z' ;
end

%
xx_a = xx ;
xx_p = xx ;
FRO_post_a = [] ;
FRO_post_p = [] ;


for k = 1:size(animal_density_without_z,1)
    animal_density_without_zero = animal_density_without_z(k,:) ;
    plant_density_without_zero = plant_density_perturb(k,:) ;
    
    if (all(animal_density_without_zero))==0 || (all(plant_density_without_zero))==0
       val_a = find(animal_density_without_zero>0) ;
       val_p = find(plant_density_without_zero>0) ;
       
       animal_density_without_zero = animal_density_without_zero(animal_density_without_zero>0) ; 
       xx_a = xx(val_a) ;
       plant_density_without_zero = plant_density_without_zero(plant_density_without_zero>0) ;
       xx_p = xx(val_p) ;
       
    end     
    
    EW_animal = [] ;
    EW_plant = [] ;
    
      
        for jj = 1:(length(xx_a)-1)
            EW_animal = [EW_animal,abs(xx_a(jj+1) - xx_a(jj))./(animal_density_without_zero(:,jj+1)+animal_density_without_zero(:,jj))];
        end     
        
        
        for jj = 1:(length(xx_p)-1)  
            EW_plant = [EW_plant,abs(xx_p(jj+1) - xx_p(jj))./(plant_density_without_zero(:,jj+1)+plant_density_without_zero(:,jj))];
        end
    
        PEW_animal = EW_animal./sum(EW_animal,2);
        PEW_plant = EW_plant./sum(EW_plant,2);

        FRO_post_a = [FRO_post_a ; sum(min(PEW_animal,1/(number_of_animals-1)))];
        FRO_post_p = [FRO_post_p ; sum(min(PEW_plant,1/(number_of_plants-1)))];
       
end 