% Prend en entree le vecteur vertical B contenant (dans l'ordre !) les
% densites des animaux, les densites des plantes, les efforts_ij (matrice
% prealablement transformee en 1 unique vecteur vertical)
% Donne en sortie la dynamique des animaux, celles des plantes et celles
% des efforts

% Adaptation de demographic_system_evol_foraging avec une nouvelle methode
% de calcul de la dynamique des efforts (20 Janvier 2020) :

% Pour la matrice des forces d'interaction delta, et pour les efforts, les
% lignes correspondent aux plantes et les colonnes aux animaux. C'est
% l'inverse par rapport à l'ancien code.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dB = demographic_system_evol_foraging_bis(t, B, number_of_animals, delta_ij, extraction_coeff, conversion_coeff,...
       animal_intrinsic_growth, animal_intraspe_compet, plant_intrinsic_growth, seuil_effort, A_foraging,...
       handling_time, number_of_plants, A_animal, A_plant, dx, dy, dz, Compet, K, seuil_abondance,...
       effort_speed_of_change, number_of_foraging, Foraging_trait, plant_intraspe_compet)

B_animal = B(1:number_of_animals*number_of_foraging);
B_plant = B(number_of_animals*number_of_foraging + 1:number_of_animals*number_of_foraging+number_of_plants);

% B_animal est maintenant une matrice, x en ordonnees, z (foraging) en abscisse
B_animal = reshape(B_animal,number_of_animals,number_of_foraging);

% Mise a zero si inferieur au seuil d'abondance
B_animal = B_animal.*(B_animal>seuil_abondance);
B_plant = B_plant.*(B_plant>seuil_abondance);

% Efforts en matrice 3D avec z (foraging) en dimension 3
Effort_ij = B(number_of_animals*number_of_foraging+number_of_plants+1:end);
Effort_ij = reshape(Effort_ij,number_of_plants,number_of_animals,number_of_foraging);

% Efforts tondeuse
Effort_sans_of = B_plant./(sum(B_plant) + (sum(B_plant)==0)) ; %pour eviter la division par 0
%Effort_sans_of = B_plant./sum(B_plant);
Effort_sans_of = repmat(Effort_sans_of,[1 number_of_animals number_of_foraging]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
   
% Agregation des efforts AF et tondeuse pondere par valeur de trait
Effort_real = Effort_ij.*Foraging_trait + (1-Foraging_trait).*Effort_sans_of;
% Effort_real = Effort_real./(sum(Effort_real,2)+(sum(Effort_real,2)==0));

% Calcul des c_ij
dc_ij = (1 + handling_time.*extraction_coeff.*sum(Effort_real.*delta_ij.*B_plant,1).*dy);
c_ij = (extraction_coeff.*delta_ij.*B_plant)./dc_ij;    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        % Mutations des traits %
                        
diffusion_animal = A_animal*B_animal;
diffusion_foraging = B_animal*A_foraging;
diffusion_plant = A_plant*B_plant;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             % Reponses fonctionnelles des animaux et des plantes %

functional_response_animal = sum(Effort_real.*c_ij,1).*dy;                         
functional_response_animal = reshape(functional_response_animal,number_of_animals,number_of_foraging);

functional_response_plant = sum((extraction_coeff.*Effort_real.*delta_ij.*reshape(B_animal,1,number_of_animals,number_of_foraging))...
    ./dc_ij,[2 3]).*dx.*dz;                   

dB_animal = diffusion_foraging + diffusion_animal + (conversion_coeff.*functional_response_animal - animal_intrinsic_growth - animal_intraspe_compet.*sum(B_animal,'all').*dx.*dz).*B_animal;
dB_animal = reshape(dB_animal,number_of_animals*number_of_foraging,1);

% VERSION DOEBELI DIECKANN
p_eff = sum(Compet.*B_plant',2).*dy;
plant_growth = plant_intrinsic_growth*(1 - p_eff./K);
dB_plant = diffusion_plant + (plant_growth - functional_response_plant).*B_plant;

% VERSION 1.20
% dB_plant = diffusion_plant + (plant_intrinsic_growth - plant_intraspe_compet.*sum(B_plant).*dy - functional_response_plant).*B_plant;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamique des efforts foraging %

% La boucle principale se fait sur les animaux (number_of_animals). Les
% efforts sont donc mis à jour en traitant les animaux un par un et leurs
% efforts sur toutes les plantes en même temps (à l'inverse de l'ancien
% code qui a une boucle principale sur les plantes).
% La boucle interne sur le number_of_foraging permet de traiter chaque
% valeur de trait z pour chaque animal de trait y.

F_effort = zeros(number_of_animals*number_of_plants,number_of_foraging);
dc_ij = reshape(dc_ij,number_of_animals,number_of_foraging); % matrice : number_of_animals rows / number_of_foraging cols
index_j = 1:number_of_animals;
Ij = logical( (sum(B_animal,2)>seuil_abondance)  );
index_j = index_j( Ij  );
for j = index_j
    % On extrait les efforts pour les animaux de trait y_j et tous les traits z 
    effort = reshape(Effort_ij(:,j,:),number_of_plants,number_of_foraging);
    c_ij_without_denominator = extraction_coeff.*delta_ij(:,j).*B_plant;
    [Cj,Ci] = meshgrid(c_ij_without_denominator,c_ij_without_denominator);
    C = ((B_plant)').*max(Cj-Ci,0);
    diagQ = diag(-sum(C,2));
    Q_phi = (C + diagQ)';
    
    index_z = 1:number_of_foraging;
    Ij = logical( (B_animal(j,:)>seuil_abondance)  );
    index_z = index_z( Ij  );
    effort_dyn = zeros(number_of_plants,number_of_foraging);
    for z = index_z
        % La matrice de diffusion est adaptée pour chaque trait z
        Q_phi_z = Q_phi.*(1./dc_ij(j,z));
        effort_dyn(:,z) = effort_speed_of_change.*B_animal(j,z).*Q_phi_z*effort(:,z);
    end
    F_effort(number_of_plants*(j-1)+1:number_of_plants*j,:) = effort_dyn;
end
deffort = reshape(F_effort,[],1);

dB = [dB_animal;dB_plant;deffort];
end
