%% POSTPROCESS NICHE OVERLAP Z FIXED

clear all; clc; close all;
Color = get(gca,'colororder');

% load('../Data/effect_of_z.mat')
%% INITIALISATION %%%
number_of_animals = 11 ; % 21;
number_of_plants  = 11 ; %21;
number_of_foraging = 101; %11
nt = 10;
NjK = 3;
RHO_MEAN = zeros(NjK,number_of_foraging);
R = zeros(NjK,nt,number_of_foraging);

for jK = 1:NjK
    if (jK == 1)
        load('../Data/effect_of_z_101.mat')
    elseif(jK == 2)
        load('../Data/effect_of_z_101_K4.mat')
    elseif(jK == 3)
        load('../Data/effect_of_z_101_KGaussian_mortality_to.mat')
    end
    
    % FONCTION K(y)
    K0 = 50; % maximal carrying capacity
    y0 = 0 ; %optimal trait
    sigmaK = 2.5;
    Kf = @(trait_plant) K0*exp(-(trait_plant - y0).^2./(2*sigmaK^2));
    
    % FONCTION C(y-y0)
    Beta = 0; % Beta = 0 : symetrical competition
    sigmaC = sigmaK-1;
    C = @(trait_plant_y,trait_plant_yi) exp(sigmaC^2*Beta^2/2).*exp(-(trait_plant_y-trait_plant_yi + sigmaC^2*Beta^2).^2./(2*sigmaC^2));
    
    % TRAITS %
    % Animal traits
    xmin = -5 ;
    xmax = 5 ;
    xx = linspace(xmin,xmax,number_of_animals);
    dx = 1;
    traits_of_animals = xx';
    
    foraging_trait = linspace(0,1,number_of_foraging);
    %dz = foraging_trait(2) - foraging_trait(1);
    dz = 1;
    Foraging_trait = reshape(foraging_trait,1,1,number_of_foraging);
    
    % Plants traits
    ymin = -5 ;
    ymax = 5 ;
    yy = linspace(ymin,ymax,number_of_plants);
    dy = 1;
    traits_of_plants = yy';
    
    [YY,XX] = meshgrid(yy,xx);
    SIGMA = .9;
    delta_ij = complementary_traits(SIGMA,XX,YY);
    delta_ij = delta_ij';
    Delta_ij =repmat(delta_ij,1,number_of_foraging);
    
    % HOMOGENEOUS INITIALISATION %
    extraction_coeff = .5; %.8;         % .*ones(number_of_animals,1);
    conversion_coeff = .3;         % .*ones(number_of_animals,1);
    animal_intrinsic_growth = .1; % .*ones(number_of_animals,1);
    animal_intraspe_compet  = .01; % .*ones(number_of_animals,1);
    plant_intrinsic_growth  = .8; % .5.*ones(number_of_plants,1);
    seuil_abondance = 1e-5;
    seuil_effort = seuil_abondance;
    plant_intraspe_compet = 0.01;
    
    % HANDLING TIME TRADE-OFF
    % alpha_h = 5; % expo
    alpha_h = 1; % lineaire
    % alpha_h = .5; % racine
    hmin = .1;
    hmax = 0.55;
    h = @(z) hmin + (hmax - hmin).*z.^alpha_h;
    hz = h(foraging_trait);
    Hz = reshape(repmat(hz,number_of_animals,1),1,number_of_foraging*number_of_animals);
    zz = repmat(foraging_trait,number_of_animals,1);
    zz = zz(:)';
    %% COMPUTE RHO

    RHO = zeros(nt,number_of_foraging);
    
    for tAF = 0:nt-1
        effort = reshape(Effort(end-tAF,:),number_of_plants,number_of_animals,number_of_foraging,number_of_foraging);
        plant_density = reshape(Output_p(end-tAF,:),number_of_plants,number_of_foraging);
        Rho = zeros(1,number_of_foraging);
        parfor iAF = 1:number_of_foraging
            plant_iAF = plant_density(:,iAF);
            effort_AF  = reshape(effort(:,:,:,iAF),number_of_plants,number_of_animals*number_of_foraging);
            effort_t   = plant_iAF./sum(plant_iAF+ (sum(plant_iAF)==0) );
            effort_iAF = effort_AF.*zz+(1-zz).*effort_t;
            
            ui = effort_iAF.* plant_iAF.*Delta_ij;
            dui = 1+Hz.*extraction_coeff.*sum(ui,1);
            uui = ui./dui;
            Ui = uui'*uui;
            norm_Ui = diag(Ui);
            N_Ui = sqrt(norm_Ui+norm_Ui');
            UUi = (Ui-diag(norm_Ui))./N_Ui;
            Rho(iAF) = mean(UUi(UUi>0));
        end
        RHO(end-tAF,:) = Rho;
    end
    R(jK,:,:) = RHO;
    RHO_MEAN(jK,:) =  mean(RHO,1);
end
%% FIGURE
nr = size(R,2);
zz = repmat(foraging_trait,nr,1);
figure(1)
clf
hold on
% scatter(zz(:),R(:),'filled')
scatter(foraging_trait,RHO_MEAN(1,:),'filled')
% scatter(foraging_trait,RHO_MEAN(2,:),'filled','d')
% plot(foraging_trait,RHO_MEAN(2,:),'--','LineWidth',4,'color',Color(1,:))
plot(foraging_trait,RHO_MEAN(3,:),'-.','LineWidth',4,'color',Color(1,:))
ax = gca;
ax.FontSize = 16;
xlabel({'foraging trait $z$'},'interpreter','latex','FontSize',20)
ylabel('Niche overlap $\rho$','interpreter','latex','FontSize',20,'color','k')











