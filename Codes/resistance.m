%% Fonction métrique résistance pour perturb pulse

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out] = resistance(biomass_pre_p,biomass_post_p, biomass_pre_a, biomass_post_a, biomass_pre_tot, ...
    biomass_post_tot, productivity_pre, productivity_post, productivity_pond_pre, ...
    productivity_pond_post, FDis_pre_p, FDis_post_p, FDis_pre_a, FDis_post_a, ...
    FRO_pre_p, FRO_post_p, FRO_pre_a, FRO_post_a, t_pre, t_post)


 %idx = find(t_pre> 0.25*max(t_pre)) ; %index pour enlever l'emergence de la commuanuté dans les caluls 
idx = find(t_pre>0) ;


%%% Biomass %%%
%%%%%%%%%%%%%%%
out.resistance_tot_B = -abs(mean(biomass_pre_tot(idx)) - mean(biomass_post_tot))/mean(biomass_pre_tot(idx)) ; %community biomass 
out.resistance_a_B = -abs(mean(biomass_pre_a(idx)) - mean(biomass_post_a))/mean(biomass_pre_a(idx)) ; % animal biomass
out.resistance_p_B = -abs(mean(biomass_pre_p(idx)) - mean(biomass_post_p))/mean(biomass_pre_p(idx)) ; %plant biomass 


%%% Productivity %%%
%%%%%%%%%%%%%%%%%%%%
out.resistance_tot_P = -abs(mean(productivity_pre(idx)) - mean(productivity_post))/mean(productivity_pre(idx)) ; %community biomass 

%%% Fonctionnal diverstiy FDis %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out.resistance_a_FD = -abs(mean(FDis_pre_a(idx)) - mean(FDis_post_a))/mean(FDis_pre_a(idx)) ; % animal biomass
out.resistance_p_FD = -abs(mean(FDis_pre_p(idx)) - mean(FDis_post_p))/mean(FDis_pre_p(idx)) ; %plant biomass 


%%% FRO evenness %%%
%%%%%%%%%%%%%%%%%%%%
out.resistance_a_FRO = -abs(mean(FRO_pre_a) - mean(FRO_post_a))/mean(FRO_pre_a) ; % animal biomass
out.resistance_p_FRO = -abs(mean(FRO_pre_p) - mean(FRO_post_p))/mean(FRO_pre_p) ; %plant biomass 

end