%% Tolerance : % de variation avant/après perturb
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [out] = tolerance(biomass_pre_p,biomass_post_p, biomass_pre_a, biomass_post_a, biomass_pre_tot, ...
    biomass_post_tot, productivity_pre, productivity_post, productivity_pond_pre, ...
    productivity_pond_post, FDis_pre_p, FDis_post_p, FDis_pre_a, FDis_post_a, ...
    FRO_pre_p, FRO_post_p, FRO_pre_a, FRO_post_a, t_pre, t_post)

 %idx = find(t_pre> 0.25*max(t_pre)); %index pour enlever l'emergence de la commuanuté dans les caluls 
idx = find(t_pre>0) ;
%%

%%% Biomass %%%%
%%%%%%%%%%%%%%%%
out.biomass.tolerance_a = (biomass_post_a - mean(biomass_pre_a(idx)))/mean(biomass_pre_a(idx)) ;
out.biomass.tolerance_p = (biomass_post_p - mean(biomass_pre_p(idx)))/mean(biomass_pre_p(idx)) ;

%%% Productivity %%%
%%%%%%%%%%%%%%%%%%%%
out.productivity.tolerance_a = (productivity_post - mean(productivity_pre(idx)))/mean(productivity_pre(idx)) ;

%%% Fonctionnal diverstiy FDis %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out.FDis.tolerance_a = (FDis_post_a -  mean(FDis_pre_a(idx)))/mean(FDis_pre_a(idx)) ;
out.FDis.tolerance_p = (FDis_post_p - mean(FDis_pre_p(idx)))/mean(FDis_pre_p(idx)) ;

%%% FRO evenness %%%
%%%%%%%%%%%%%%%%%%%%
out.FRO.tolerance_a = (FRO_post_a - mean(FRO_pre_a))/mean(FRO_pre_a) ; %idx dejà fait dans le proxi
out.FRO.tolerance_p = (FRO_post_p - mean(FRO_pre_p))/mean(FRO_pre_p) ;

end