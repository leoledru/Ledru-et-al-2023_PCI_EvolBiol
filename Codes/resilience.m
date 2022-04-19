%% CALCUL RESILIENCE POUR PULSE PERTURB

%ATTENTION REMETTRE IDX >0.25 quand on utilise pas le process car emergence
%sinon enlevé dans les proxis

function [out] = resilience(biomass_pre_p,biomass_post_p, biomass_pre_a, biomass_post_a, biomass_pre_tot, ...
    biomass_post_tot, productivity_pre, productivity_post, productivity_pond_pre, ...
    productivity_pond_post, FDis_pre_p, FDis_post_p, FDis_pre_a, FDis_post_a, ...
    FRO_pre_p, FRO_post_p, FRO_pre_a, FRO_post_a, t_pre_eq, t_pre, t_post)

%%%% ATTENTION il faut bien avoir le bon t_pre
 % si on utilise pas proxis_process
    %idx = find(t_pre> 0.25*max(t_pre)); %index pour enlever l'emergence de la commuanuté dans les calculs
    %t_pre_eq = t_pre(idx) ;

 % si on utilise proxis_process   
 idx = find(t_pre_eq>0) ; % idx ne sert à rien dans ce cas, cette ligne sert juste à l'annuler mais idx est utile si on utilise pas proxis process

%%% Biomass %%
%%%%%%%%%%%%%%

%%% Plant %%%
test_res_p = round(biomass_post_p,2) >= round(min(biomass_pre_p(idx)),2) & round(biomass_post_p,2) <= round(max(biomass_pre_p(idx)),2) ;

% trouver le premier indice pour lequel on ne sort plus de l'intervalle
test_res_p = flipud(test_res_p); % retourner le vecteur horizontalement pour partir de la fin

return_value_p = find(~test_res_p,1) ; %find trouve le premier "true", on cherche le premier "false", donc on prendre "non resilience"

% resilient si revenu a l'eq pendant au moins les 10 derniers "pas de temps" 
% test_res est flip dont derniers pas de temps sont en fait les 10 premiers
if sum(test_res_p(1:10))==10 && isempty(return_value_p) == 0 
    return_value_p = length(test_res_p) - return_value_p +1 +1 ; %vrai indice dans le vecteur non retourné ; +1 pour le calcul du symétrique et encore +1 car on avait la valeur avant laquelle c'était bon, on veut la première à partir de laquelle c'est bon donc +1. 
    out.resilience_p_B = t_post(return_value_p) ; %savoir a quel pas de temps post_perturb ça correspond

elseif sum(test_res_p(1:10))==10 && isempty(return_value_p) == 1
    out.resilience_p_B = "never leave equilibrum" ;

else     
        out.resilience_p_B = "not resilient";
    
end

%%% Animal %%%

test_res_a = round(biomass_post_a,2) >= round(min(biomass_pre_a(idx)),2) & round(biomass_post_a,2) <= round(max(biomass_pre_a(idx)),2) ;

% trouver le premier indice pour lequel on ne sort plus de l'intervalle
test_res_a = flipud(test_res_a); % retourner le vecteur horizontalement pour partir de la fin
return_value_a = find(~test_res_a,1) ; %find trouve le premier "true", on cherche le premier "false", donc on prendre "non resilience"

if sum(test_res_a(1:10))==10 && isempty(return_value_a) == 0 
    return_value_a = length(test_res_a) - return_value_a +1 +1 ; %vrai indice dans le vecteur non retourné ; +1 pour le calcul du symétrique et encore +1 car on avait la valeur avant laquelle c'était bon, on veut la première à partir de laquelle c'est bon donc +1. 
    out.resilience_a_B = t_post(return_value_a) ; %savoir a quel pas de temps post_perturb ça correspond

elseif sum(test_res_a(1:10))==10 && isempty(return_value_a) == 1
    out.resilience_a_B = "never leave equilibrum" ;
    
else 
    out.resilience_a_B = "not resilient";
end

%%% Community %%%

test_res_tot = round(biomass_post_tot,2) >= round(min(biomass_pre_tot(idx)),2) & round(biomass_post_tot,2) <= round(max(biomass_pre_tot(idx)),2) ;

%trouver le premier indice pour lequel on ne sort plus de l'intervalle
test_res_tot = flipud(test_res_tot); % retourner le vecteur horizontalement pour partir de la fin
return_value_tot = find(~test_res_tot,1) ; %find trouve le premier "true", on cherche le premier "false", donc on prendre "non resilience"

if sum(test_res_tot(1:10))==10 && isempty(return_value_tot) == 0 
    return_value_tot = length(test_res_tot) - return_value_tot +1 +1 ; %vrai indice dans le vecteur non retourné ; +1 pour le calcul du symétrique et encore +1 car on avait la valeur avant laquelle c'était bon, on veut la première à partir de laquelle c'est bon donc +1. 
    out.resilience_tot_B = t_post(return_value_tot) ; %savoir a quel pas de temps post_perturb ça correspond

elseif sum(test_res_tot(1:10))==10 && isempty(return_value_tot) == 1 
    out.resilience_tot_B = "never leave equilibrum" ;

else 
    out.resilience_tot_B = "not resilient";
end


%%% Productivity %%%
%%%%%%%%%%%%%%%%%%%%

%Que pour la commu car pas de productivité u système entier

%%% Community %%%

test_res = round(productivity_post,2) >= round(min(productivity_pre(idx)),2) & round(productivity_post,2) <= round(max(productivity_pre(idx)),2) ;

%trouver le premier indice pour lequel on ne sort plus de l'intervalle
test_res = flipud(test_res); % retourner le vecteur horizontalement pour partir de la fin
return_value = find(~test_res,1) ; %find trouve le premier "true", on cherche le premier "false", donc on prendre "non resilience"

% resilient si revenu a l'eq pendant au moins les 10 derniers "pas de temps" 
if sum(test_res(1:10))==10 && isempty(return_value) == 0
    return_value = length(test_res)-return_value +1 +1 ; %vrai indice dans le vecteur non retourné ; +1 pour le calcul du symétrique et encore +1 car on avait la valeur avant laquelle c'était bon, on veut la première à partir de laquelle c'est bon donc +1. 
    out.resilience_P = t_post(return_value) ;%savoir a quel pas de temps post_perturb ça correspond

elseif sum(test_res(1:10))==10 && isempty(return_value) == 1 
    out.resilience_P = "never leave equilibrum" ;

else
    out.resilience_P = "not resilient";
end

%%% Fonctionnal diversity FDis %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Plant %%%
test_res_p = round(FDis_post_p,2) >= round(min(FDis_pre_p(idx)),2) & round(FDis_post_p,2) <= round(max(FDis_pre_p(idx)),2) ;

% trouver le premier indice pour lequel on ne sort plus de l'intervalle
test_res_p = flipud(test_res_p); % retourner le vecteur horizontalement pour partir de la fin

return_value_p = find(~test_res_p,1) ; %find trouve le premier "true", on cherche le premier "false", donc on prendre "non resilience"

if sum(test_res_p(1:10))==10 && isempty(return_value_p) == 0
    return_value_p = length(test_res_p) - return_value_p +1 +1 ; %vrai indice dans le vecteur non retourné ; +1 pour le calcul du symétrique et encore +1 car on avait la valeur avant laquelle c'était bon, on veut la première à partir de laquelle c'est bon donc +1. 
    out.resilience_p_FD = t_post(return_value_p) ; %savoir a quel pas de temps post_perturb ça correspond

elseif sum(test_res_p(1:10))==10 &&isempty(return_value_p) == 1 
    out.resilience_p_FD = "never leave equilibrum" ;

else 
    out.resilience_p_FD = "not resilient";
end

%%% Animal %%%

test_res_a = round(FDis_post_a,2) >= round(min(FDis_pre_a(idx)),2) & round(FDis_post_a,2) <= round(max(FDis_pre_a(idx)),2) ;

%trouver le premier indice pour lequel on ne sort plus de l'intervalle
test_res_a = flipud(test_res_a); % retourner le vecteur horizontalement pour partir de la fin
return_value_a = find(~test_res_a,1) ; %find trouve le premier "true", on cherche le premier "false", donc on prendre "non resilience"

% resilient si revenu a l'eq pendant au moins les 10 derniers "pas de temps" 
if sum(test_res_a(1:10))==10 && isempty(return_value_a) == 0
    return_value_a = length(test_res_a) - return_value_a +1 +1 ; %vrai indice dans le vecteur non retourné ; +1 pour le calcul du symétrique et encore +1 car on avait la valeur avant laquelle c'était bon, on veut la première à partir de laquelle c'est bon donc +1. 
    out.resilience_a_FD = t_post(return_value_a) ;%savoir a quel pas de temps post_perturb ça correspond

elseif sum(test_res_a(1:10)) && isempty(return_value_a) == 1 
    out.resilience_a_FD = "never leave equilibrum" ;

else
    out.resilience_a_FD = "not resilient";
end


%%% FRO evenness %%%
%%%%%%%%%%%%%%%%%%%%

%%% Plant %%%
test_res_p = round(FRO_post_p,2) >= round(min(FRO_pre_p),2) & round(FRO_post_p,2) <= round(max(FRO_pre_p),2) ; %pas idx, emerfence déjà enlevée dans calcul du proxi

% trouver le premier indice pour lequel on ne sort plus de l'intervalle
test_res_p = flipud(test_res_p); % retourner le vecteur horizontalement pour partir de la fin

return_value_p = find(~test_res_p,1) ; %find trouve le premier "true", on cherche le premier "false", donc on prendre "non resilience"

if sum(test_res_p(1:10))==10 && isempty(return_value_p) == 0 
    return_value_p = length(test_res_p) - return_value_p +1 +1 ; %vrai indice dans le vecteur non retourné ; +1 pour le calcul du symétrique et encore +1 car on avait la valeur avant laquelle c'était bon, on veut la première à partir de laquelle c'est bon donc +1. 
    out.resilience_p_FRO = t_post(return_value_p) ; %savoir a quel pas de temps post_perturb ça correspond

elseif sum(test_res_p(1:10))==10 && isempty(return_value_p) == 1 
    out.resilience_p_FRO = "never leave equilibrum" ;

else
    out.resilience_p_FRO = "not resilient";
end

%%% Animal %%%

test_res_a = round(FRO_post_a,2) >= round(min(FRO_pre_a),2) & round(FRO_post_a,2) <= round(max(FRO_pre_a),2) ; %pas idx car emergence déjà enlevée dans calcul proxi

%trouver le premier indice pour lequel on ne sort plus de l'intervalle
test_res_a = flipud(test_res_a); % retourner le vecteur horizontalement pour partir de la fin
return_value_a = find(~test_res_a,1) ; %find trouve le premier "true", on cherche le premier "false", donc on prendre "non resilience"

if sum(test_res_a(1:10))==10 && isempty(return_value_a) ==0 
    return_value_a = length(test_res_a) - return_value_a +1 +1 ; %vrai indice dans le vecteur non retourné ; +1 pour le calcul du symétrique et encore +1 car on avait la valeur avant laquelle c'était bon, on veut la première à partir de laquelle c'est bon donc +1. 
    out.resilience_a_FRO = t_post(return_value_a) ; %savoir a quel pas de temps post_perturb ça correspond

elseif sum(test_res_a(1:10))==10 && isempty(return_value_a) == 1 
    out.resilience_a_FRO = "never leave equilibrum" ;

else
    out.resilience_a_FRO = "not resilient";
end

end