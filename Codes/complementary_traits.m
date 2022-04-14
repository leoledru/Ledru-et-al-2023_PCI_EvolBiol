function delta_ijk_complementary = complementary_traits(SIGMA,Tn_ik,Tm_jk)

% Prend comme argument les valeurs des traits pour l'animal i 
% et la plante j, et retourne la valeur de la fonction de complemetarite
% entre les traits pour chaque couple de traits.
% Si on considere 5 traits, Tn_ik et Tm_jk seront des vecteurs (5,1)

MU = 0;
x = Tn_ik; y = Tm_jk;
delta_ijk_complementary = 1./(SIGMA.*sqrt(2.*pi)).*exp(-(x-y-MU).^2/(2.*SIGMA^2));