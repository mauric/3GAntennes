function TracePolar(Phi,Gain,SeuilDB)
%TRACEPOLAR trace le diagramme polaire d'une fonction de gain
%   TRACEPOLAR(PHI,GAIN,SEUILDB) trace une puissance GAIN exprimée
%   sur une échelle linéaire en fonction de la valeur angulaire PHI
%   exprimée en radians.
%   La valeur SEUILDB correspond à la valeur du seuil minimal
%   qui sera tracé. Si l'amplitude maximale du gain est normalisée
%   à une valeur unité, alors une valeur de ce seuil est de -50 dB.
%
%   Date : novembre 2005 - Pascal SCALART
%

if size(Phi,1) == 1
    Phi=Phi';
end;
if size(Gain,1) == 1
    Gain=Gain'; 
end;


valmax = max(10*log10(Gain));
val = find(10*log10(Gain) < SeuilDB);Gain(val)=10^(SeuilDB/10);
valmin = SeuilDB;

Phi2=[Phi;pi+Phi(length(Phi):-1:1)];
G2=[Gain;Gain(length(Gain):-1:1)];
polar(Phi2,(10*log10(G2)-valmin)/(valmax-valmin));




