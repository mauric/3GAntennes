function TracePolar(Phi,Gain,SeuilDB)
%TRACEPOLAR trace le diagramme polaire d'une fonction de gain
%   TRACEPOLAR(PHI,GAIN,SEUILDB) trace une puissance GAIN exprim�e
%   sur une �chelle lin�aire en fonction de la valeur angulaire PHI
%   exprim�e en radians.
%   La valeur SEUILDB correspond � la valeur du seuil minimal
%   qui sera trac�. Si l'amplitude maximale du gain est normalis�e
%   � une valeur unit�, alors une valeur de ce seuil est de -50 dB.
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




