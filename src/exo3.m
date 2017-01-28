%% ANTENNES INTELIGENTES TP


%% VARIABLES GLOBALES


clc
clear all; % effacement de toutes les variables de lï¿½espace travail
close all; % fermeture de tous les fichiers (ï¿½ventuellement) ouverts
global NOMBRE_ANTENNES; % nombre total de capteurs de lï¿½antenne
global BINARY_DATA_RATE; % débit de la source binaire transmise
global FACTEUR_SURECH; % facteur de sur-échantillonnage au récepteur
global ROLL_OFF_FACTOR; % facteur de retombï¿½e des filtres en cosinus sur-ï¿½lï¿½vï¿½
global SAMPLING_FREQ; % fréquence déchantillonage du signal au récepteur
global BAUD_RATE; % rapidité de modulation des donniées transmises

%% INIT PARAMETRES

%=============================================================================
% 1- Exemple dï¿½initialisation des ces paramï¿½tres
%=============================================================================
ROLL_OFF_FACTOR=0.3;
NOMBRE_ANTENNES=16;
FACTEUR_SURECH=2;
BANDWIDTH=200e3;
DUREE_SYMBOLE=1/BANDWIDTH;
BAUD_RATE=1/DUREE_SYMBOLE;
SAMPLING_FREQ=FACTEUR_SURECH*BAUD_RATE;
d_sur_lambda =.5;
M = 16;

%toujours un vecteur c'est un vecteur colonne

%% ALGORITHME SIMPLE DE 
 

nombre_points = 200;    %nombre de points de la discretization
phi = linspace(-pi/2,pi/2,nombre_points); 
v = zeros(M,1);%

angle_direction = pi/5;
M_const = 1/M^0.5;
w = ones(16,1);
for m = 1:size(v,1) 
 w(m) = M_const*exp(-1i*2*pi*d_sur_lambda*(m-1)*sin(angle_direction));
 w(m) = w(m)/abs(w(m));
end  

for i = 1:size(phi,2)
    for m = 1:size(v,1) 
     v(m) = M_const*exp(-1i*2*pi*d_sur_lambda*(m-1)*sin(phi(i)));  
    end  
    C(i) = w'*v;
end
C = abs(C).^2;

%plot some data
figure(); 
plot(phi,(C),'LineWidth',1);
grid on
legend('C(phi)', 'Location', 'SouthEast');
xlabel('\phi');
ylabel('|C(\phi)|');
title('Diagramme de Rayonnement');

figure();
TracePolar(phi,(C), -50);
legend('C(phi)', 'Location', 'SouthEast');

title('Diagramme de Rayonnement');

%% GENERATION DE SIGNALES

%=============================================================================
% 2- Gï¿½nï¿½ration des signaux sur lï¿½ensemble des capteurs
%=============================================================================

Phis = pi/2;%20*pi/180; % dangle dï¿½incidence du signal utile
Phi1 = -30*pi/180 % angle dï¿½incidence du premier interfï¿½rent
Phi2 = 60*pi/180 %  angale dï¿½incidence du second interfï¿½rent
RSB = 30; % rapport de puissance (en dB) entre le signal utile
% et le bruit au niveau de chaque capteur
RSI1 = 300 % rapport de puissance (en dB) entre le signal utile et lï¿½interfï¿½rent n?1
RSI2 = 300 % rapport de puissance (en dB) entre le signal utile et lï¿½interfï¿½rent n?2

[MatriceR,Sig,BinaireIn,PenteSCurve] = GeneSignaux(Phis,Phi1,Phi2,RSB,RSI1,RSI2);


%% RAPPORT  PUISSANCE SIGNAl - BRUIT

%output
a = 3.0075;
b = 3.03525;
%input
aa = 0.73625;
bb = 0.7375;
Ps_out = a*a+b*b;
Ps_in  = aa*aa+bb*bb;

G = Ps_out/Ps_in




%% DOCUMENTATION

h = get(0,'children');
for i=length(h):-1:1
  saveas(h(i), ['exo3_' num2str(length(h)+1-i)], 'png');
end






