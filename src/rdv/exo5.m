%% ANTENNES INTELIGENTES TP FILTRE DE CAPON
%   Vecteur optimale :
%
%                 
%   w(m) =  inv(Rib)P(P^h inv(Rib)P)^-1 f
%        
%   Comparaison avec filtre MV
%
%  Faire le meme analyse que auparavant
%
%% VARIABLES GLOBALES
clc
clear all; % effacement de toutes les variables de lï¿½espace travail
close all; % fermeture de tous les fichiers (ï¿½ventuellement) ouverts
global NOMBRE_ANTENNES; % nombre total de capteurs de lï¿½antenne
global BINARY_DATA_RATE; % dï¿½bit de la source binaire transmise
global FACTEUR_SURECH; % facteur de sur-ï¿½chantillonnage au rï¿½cepteur
global ROLL_OFF_FACTOR; % facteur de retombï¿½e des filtres en cosinus sur-ï¿½lï¿½vï¿½
global SAMPLING_FREQ; % frï¿½quence dï¿½ï¿½chantillonage du signal au rï¿½cepteur
global BAUD_RATE; % rapiditï¿½ de modulation des donnï¿½es transmises

%% INIT PARAMETRES

%=============================================================================
% 1- Exemple dï¿½initialisation des ces paramï¿½tres
%=============================================================================
ROLL_OFF_FACTOR=0.3;
NOMBRE_ANTENNES=16;
FACTEUR_SURECH=6;
BANDWIDTH=200e3;
DUREE_SYMBOLE=1/BANDWIDTH;
BAUD_RATE=1/DUREE_SYMBOLE; 
SAMPLING_FREQ=FACTEUR_SURECH*BAUD_RATE;
d_sur_lambda =.5;
M = NOMBRE_ANTENNES;

%toujours un vecteur c'est un vecteur colonne

Phis = pi/5%20*pi/180; % dangle dï¿½incidence du signal utile
Phi1 = -30*pi/180 % angle dï¿½incidence du premier interfï¿½rent
Phi2 = +60* pi/180 %  angale dï¿½incidence du second interfï¿½rent
RSB = 10; % rapport de puissance (en dB) entre le signal utile
% et le bruit au niveau de chaque capteur
RSI1 = 800 % rapport de puissance (en dB) entre le signal utile et lï¿½interfï¿½rent n?1
RSI2 = 800 % rapport de puissance (en dB) entre le signal utile et lï¿½interfï¿½rent n?2
%% CONSTRUCTION OF VECTEUR OPTIMALE
 
nombre_points = 500;    %nombre de points de la discretization
phi = linspace(-pi/2,pi/2,nombre_points); 
v = zeros(M,1);%
vs = zeros(M,1);%
v1 = zeros(M,1);%
v2 = zeros(M,1);%

M_const = 1/sqrt(M);
w = ones(size(v,1),1);
for m = 1:size(v,1)
 v1(m) = M_const*exp(-1i*2*pi*d_sur_lambda*(m-1)*sin(Phi1));
 v2(m) = M_const*exp(-1i*2*pi*d_sur_lambda*(m-1)*sin(Phi2));
 vs(m) = M_const*exp(-1i*2*pi*d_sur_lambda*(m-1)*sin(Phis));
end  
Rib = ((M/RSI1).*v1*v1'+(M/RSI2).*v2*v2'+(1/RSB).*eye(M));
Rib_inv = inv(Rib);
P = [vs v1 v2];
f = [1 0 0]';
w =  inv(Rib)*P*inv((P'* inv(Rib)*P))*f
 
for i = 1:size(phi,2)
    for m = 1:size(v,1) 
     v(m) = M_const*exp(-1i*2*pi*d_sur_lambda*(m-1)*sin(phi(i))); %li:i complex 
    end  
    C(i) = w'*v*v'*w;
end


%% plot some data
figure(); 
plot(phi,abs(C),'LineWidth',1);
grid on
legend('C(phi)', 'Location', 'SouthEast');
xlabel('\phi');
ylabel('|C(\phi)|');
title('Diagramme de Rayonnement');

figure();
TracePolar(phi,abs(C), -50);
legend('C(phi)', 'Location', 'SouthEast');
title('Diagramme de Rayonnement');

semilogy(phi,abs(C),'-b','LineWidth',1);
grid on
legend('C(phi)', 'Location', 'SouthEast');
xlabel('log \phi');
ylabel('log |C(\phi)|');    
title('Diagramme de Rayonnement (log)');
%% GENERATION DE SIGNAUX

%=============================================================================
% 2- Generation de signaux pour l´ensemble des capteurs
%=============================================================================

[MatriceR,MatriceS,x,BinaireIn,PenteSCurve] = GeneSignaux(Phis,Phi1,Phi2,RSB,RSI1,RSI2);

y = zeros(1,size(x,2));
for id =1:size(x,2)
    y(id) = w'* x(:,id);
end

%diagramme de l'oeil pour l' entree
 eyediagram(x(1,1000:10000),FACTEUR_SURECH) ;
% %diagramme de l'oeil pour la sortie
 eyediagram(y(1000:10000),FACTEUR_SURECH);

%% RAPPORT  PUISSANCE SIGNAl - BRUIT

%output
a = 3.0075;
b = 3.03525;
%input
aa = 0.73625;
bb = 0.7375;
Ps_out = a^2+b^2;
Ps_in  = aa^2+bb^2;
G = Ps_out/Ps_in
rapport = 10*log10(Ps_out/Ps_in)

%% Estimation Rx

Rx = 1*(x*x');

for i = 1:size(phi,2)
    for m = 1:size(v,1) 
     v(m) = M_const*exp(-1i*2*pi*d_sur_lambda*(m-1)*sin(phi(i))); %li:i complex 
    end  
    sigma_y(i) = inv(v'*inv(Rx)*v);
end



figure(); 
plot(phi,abs(sigma_y),'LineWidth',1);
grid on
legend('\sigma_y(phi)', 'Location', 'SouthEast');
xlabel('\phi');
ylabel('|\sigma_y(\phi)|');
title('pseudo_spectre');


%% DOCUMENTATION

h = get(0,'children');
for i=length(h):-1:1
  saveas(h(i), ['exo4_' num2str(length(h)+1-i)], 'png');
end






