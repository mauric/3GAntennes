%% ANTENNES INTELIGENTES TP
%   Vecteur optimale :
%
%               inv(Rib)*v(phi)    
%   w(m) = -------------------------
%           v(phi)'inv(Rib)v(phi)
%           
% RSB = 7; % rapport de puissance (en dB) entre le signal utile
% RSI1 = 300 % rapport de puissance (en dB) entre le signal utile et 
% RSI2 = 300
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
FACTEUR_SURECH=2;
BANDWIDTH=200e3;
DUREE_SYMBOLE=1/BANDWIDTH;
BAUD_RATE=1/DUREE_SYMBOLE; 
SAMPLING_FREQ=FACTEUR_SURECH*BAUD_RATE;
d_sur_lambda =.5;
M = 16;

%toujours un vecteur c'est un vecteur colonne

Phis = pi/5;%20*pi/180; % dangle dï¿½incidence du signal utile
Phi1 = -30*pi/180 % angle dï¿½incidence du premier interfï¿½rent
Phi2 = 60*pi/180 %  angale dï¿½incidence du second interfï¿½rent
RSB = 7; % rapport de puissance (en dB) entre le signal utile
% et le bruit au niveau de chaque capteur
RSI1 = 300 % rapport de puissance (en dB) entre le signal utile et lï¿½interfï¿½rent n?1
RSI2 = 300 % rapport de puissance (en dB) entre le signal utile et lï¿½interfï¿½rent n?2
%% CONSTRUCTION OF VECTEUR OPTIMALE
 
nombre_points = 200;    %nombre de points de la discretization
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

Rib = ((M/RSI2).*v1*v1'+(M/RSI2).*v2*v2'+(1/RSB).*eye(M));
Rib_inv = inv(Rib);
w = (Rib_inv*vs)/(vs'*Rib_inv*vs);
 
for i = 1:size(phi,2)
    for m = 1:size(v,1) 
     v(m) = M_const*exp(-1i*2*pi*d_sur_lambda*(m-1)*sin(phi(i))); %li:i complex 
    end  
    C(i) = w'*v*v'*w;
end


%% plot some data
    %% plot some data
    figure(k); 
    subplot(1,3,1)
    plot(phi,abs(C),'LineWidth',1);
    grid on
    legend('C(phi)', 'Location', 'SouthEast');
    xlabel('\phi');
    ylabel('|C(\phi)|');
    title('Diagramme de Rayonnement');
    
    subplot(1,3,2)
    TracePolar(phi,abs(C), -50);
    legend('C(phi)', 'Location', 'SouthEast');
    title('Diagramme de Rayonnement');
    
    subplot(1,3,3)
    loglog(phi,abs(C),'-b','LineWidth',1);
    grid on
    legend('C(phi)', 'Location', 'SouthEast');
    xlabel('log \phi');
    ylabel('log |C(\phi)|');    
    title('Diagramme de Rayonnement (log-log)');
%% GENERATION DE SIGNALES

%=============================================================================
% 2- Generation de signaux pour l´ensemble des capteurs
%=============================================================================

[MatriceR,MatriceS,x,BinaireIn,PenteSCurve] = GeneSignaux(Phis,Phi1,Phi2,RSB,RSI1,RSI2);

y = zeros(1,size(x,2));
for id =1:size(x,2)
    y(id) = w'* x(:,id);
end

%diagramme de l'oeil pour l' entree
 eyediagram(x(1,:),FACTEUR_SURECH) ;
% %diagramme de l'oeil pour la sortie
 eyediagram(y,FACTEUR_SURECH);

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

%% Calcul du vecteur de ponderations optimales





%% DOCUMENTATION

h = get(0,'children');
for i=length(h):-1:1
  saveas(h(i), ['exo42_' num2str(length(h)+1-i)], 'png');
end






