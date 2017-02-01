%% ANTENNES INTELIGENTES TP


%% VARIABLES GLOBALES


clc
clear all; % effacement de toutes les variables de l�espace travail
close all; % fermeture de tous les fichiers (�ventuellement) ouverts
global NOMBRE_ANTENNES; % nombre total de capteurs de l�antenne
global BINARY_DATA_RATE; % d�bit de la source binaire transmise
global FACTEUR_SURECH; % facteur de sur-�chantillonnage au r�cepteur
global ROLL_OFF_FACTOR; % facteur de retomb�e des filtres en cosinus sur-�l�v�
global SAMPLING_FREQ; % fr�quence d��chantillonage du signal au r�cepteur
global BAUD_RATE; % rapidit� de modulation des donn�es transmises

%% INIT PARAMETRES

%=============================================================================
% 1- Exemple d�initialisation des ces param�tres
%=============================================================================
ROLL_OFF_FACTOR=0.3;
NOMBRE_ANTENNES=16;
FACTEUR_SURECH=2;
BANDWIDTH=200e3;
DUREE_SYMBOLE=1/BANDWIDTH;
BAUD_RATE=1/DUREE_SYMBOLE;
SAMPLING_FREQ=FACTEUR_SURECH*BAUD_RATE;
d_sur_lambda = [.125 .25 .5];
M = 16;

%% ALGORITHME SIMPLE DE 
 

nombre_points = 1000;
phi = linspace(-pi/2,pi/2,nombre_points);
CC = zeros(3,nombre_points);
v = zeros(M,1);%

M_const = 1/sqrt(M);
w = ones(16,1)*M_const;
for j = 1:3
    for i = 1:size(phi,2)
        for m = 1:size(v,1)
         v(m) = M_const*exp(-1i*2*pi*d_sur_lambda(j)*(m-1)*sin(phi(i)));  
        end
       C(i) = abs(w'*v*v'*w);
    end
    CC(j,:) = C;
    

figure('Units', 'pixels', ...
    'Position', [100 100 500 375]);
plot(phi,(C),'LineWidth',1);
grid()
title('Initial Weights ','FontSize',12);


figure()
TracePolar(phi,(C), -50);
grid()
title('Final Weights ','FontSize',12);


end



figure()
hold on
plot(phi,CC(1,:),'-b','LineWidth',1);
plot(phi,CC(2,:),'-m','LineWidth',1);
plot(phi,CC(3,:),'-r','LineWidth',1);
hold off
grid on
title('Initial Weights ','FontSize',12);



grid on

%% GENERATION DE SIGNALES

%=============================================================================
% 2- G�n�ration des signaux sur l�ensemble des capteurs
%=============================================================================
% 
Phis = 0;%20*pi/180; % dangle d�incidence du signal utile
Phi1 = -30*pi/180 % angle d�incidence du premier interf�rent
Phi2 = 60*pi/180 %  angale d�incidence du second interf�rent
RSB = 30; % rapport de puissance (en dB) entre le signal utile
% et le bruit au niveau de chaque capteur
RSI1 = 300 % rapport de puissance (en dB) entre le signal utile et l�interf�rent n?1
RSI2 = 300 % rapport de puissance (en dB) entre le signal utile et l�interf�rent n?2


[MatriceR,MatriceS,Sig,BinaireIn,PenteSCurve] = ...
    GeneSignaux(Phis,Phi1,Phi2,RSB,RSI1,RSI2);


%code rajoute pour moi
%========================================================
%calcule de la sortie y 
M = 16;
w = ones(16,1)*1/M^0.5;
y = zeros(1,size(Sig,2));
for id =1:size(Sig,2)
    y(id) = w'* Sig(:,id);
end

%diagramme de l'oeil pour l' entree
 eyediagram(Sig(1,:),FACTEUR_SURECH) ;
% %diagramme de l'oeil pour la sortie
 eyediagram(y,FACTEUR_SURECH) ;






%% RAPPORT  PUISSANCE SIGNAl - BRUIT

% %output
% a = 3.0075;
% b = 3.03525;
% %input
% aa = 0.73625;
% bb = 0.7375;
% Ps_out = a*a+b*b;
% Ps_in  = aa*aa+bb*bb;
% 
% G = Ps_out/Ps_in
% 
% %% DOCUMENTATION
% 


h = get(0,'children');
for i=length(h):-1:1
  saveas(h(i), ['exo2_' num2str(length(h)+1-i)], 'png');
end











