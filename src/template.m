clear all; % effacement de toutes les variables de l�espace travail
close all; % fermeture de tous les fichiers (�ventuellement) ouverts
global NOMBRE_ANTENNES; % nombre total de capteurs de l�antenne
global BINARY_DATA_RATE; % d�bit de la source binaire transmise
global FACTEUR_SURECH; % facteur de sur-�chantillonnage au r�cepteur
global ROLL_OFF_FACTOR; % facteur de retomb�e des filtres en cosinus sur-�l�v�
global SAMPLING_FREQ; % fr�quence d��chantillonage du signal au r�cepteur
global BAUD_RATE; % rapidit� de modulation des donn�es transmises
%=============================================================================
% 1- Exemple d�initialisation des ces param�tres
%=============================================================================
ROLL_OFF_FACTOR=0.3;
NOMBRE_ANTENNES=10;
FACTEUR_SURECH=2;
BANDWIDTH=200e3;
DUREE_SYMBOLE=1/BANDWIDTH;
BAUD_RATE=1/DUREE_SYMBOLE;
SAMPLING_FREQ=FACTEUR_SURECH*BAUD_RATE;
%=============================================================================
% 2- G�n�ration des signaux sur l�ensemble des capteurs
%=============================================================================
Phis = 20*pi/180; % angle d�incidence du signal utile
Phi1 = -30*pi/180 % angle d�incidence du premier interf�rent
Phi2 = 60*pi/180 % angle d�incidence du second interf�rent
RSB = 15; % rapport de puissance (en dB) entre le signal utile
% et le bruit au niveau de chaque capteur
RSI1= 2 % rapport de puissance (en dB) entre le signal utile et l�interf�rent n?1
RSI2= 3 % rapport de puissance (en dB) entre le signal utile et l�interf�rent n?2
[MatriceR,Sig,BinaireIn,PenteSCurve]=GeneSignaux(Phis,Phi1,Phi2,RSB,RSI1,RSI2);


%il manque un argument dans l'appel de la fonction. Il faut lire :
%[MatriceIB,MatriceS,Sig,BinaireIn,PenteSCurve]=GeneSignaux(Phis,Phi1,Phi2,RSB,RSI1,RSI2)