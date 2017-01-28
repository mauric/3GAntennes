clear all; % effacement de toutes les variables de l’espace travail
close all; % fermeture de tous les fichiers (éventuellement) ouverts
global NOMBRE_ANTENNES; % nombre total de capteurs de l’antenne
global BINARY_DATA_RATE; % débit de la source binaire transmise
global FACTEUR_SURECH; % facteur de sur-échantillonnage au récepteur
global ROLL_OFF_FACTOR; % facteur de retombée des filtres en cosinus sur-élévé
global SAMPLING_FREQ; % fréquence d’échantillonage du signal au récepteur
global BAUD_RATE; % rapidité de modulation des données transmises
%=============================================================================
% 1- Exemple d’initialisation des ces paramètres
%=============================================================================
ROLL_OFF_FACTOR=0.3;
NOMBRE_ANTENNES=10;
FACTEUR_SURECH=2;
BANDWIDTH=200e3;
DUREE_SYMBOLE=1/BANDWIDTH;
BAUD_RATE=1/DUREE_SYMBOLE;
SAMPLING_FREQ=FACTEUR_SURECH*BAUD_RATE;
%=============================================================================
% 2- Génération des signaux sur l’ensemble des capteurs
%=============================================================================
Phis = 20*pi/180; % angle d’incidence du signal utile
Phi1 = -30*pi/180 % angle d’incidence du premier interférent
Phi2 = 60*pi/180 % angle d’incidence du second interférent
RSB = 15; % rapport de puissance (en dB) entre le signal utile
% et le bruit au niveau de chaque capteur
RSI1= 2 % rapport de puissance (en dB) entre le signal utile et l’interférent n?1
RSI2= 3 % rapport de puissance (en dB) entre le signal utile et l’interférent n?2
[MatriceR,Sig,BinaireIn,PenteSCurve]=GeneSignaux(Phis,Phi1,Phi2,RSB,RSI1,RSI2);


%il manque un argument dans l'appel de la fonction. Il faut lire :
%[MatriceIB,MatriceS,Sig,BinaireIn,PenteSCurve]=GeneSignaux(Phis,Phi1,Phi2,RSB,RSI1,RSI2)