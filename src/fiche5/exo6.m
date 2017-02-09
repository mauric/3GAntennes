%% ANTENNES INTELIGENTES TP FILTRE DE CAPON
%   Vecteur optimale :
%
% ALGORITHME DE FROST
%
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
FACTEUR_SURECH=6;
BANDWIDTH=200e3;
DUREE_SYMBOLE=1/BANDWIDTH;
BAUD_RATE=1/DUREE_SYMBOLE; 
SAMPLING_FREQ=FACTEUR_SURECH*BAUD_RATE;
d_sur_lambda =.5;
M = NOMBRE_ANTENNES;

%toujours un vecteur c'est un vecteur colonne

Phis = pi/5%20*pi/180; % dangle d�incidence du signal utile
Phi1 = -30*pi/180 % angle d�incidence du premier interf�rent
Phi2 = +60* pi/180 %  angale d�incidence du second interf�rent
RSB = 10; % rapport de puissance (en dB) entre le signal utile
% et le bruit au niveau de chaque capteur
RSI1 = 800 % rapport de puissance (en dB) entre le signal utile et l�interf�rent n?1
RSI2 = 800 % rapport de puissance (en dB) entre le signal utile et l�interf�rent n?2
%% CONSTRUCTION OF VECTEUR OPTIMALE
 
nombre_points = 500;    %nombre de points de la discretization
phi = linspace(-pi/2,pi/2,nombre_points); 
v = zeros(M,1);%
vs = zeros(M,1);%
v1 = zeros(M,1);%
v2 = zeros(M,1);%

M_const = 1/sqrt(M);
w = zeros(size(v,1),1);

 m = transpose(1:size(v,1));
 v1 = M_const*exp(-1i*2*pi*d_sur_lambda*(m-1)*sin(Phi1));
 v2 = M_const*exp(-1i*2*pi*d_sur_lambda*(m-1)*sin(Phi2));
 vs = M_const*exp(-1i*2*pi*d_sur_lambda*(m-1)*sin(Phis));
 
P = [vs v1 v2];
f = [1 0 0];
wfixe = P*(P'*P)*f';
Q = eye(M) - P*inv(P'*P)*P'

%% GENERATION DE SIGNAUX

%=============================================================================
% 2- Generation de signaux pour l�ensemble des capteurs
%=============================================================================

[MatriceR,MatriceS,x,BinaireIn,PenteSCurve] = GeneSignaux(Phis,Phi1,Phi2,RSB,RSI1,RSI2);
mu= 0.003;
y = zeros(1,size(x,2));

    clear v;

for i = 1:size(phi,2)
m = 1:M ;
v(:,i) =  M_const*exp(-1i*2*pi*d_sur_lambda*(m-1)*sin(phi(i))); %li:i complex 
end


for n =1:size(x,2)
    
    y(n) = w'* x(:,n);
    w = Q*(w - 2*mu*x(:,n)*conj(y(n)))+ wfixe;
    
    if (rem(n,200)==0)
    for i = 1:size(phi,2)
     C(i) = w'*v(:,i)*v(:,i)'*w;
    end
    figure()
  %  plot(phi,10*log10(abs(C)))
    semilogy(phi,abs(C),'-b','LineWidth',1);
    grid on
    legend('C(phi)', 'Location', 'SouthEast');
    xlabel('log \phi');
    ylabel('log |C(\phi)|');    
    title('Diagramme de Rayonnement (log)');
    end
    
%tous les 200 echantillon
 %diagramme de l'oeil pour l' entree
 %eyediagram(x(1,1000:10000),FACTEUR_SURECH) ;
% %diagramme de l'oeil pour la sortie
 %eyediagram(y(1000:10000),FACTEUR_SURECH);   
    
    
end



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
rapport = 10*log10(Ps_out/Ps_in)

%% DOCUMENTATION

h = get(0,'children');
for i=length(h):-1:1
  saveas(h(i), ['exo4_' num2str(length(h)+1-i)], 'png');
end






