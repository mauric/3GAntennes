function    [MatriceIB,MatriceS,Sig,BinaireIn,PenteSCurve]=GeneSignaux(Phis,Phi1,Phi2,RSB,RSI1,RSI2)
% .
%
%	[Y] = DIFF_ENC(X,DELAY,INITIAL) generates the output sequence Y from X
%		such that {Y(n) = X(n) (xor) X(n-DELAY), n = 1,2,...}
%		using initial values in INITIAL for X(0),X(-1),...,X(1-DELAY).
%		If X is the input sequence with N data points, Y will be the
%		differentially encoded sequence with (N+DELAY) data points.
%		The input sequence X must be a binary sequence.
%
%	DIFF_ENC(X,DELAY) encodes the sequence using initial conditions all
%		set to 1.
%
%	DIFF_ENC(X) encodes the sequence X with DELAY=1 and with all initial
%		conditions set to 1.
%
%	See also DIFF_DEC.

%	AUTHORS : M. Zeytinoglu & N. W. Ma
%             Department of Electrical & Computer Engineering
%             Ryerson Polytechnic University
%             Toronto, Ontario, CANADA
%
%	DATE    : August 1991.
%	VERSION : 1.0

%===========================================================================
% Modifications history:
% ----------------------
%	o	Tested (and modified) under MATLAB 4.0/4.1 08.16.1993 MZ
%===========================================================================

%------------------------------------------------------------------------------
%	Check input sequence and number of input parameters
%------------------------------------------------------------------------------
if nargin == 0
    error('Requires six input arguments.')
end

global FACTEUR_SURECH;
global ROLL_OFF_FACTOR;
global SAMPLING_FREQ;
global NOMBRE_ANTENNES;
global BAUD_RATE;


%==========================================================================
% Param tres li s aux mobiles et interf rents
%==========================================================================
% simule 0.25 ms de signal ce qui g n re un nombre d' chantillons complexes
%  gal    0.25 ms *Fe
%
nbbits = 50000;
BW = 200e3;
Ts = 1/BW;
R=1/Ts;
Tb=Ts/2;
D=1/Tb;
Fe=FACTEUR_SURECH*R;

% cr ation du train binaire
BinaireIn=round(rand(nbbits,1));
trbinRe=diff_enc(BinaireIn(1:2:end),1);
trbinIm=diff_enc(BinaireIn(2:2:end),1);
% transcodage
symb =  (1-2*trbinRe(1:end-1)) + sqrt(-1)*(1-2*trbinIm(1:end-1));
nbsymb = length(symb)

% sur chantillonnage
symbsurech = zeros(FACTEUR_SURECH*nbsymb,1);
size(symb)
size(symbsurech)
symbsurech(1:FACTEUR_SURECH:end) = symb;

% filtrage
if (FACTEUR_SURECH ==1)
    filtre = 1;
    signalbdb = symbsurech;
else
    filtre = rcosine(R,Fe,'fir/sqrt',ROLL_OFF_FACTOR);
    signalbdb = filter(filtre, 1, symbsurech);
end;
nbech = length(signalbdb);

%=================================
% calcul parametres du r cepteur
%=================================
z=zeros(1,5000);
z(2500)=1+1*i;
h2=conv(filtre,filtre);
SortieFA=filter(h2,1,z);
tech=50;
TED=real(SortieFA(1+FACTEUR_SURECH/2:end-FACTEUR_SURECH/2)).*real(SortieFA(1+FACTEUR_SURECH:end)-SortieFA(1:end-FACTEUR_SURECH)); % decalage de T/2 de part et d'autre
TED=TED+imag(SortieFA(1+FACTEUR_SURECH/2:end-FACTEUR_SURECH/2)).*imag(SortieFA(1+FACTEUR_SURECH:end)-SortieFA(1:end-FACTEUR_SURECH)); % decalage de T/2 de part et d'autre
[A1,P1]=max(TED);
[A2,P2]=min(TED);
PenteSCurve=(A1-A2)/((P1-P2)/FACTEUR_SURECH)  % pente de la S-courbe P=10.64/(18Te)-> P=10.64/(18*T/K) - > PT=10.64/(18/K)
%figure;
%stem((2500+length(filtre)-40:2500+length(filtre)+40),TED(2500+length(filtre)-40:2500+length(filtre)+40))

pause (2)
% calcul des signaux sur les antennes

dphis=1;
DistCapteur=0.5; %d=lambda/4;
for antenne=1:NOMBRE_ANTENNES-1
    dphis=[dphis; exp(-j*2*pi*DistCapteur*antenne*sin(Phis))];
end;

dphi1=1;
for antenne=1:NOMBRE_ANTENNES-1
    dphi1=[dphi1; exp(-j*pi*antenne*sin(Phi1))];
end;

dphi2=1;
for antenne=1:NOMBRE_ANTENNES-1
    dphi2=[dphi2; exp(-j*pi*antenne*sin(Phi2))];
end;


% sur chantillonnage
K=20;
nbre=nbsymb/K; % BP= 1 / dur e symbole de l'interf rent = K* BP du signal utile
interfer1=randn(1,nbre)+j*randn(1,nbre);
interfer1surech = zeros(FACTEUR_SURECH*nbsymb,1);
interfer1surech(1:K*FACTEUR_SURECH:end) = interfer1;

K=5;
nbre=nbsymb/K; % BP= 1 / dur e symbole de l'interf rent = K* BP du signal utile
interfer2=randn(1,nbre)+j*randn(1,nbre);
interfer2surech = zeros(FACTEUR_SURECH*nbsymb,1);
interfer2surech(1:K*FACTEUR_SURECH:end) = interfer2;
% filtrage
filtre = rcosine(R/K,Fe,'fir/sqrt',ROLL_OFF_FACTOR);
int1 = filter(filtre, 1, interfer1surech);
int2 = filter(filtre, 1, interfer2surech);



interfer1=randn(size(signalbdb))+j*randn(size(signalbdb));
interfer2=randn(size(signalbdb))+j*randn(size(signalbdb));

% generation du bruit %
RSB=10^(RSB/10);
Psignal=mean(signalbdb.^2);
sigma2=Psignal/RSB;

% generation interferent   -30  %
RSI1=10^(RSI1/10);
PowerI1=mean(interfer1.^2);
CAG1=sqrt((Psignal/RSI1)/PowerI1);;
interfer1 = interfer1 * CAG1;

% generation interferent   +60  %
RSI2=10^(RSI2/10);
PowerI2=mean(interfer2.^2);
CAG2=sqrt((Psignal/RSI2)/PowerI2);;
interfer2 = interfer2 * CAG2;

abs(dphis)

size(interfer2)
Sig=zeros(NOMBRE_ANTENNES,length(signalbdb));
SigUtile=zeros(NOMBRE_ANTENNES,length(signalbdb));
SigBuitInterferent=zeros(NOMBRE_ANTENNES,length(signalbdb));



for antenne=1:NOMBRE_ANTENNES
    SigUtile(antenne,:)=transpose(signalbdb)*dphis(antenne);
    Sig(antenne,:)=transpose(signalbdb)*dphis(antenne)+interfer1'*dphi1(antenne)+interfer2'*dphi2(antenne)+sqrt(sigma2)*randn(1,length(signalbdb));
    SigBruitInterferent(antenne,:)=interfer1'*dphi1(antenne)+interfer2'*dphi2(antenne)+sqrt(sigma2)*randn(1,length(signalbdb));
end;

%  MatriceR=zeros(NOMBRE_ANTENNES,NOMBRE_ANTENNES);
%  MatriceR=MatriceR+sigma2*diag(ones(1,NOMBRE_ANTENNES));
%  MatriceR=MatriceR+mean(abs((signalbdb).^2))*dphis*dphis';
%  MatriceR=MatriceR+mean(abs((interfer2).^2))*dphi2*dphi2'+mean(abs((interfer1).^2))*dphi1*dphi1';
%
MatriceR=(1/length(signalbdb))*(Sig*Sig');
MatriceS=(1/length(signalbdb))*(SigUtile*SigUtile');
MatriceIB=(1/length(signalbdb))*(SigBruitInterferent*SigBruitInterferent');

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
 eyediagram(Sig(1,:),2*FACTEUR_SURECH) ;
% %diagramme de l'oeil pour la sortie
 eyediagram(y,2*FACTEUR_SURECH) ;
