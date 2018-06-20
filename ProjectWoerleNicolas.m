clc
close all
clear all 

format short

% variables globales à utiliser dans la fonction dynamique

global K V m N P lrepos k M nu AD muk kt UNS 

% DEFINITION DES PI CORRESPONDANTS AU VARIABLES DU MODEL

PIk = 0.01 ; % PI-ressort k entre les blocs 
PIK =  0.0032 ; % PI-ressort K d'introduction de la force 
PIkt =  0.000404 ; % PI-ressort kt reliant les blocs au support 
PIs = 0.01 ; % PI-section du solide
PIpoids = 0.000016 ; % PI-charge normale appliquée sur le solide
PIV =  5 * 6.9282*1e-10 ; % PI-Vitesse appliqué à l'extremité 
                               % (non reliée au premier bloc) du ressort K 
PIdt = 0.3 ; % PI-pas de temps
PInu = sqrt(0.1) ; % PI-coefficient de viscosité 

% DEFINITION DES CARACTERISTIQUES PHYSIQUE ET GEOMETRIQUE DU SOLIDE  

% longueur du solide
L = 0.1 ; %[m]

% masse du solide
Ma = 0.0012 ; %[kg]

% module elastique (young)
E = 2.5*1e9 ; %[Pa]

% section du solide
s = PIs * L^2 ; %[m^2]

% charge normal appliqué sur le solide
poids = PIpoids * E * L^2 ; %[N]

% coefficient de frottement dynamique 
muk = 0.45 ;

% coefficient de frottement statique
mus = 0.7 ;

% DEFINITION DES CARACTERISTIQUES DU MODEL

% nombre de blocs
N = 50 ;

% masse de chaque bloc
m = Ma / N ;

% constante de rigité des ressorts entre les blocs
k = PIk * (N-1) * E * L ; 

% constante de rigidité du ressort d'introduction de la force
K = PIK * E * L ; 

% constante de rigidité des ressorts entre les blocs blocs et le support
kt = PIkt * E * L * (N-1) ; 

% vitesse de déplacement de l'exremite du premier ressort K
V = PIV * E^(0.5) * Ma^(-0.5) * L ;

% vecteur : force normal agissant sur chaque bloc 
P = (poids/N) * ones(N,1) ;

% coefficient visqueux  
nu = PInu * sqrt(k*m) ;

% distance entre les blocs au temps t = 0
dist = L / (N-1) ;

% DEFINITION DES PARAMETRES TEMPORELLES

t = 0 ; % temps initial
dt = PIdt * ( 2 * pi * sqrt( m / ( 2*k + kt ) ) ) ; % pas de temps [sec]
Tf = 2e6  * dt; % temps final [sec]
 
% definition du nombre d'iterations totale (pour mener à bien l'integration) 
I = floor(Tf/dt) ;
n = 0 ;

% definition d'un nombre d'iterations de stockage   
pasdestockage = 1 ; % toutes les [pasdestockage] itérations je stocke les valeurs 
i = 1 / pasdestockage ; 
% condition pour permettre de stocker les valeurs à partir de l'indice 1
if i == 1 
    a = 0;
else 
    a = 1;
end
% nombre d'iteration stockées 
I1 = floor(i*I);

% INITIALISATION DES VECTEURS ET MATRICES 

% vecteur de stockage du temps 
T = zeros(1,I1) ;

% vecteur des positions initiales des blocs
X0 = zeros ( N , 1 );
X0 (1:end) = linspace ( 0 , (N-1) * dist , N ) ;

% vecteur des vitesses initiales des blocs
V0 = zeros(N,1);

% vecteur de stockage position - vitesse 
XV=zeros(2*N,1);
XV(1:N)=X0;
XV(N+1:2*N)=V0;

% vecteur des variales binaire : = 1 si le ressort (kt) est attaché au bloc 
%                                = 0 si le ressort (kt) est détaché du bloc
AD = ones(N,1); 

% définition de la matrice de diagonale (-1, -2,..., -2, -1) et de sous
% diagonale (1,...,1) pour la définition des forces induites par les
% ressorts k et les forces visqueuses (nu) 
M=diag(-2*ones(N, 1),0)+diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);
M(1,1)=-1;
M(end,end)=-1;

% vecteur de stockage des variables u,n(stick) = u,n(r) - ((tau,n(r))/kt),
% avec u,n(r): la position du bloc lors de son dernier rattachement, et
% tau,n(r): la force tangentielle agissant sur le bloc lors de son dernier
% rattachement
UNS=X0;

% vecteur des longueurs au repos des ressorts
lrepos=zeros(N,1);
lrepos(1)=dist*(-1);
lrepos(N)=dist*1;

% matrice de stockage tau/p
RAP=zeros(N+2,I1);

% matrice de stockage position et vitesse
S=zeros(2*N,I1);

% matrice de stockage des AD
SAD=zeros(N,I1);

% matrice de stockage des tau
Stau=zeros(N,I1);

% matrice de stockage des contraintes
Scont=zeros(N-1,I1);

% matrice de stockage des déformations
Sdefo=zeros(N-1,I1);

% matrice de stockage des energies cinetique
Scin=zeros(1, I1);

% vecteur de stockage des energies cinetique
cin=zeros(N,1);

% matrice de stockage des energies potentielle des ressorts entre bloc
Spot=zeros(1, I1);

% vecteur de stockage des energies potentielle des ressorts entre les blocs
pot=zeros(N-1,1);

% vecteur de stockage de l'energie total du système 
etot=zeros(I1,1);

% vecteur de stockage de la longueur du système 
longueur = zeros(I1, 1);

% vecteur de stockage des énergies de rattachement des ressorts kt  
epktr = zeros(N,1);

% matrice de stockage de la somme des epktr à chaque itération 
Sepktr = zeros(1,I1);

% vecteur de stockage des énergies de détachement des ressorts kt 
epktd = zeros(N,1);

% matrice de stockage de la somme des epktd à chaque itération 
Sepktd = zeros(1,I1);

% initialisation de la valeur de l'énergie potentielle stocker dans le
% ressort K (introduction de la force) 
potK =0;

% matrice de stockage de potK à chaque itération 
SpotK = zeros(1, I1);

% initialisation du travail efféctué par les forces visqueuses
enu = 0;

% matrice de stockage de la somme du travail effectue depuis t = 0 par les
% forces visqueuses pour chaque itération 
Senu = zeros(1,I1);

% vecteur de stockage du travail effectue par les forces de frottement
% dynamique à chaque itération 
wmuk = zeros(N,1);

% matrice de stockage de la somme des travaux effectue par les forces de
% frottement dynamique depuis t = 0 pour chaque itération
Swmuk = zeros(1,I1);

% matrice de stockage des énergies potentielles dans les ressorts k entre
% les blocs
Sepkt = zeros(1, I1);

% initialisation du travail apporté au système par l'intermédiaire du
% ressort K
WKa=0;

% matrice de stockage de la somme du travail apporté au système
SW = zeros(1, I1);

% vecteur de stockage des longueur de front
FR=zeros(1,I1);

% Boucle d'itération avec la méthode de Runge-Kutta : RK4
while t<=Tf
    
    a = a + i ; 
    c = floor(a); % numéro de stockage de l'itération 
    n = n+1 ; % numéro de l'itération 
    T(c) = t;
    
%     calcul des coefficient de RK4  
       k1 = dyn(t,XV);
       k2 = dyn(t + dt/2, XV + (dt/2 * k1));
       k3 = dyn(t + dt/2, XV + (dt/2 * k2));
       k4 = dyn(t + dt, XV + (dt * k3));
       
%      vecteur de la somme des forces sur les blocs : tausomme 
       tauintersomme = dyn(t,XV);
       tausomme = (tauintersomme(N+1:2*N))*m;
       
%      vecteur de la somme des forces sur les blocs, sans les forces
%      de frottement ( kt , frottement dynamique ) : tau
       tauinter = dyntau(t,XV);
       tau = (tauinter(N+1:2*N))*m;

%      rapport: force tangentielle / force normal pour chaque bloc 
       rap = tau ./ P;
       
%      vecteur de condition de détachement Cd = 0 si détachement 
       Cd = (1-AD) | (abs(-kt*(XV(1:N)-UNS)) < mus*P);
       AD = AD & Cd; % mise à jour du vecteur binaire (attaché = 1 ; détaché = 0)
       Cdn = 1 - Cd; % composante vecteur = 1 si détachement 
       
%      vecteur de condition de ratachement
       Cr = (1-AD) & Cd & (XV(N+1:2*N) < 0.05 * V);
       AD = AD | Cr;
       
%      stockage dans UNS (position de rattachement des ressorts kt) 
       Crn = 1 - Cr;
       UNS = Cr.*((XV(1:N))-(tau/kt))+Crn.*UNS;
       
       % stockage des positions et vitesses des blocs 
       S(1:2*N,c) = XV;
          
       % stockage des rapports tau / force normale
       RAP(1:N, c) = rap;
       RAP(N+1, c) = mus; % ajout de la valeur du coefficient de frottement statique 
       RAP(N+2, c) = muk; % ajout de la valeur du coefficient de frottement dynamique 
       Stau(1:N, c) = tau ./ P;
       
%      calcul et stockage des déformations des ressort k entre les blocs
       Xplus=XV(2:N);
       Xmoins=XV(1:N-1);
       Sdefo(1:N-1, c) = ((Xplus - Xmoins) - (dist .* ones(N-1, 1)))./(dist .* ones(N-1, 1));
       
%      calcul et stockage des contraintes dans les ressorts k 
       Scont(1:N-1, c) = k .* ((Xplus - Xmoins) - dist .* ones(N-1, 1)) / s;
       
%      calcul et stockage des energies cinetique des blocs
       cin(1:N) = 0.5 * m .* ((XV(N+1:2*N)).^2);
       
%      calcul energie potentielle stockée par le ressort K (force exterieur)
       potK = 0.5 * K * (abs( (V * t) - XV(1) ))^2;
       
%      calcul du travail apporté lors de l'itération 
       WKa = WKa + K* abs((V*t - (XV(1))))*V*dt;
       
%      calcul et stockage des energies potentielle dans les ressorts k
%      entre les blocs
       pot(1:N-1) = 0.5 *  k * abs((Xplus - Xmoins) - (dist .* ones(N-1, 1)) ).^2;
       
%      calcul de l'energie potentiel qui disparait lorsque les kt se rattache 
       epktr = epktr + Cr.* ( 0.5 * kt * (abs(XV(1:N)-UNS)).^2);
       
%      calcul de l'energie potentiel qui apparaît lorsque les kt se détache 
       epktd = epktd + Cdn.* ( 0.5 * kt * (abs(XV(1:N)-UNS)).^2);
 
 %     calcul de l'energie potentielle stockés dans les ressort kt si ils sont attchés     
       epkt = AD.*(0.5 * kt * (abs(XV(1:N)-UNS)).^2);
       
%      calcul de l'energie dissipé par le frottement visqueux
       enu = (nu*(M*(XV(N+1:2*N)))) .* XV(N+1:end) * dt ; 
       
%      calcul du travail de la force de frottement dynamique 
       wmuk = wmuk + abs((ones(N,1) - AD) .* (muk*P) .* XV(N+1:end) .* dt);
       
%      calcul de la somme des energies à chaque itéraion 
       etot(c) = sum(cin) + sum(pot) - sum(epktr) + sum(epktd) + potK ...
           - WKa + sum(enu) + sum(epkt) + sum(wmuk);
       
%      stockage des valeurs calculées plus haut 
       Spot(c)=sum(pot);
       Scin(c)=sum(cin);
       Sepktr(c)=sum(epktr);
       Sepktd(c)=sum(epktd);
       SpotK(c)=potK;
       SW(c)=WKa;
       Sepkt(c)=sum(epkt);
       Senu(c) = sum(enu);
       Swmuk(c) = sum(wmuk);
       
       % stcokage de la longueur du systeme
       longueur (c) = XV(N) - XV(1); 
       
       % calcul et stockage de la position du front d'onde 
       fr = max(((1-Cd).*(X0))-X0(1));
       FR(c)=fr/L + max(1-Cd);
       
        % stockage des variables binaires (attaché / détaché)
       SAD(1:N, c) = AD;  
       
       % nouveau vecteur position - vitesse 
       XV = XV + (dt/6 * (k1 + 2*k2 + 2*k3 + k4));
       
       % nouveau temps 
       t = t+dt;
end

% AFFICHAGE DES RESULTATS

% sauvegarde et affichage des variables 
SAUV=zeros(1,15);
SAUV(1)=dt;
SAUV(2)=Tf;
SAUV(3)=N;
SAUV(4)=muk;
SAUV(5)=mus;
SAUV(6)=E;
SAUV(7)=s;
SAUV(8)=k;
SAUV(9)=K;
SAUV(10)=kt;
SAUV(11)=V;
SAUV(12)=poids;
SAUV(13)=nu;
SAUV(14)=m;
SAUV(15)=L;
f = figure;
t = uitable(f,'Data',SAUV,'Position', [100 100 1000 50]);
t.ColumnName = {'dt','T','N','muk','mus','E','s','k','K','kt','V','poids','nu','m','L'};
t.ColumnEditable = true;

% tracer de la courbe de la position du front d'onde
figure
bar(T,FR,'barwidth',300,'FaceAlpha',1)
xlabel('time [sec]')
ylabel('xf/L + 1')
legend show
title('front propagation as a function of time')

% graphique de la longueur du systeme
figure
plot(T,longueur)
xlabel('time [sec]')
ylabel('system length [m]')
title('system length as a function of time')

% graphique de l'energie totale du systeme
figure
plot(T,etot)
xlabel('time [sec]')
ylabel('system energy [J]')
title('system energy as a function of time')

% graphique de l'energie potentielle (k)
figure
plot(T,Spot)
xlabel('time [sec]')
ylabel('system potential energy [J]')
title('system potential energy as a function of time')

% graphique de l'energie cinétique des blocs
figure
plot(T,Scin)
xlabel('time [sec]')
ylabel('system kinetc energy [J]')
title('system kinetic energy as a function of time')

% graphique de l'energie de rattachement des ressorts kt
figure
plot(T,Sepktr)
xlabel('time [sec]')
ylabel('system epktr energy [J]')
title('system epktr as a function of time')

% graphique de l'energie de detachement des ressorts kt 
figure
plot(T,Sepktd)
xlabel('time [sec]')
ylabel('system epktd energy [J]')
title('system epktd energy as a function of time')

% graphique de l'energie potentielle dans les ressort kt 
figure
plot(T,Sepkt)
xlabel('time [sec]')
ylabel('system epkt energy [J]')
title('system epkt energy as a function of time')

% graphique des l'energie potentiel stockée dans le ressort K 
figure
plot(T,SpotK)
xlabel('time [sec]')
ylabel('system potK energy [J]')
title('system potK energy as a function of time')

% graphique regroupant toutes les energies du systeme ( avec leur signe ) 
figure
plot(T,-SW)
xlabel('time [sec]')
ylabel('system energy [J]')
title('system energy as a function of time')
hold on 
plot (T, etot);
plot (T, Scin);
plot (T, Spot);
plot (T, -Sepktr);
plot (T, Sepktd);
plot (T, SpotK);
plot (T,  Senu);
plot (T, Sepkt);
plot (T , Swmuk);
legend ('W','etot','cin','potk','epktr','epktd','potK','enu','epkt','wmuk');

% definition de la force dans le ressort K en fonction du numero
% d'iteration 
Ft=(K*(V*T-(S(1,1:end))));

% graphique de la vitesse des blocs en fonction du temps 
figure
plot(T,S(N+1:2*N,1:end))
xlabel('time [sec]')
ylabel('bloc speed [m/sec]')
title('bloc speed as a function of time')


% graphique du déplacement des blocs en fonction du temps
figure
plot(T,S(1:N,1:end) - X0)
xlabel('time [sec]')
ylabel('bloc dispacement [m]')
title('blocs dispacements as a function of time')

% graphique de l'energie cinetique des blocs en fonction du temps
figure
plot(T,Scin)
xlabel('time')
ylabel('bloc kinetic energy')
title('blocs kinetic energy as a function of time')

% graphique de l'energie potentiel des ressorts k en fonction du temps
figure
plot(T,Spot)
xlabel('time')
ylabel('bloc potential energy')
title('blocs potential energy as a function of time')

% graphique du rapport (force d'introduction / force normale) dans le systeme en fonction du temps 
figure()
plot(T,Ft/poids)
xlabel('time [sec]')
ylabel('Ft/Fn')
title('imposed force normalized by total normal loading as a function of time')

% graphique du rapport tau/p pour un temps donné 
figure()
plot(linspace(1,N,N),RAP(1:N,(I1)-100))
hold on
plot(linspace(1,N,N), ones(1,N)*mus)
plot(linspace(1,N,N), ones(1,N)*muk)
title('tangential force normalized by normal loading for a given time')
xlabel('bloc number')
ylabel('tangential force normalized by normal loading')

% image du rapport (contraintes tangentielles / poids) 
figure()
image(Stau', 'CDataMapping', 'scaled');
colorbar
colormap(jet)
axis xy
title('tangential force normalized by normal loading on each bloc')
xlabel('bloc number')
ylabel('iteration number')

% image des variables binaires de détachement et rattachement des ressorts
% kt 
figure()
image(SAD', 'CDataMapping', 'scaled');
colorbar
colormap(jet)
axis xy
title('SAD')
xlabel('bloc number')
ylabel('iteration number')

% graphique des déformations
figure()
image(Sdefo', 'CDataMapping', 'scaled');
colorbar
colormap(jet)
axis xy
title('deformation')
xlabel('spring number')
ylabel('iteration number')

% graphique des contraintes
figure()
image(Scont', 'CDataMapping', 'scaled');
colorbar
colormap(jet)
axis xy
title('contraintes [Pa]')
xlabel('bloc number')
ylabel('iteration number')

% FONCTION DYNAMIQUE
 
% définition de la fonction dynamique: somme des forces sur les blocs
function qd = dyn(t,XV)
global K V m muk N P M k lrepos AD nu kt UNS 
un = zeros(N,1);
un(1) = 1;
tau = (K*(V*t-(XV(1))))*un + (k*(M*(XV(1:N))+(lrepos))) + (nu*(M*(XV(N+1:2*N)))) ...
    + AD.*(-kt*(XV(1:N)-UNS)) + (ones(N,1) - AD).*(-sign(XV(N+1:2*N))*muk.*P);
qd = [XV(N+1:2*N);tau/m];
end
% définition de la fonction dynamique: somme des forces (sans les forces de
% frottements) sur les blocs
function qd = dyntau(t,XV)
global K V m N M k lrepos nu 
un = zeros(N,1);
un(1) = 1;
tau = (K*(V*t-(XV(1))))*un + (k*(M*(XV(1:N))+(lrepos))) + (nu*(M*(XV(N+1:2*N))));
qd = [XV(N+1:2*N);tau/m];
end
