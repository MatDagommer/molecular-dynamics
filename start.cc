#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <assert.h>
#include "Graphics.h"  // Contains the definition of "Particule" and "Graphics"

char env[]="DISPLAY=:0";



/*┼┼┼┼┼┼┼▄▀▀▀▄▄▄▄▄▄▄▀▀▀▄┼┼┼┼┼┼
  ┼┼┼┼┼┼┼█▒▒░░░░░░░░░▒▒█┼┼┼┼┼┼
  ┼┼┼┼┼┼┼┼█░░█░░░░░█░░█┼┼┼┼┼┼┼
  ┼┼┼┼─▄▄──█░░░▀█▀░░░█──▄▄─┼┼┼
  ┼┼┼┼█░░█─▀▄░░░░░░░▄▀─█░░█┼┼┼
  ┼██░██░████░██░░░██░░░█████┼
  ┼██▄██░██▄▄░██░░░██░░░██░██┼
  ┼██▀██░██▀▀░██░░░██░░░██░██┼
  ┼██░██░████░████░████░█████┼*/



  /////////////////////////////Definition des differents types de collision/////////////////////////////////////////////////////////////////////
  enum  col_type {
    bottom,
    right,
    top,
    left,
    animation,
    particle
  };


  //////////////////////////////////////Définitions des structures décrivant les collisions dans Event///////////////////////////////////////////////////
  typedef struct {
    enum col_type type; /*Collision avec quel type de mur*/
    int ia; /*indice particule*/
    int ib;/*indice de l'autre particule lorsqu'il y a un choc particule-particule*/
    double time; /*Temps de colision*/
  } Event;



  ////////////////////////////////////////////Renvoit un signe aléatoire///////////////////////////////////////////////////////////////////////////////////
  int SigneAleatoire(void){
    double a=drand48();
    if (a>0.5){
      return(-1);
    }
    else{
      return(1);
    }
  }


  ////////////////////////////////////////Positionnement et teste de superposition///////////////////////////////////////////////////////////////

  void InitialisationPositionParticuleMobile(Particle * p, double l1x, double l2x, double l1y, double l2y, int i){
    int test;
    float dx,dy,Rayon;
    test=1;/*Pour entrer dans la boucle*/

    while (test==1) {/*Test pour empêcher la superposition des boules*/
      test = 0;
      p[i].x = l1x + p[i].Rad + (l2x-l1x-2*p[i].Rad)*drand48(); /*Random initial positions*/
      p[i].y = l1y + p[i].Rad + (l2y-l1y-2*p[i].Rad)*drand48();

      for(int j=0;j<i;j++){/*Pour toutes les boules créées antérieurement, on test si la boule i ne s'y superpose pas*/
        dx=abs(p[j].x-p[i].x);
        dy=abs(p[j].y-p[i].y);
        Rayon=abs(p[i].Rad+p[j].Rad);

        if (sqrt(pow(dx,2)+pow(dy,2))<Rayon){
          test=1;/*Si c'est le cas on recalcule une position*/
        }

      }
    }
  }



/////////////////////////////////////////Initialisation de la vitesse des particules/////////////////////////////////////////////////////
  void InitialisationVitesseParticuleMobile(Particle *p, int **InfoCases, int i, int v){
    if (InfoCases[i][2]==0){ /*Aucun confinement, vitesse élevée*/
      p[v].vx = 0.5*SigneAleatoire();
      p[v].vy = 0.5*SigneAleatoire();
    }

    else if (InfoCases[i][2]==1){/*Confinement peu respecté*/
      p[v].vx = 0.2*SigneAleatoire();
      p[v].vy = 0.2*SigneAleatoire();
    }

    else if (InfoCases[i][2]==2){/*Confinement stricte, presque à l'arrêt*/
      p[v].vx = 0.1*SigneAleatoire();
      p[v].vy = 0.1*SigneAleatoire();
    }
  }




////////////////////////////////////////////Initialisation paramêtres d'une particule fixe//////////////////////////////////////////////////////////

  void InitialisationParticuleFixe(Particle *particule, int IndiceParticule, double RayonParticuleFixe){

    particule[IndiceParticule].Rad=RayonParticuleFixe;
    particule[IndiceParticule].Mass=100000;
    particule[IndiceParticule].vx=0;
    particule[IndiceParticule].vy=0;
    particule[IndiceParticule].State=3;
    particule[IndiceParticule].Zone=-1;
    particule[IndiceParticule].Contamination_time = -1;
  }





//////////////////////////////////////////Initialisation de la position des particules fixes/////////////////////////////////////////////////
  int InitialisationDesParticulesFixes(Particle *p, double Lmin_x, double Lmin_y, double LCases, double DistInterPart, double RayonPartFixe, int NombreCases, int NombrePartFixesParCote){

    int IndicePartFixe = 0;

    for(int i = 1; i < NombreCases; i++){  /*On place d'abord les particules fixes en lignes horizontales*/
      for(int j = 0; j < NombrePartFixesParCote; j++){/*Plusieurs lignes en même temps*/
        p[IndicePartFixe].x = Lmin_x + j*DistInterPart;
        p[IndicePartFixe].y = LCases*i;
        InitialisationParticuleFixe(p, IndicePartFixe, RayonPartFixe);
        IndicePartFixe ++;
      }
    }

    for(int i = 1; i < NombreCases; i++){/*On crée ensuite les lignes verticales*/
      for(int j = 0; j < NombrePartFixesParCote; j++){/*Plusieurs lignes en même temps*/
        if (j%((NombrePartFixesParCote-1)/NombreCases) != 0 || j == 0 || j == NombrePartFixesParCote-1){
          p[IndicePartFixe].x = LCases*i;
          p[IndicePartFixe].y = Lmin_y + j*DistInterPart;
        //On regarde si l'itération est proportionnelle au nombre de trous par arete,si c'est le cas la particule a déjà été initialisée dans la boucle d'avant, donc on ne le refait pas
          InitialisationParticuleFixe(p, IndicePartFixe, RayonPartFixe);
          IndicePartFixe ++;
        }
      }
    }
    return IndicePartFixe;
  }






  //////////////////////////////////////////////////Initialisation des particules//////////////////////////////////////////////////////////////////////
  void InitialisationParticules(Particle *p, double Radius, double ParticleMass, double Lmin_x, double Lmin_y, double LCases, double DistInterPart, double RayonPartFixe, int NombreCases, int NombrePartFixesParCote, int** InfoCases, double** DimCases, double delai_guerison) {

    double random;
    int etat_aleatoire=0;
    int v=0;
    int IndicePartFixe = InitialisationDesParticulesFixes(p, Lmin_x, Lmin_y, LCases, DistInterPart, RayonPartFixe, NombreCases, NombrePartFixesParCote);

  ///////////////////////////////////////////Initialisation des particules mobiles
    v = IndicePartFixe;
    for (int i=0; i<NombreCases*NombreCases;i++){//Pour chaque case
      for (int j=0; j<InfoCases[i][0]; j++){ //Pour chaque particule de la case i

        random = drand48();/*Etat de la bille, définie aléatoirement*/

        if(random < 0.02*pow(2, InfoCases[i][1])){ // si InfoCases[i][1] == 0 -> taux d'infectés = 0.02; resp. 1 -> 0.04; 2 -> 0.08
          etat_aleatoire = 1; // etat infecté
        }else{
          etat_aleatoire = 0; // etat sain
        }
        p[v].State = etat_aleatoire;
        p[v].Contamination_time = drand48()*delai_guerison;
        p[v].Zone = i+1;/*Zone où est initialisée la particule*/
        p[v].Mass = ParticleMass;
        p[v].Rad = Radius;

        InitialisationPositionParticuleMobile(p, DimCases[i][0], DimCases[i][1], DimCases[i][2], DimCases[i][3], v);
        InitialisationVitesseParticuleMobile(p, InfoCases, i, v);

        v++;
      }
    }
  }



//////////////////////////Correction des imprecisions sur les positions////////////////////////////////////////////////////////////////////////
float Recadrage(float a, int b, float Rad, float Lmax, float Lmin){

  if (b==0){/*Recadrage sur x*/
    if (a<Lmin+Rad-10e-15){
      a=Lmin+Rad+ 10e-15;
    }
    else if (a> Lmax-Rad+ 10e-15){
      a=Lmax-Rad-10e-15;
    }
  }
  else {/*Recadrage sur y*/
    if (a<Lmin+Rad-10e-15){
      a=Lmin+Rad+ 10e-15;
    }
    else if (a>Lmax-Rad+ 10e-15){
      a=Lmax-Rad-10e-15;
    }
  }
return(a);
}



///////////////////////////////////////////////Remplissage du tableau des dimensions des zones///////////////////////////////////////////////////
void RemplissageTableauDimensionsCases(double ** DimCases, int NombreCases, double LCases, double Lmax_x, double Lmin_x, double Lmax_y, double Lmin_y){
  for(int i=0; i<NombreCases; i++){/*Ligne*/
    for(int j=0; j<NombreCases; j++){/*Colonne*/
      printf("Dimensions de la zone %d", j+ NombreCases*i);

      DimCases[j+(NombreCases*i)][0] = j*LCases+Lmin_x;
      printf("lmin_x= %g,  ",   DimCases[j+NombreCases*i][0]);

      if (j==NombreCases-1){
        DimCases[j+(NombreCases*i)][1] =Lmax_x;
      }else {
      DimCases[j+NombreCases*i][1] = (j+1)*LCases+Lmin_x;
      }
      printf("lmax_x= %g,  ",   DimCases[j+NombreCases*i][1]);

      DimCases[j+NombreCases*i][2] = i*LCases+Lmin_y;
      printf("lmin_y= %g,  ",   DimCases[j+NombreCases*i][2]);

      if (i==NombreCases-1){
        DimCases[j+NombreCases*i][3] =Lmax_y;
      }else {
      DimCases[j+NombreCases*i][3] = (i+1)*LCases+Lmin_y;
      }
      printf("Lmax_y= %g",   DimCases[j+NombreCases*i][3]);
      printf("\n");

    }
    printf("\n");
  }
}


/////////////////////////////////////Remplissage du tableau des informations de cases/////////////////////////////////////////////////////////
int RemplissageTableauInformationsCases(int ** InfoCases, int NombreCases, int NombrePartFixesParCote){

  int Np=0;
  int a=0, b=0, c=0;
  printf("Vous allez devoir donner le nombre d'habitants par region \n  le taux d'infection : entre 0 et 2 \n et la qualité du confinement : 0 = Nulle ou 1 = intermédiaire ou 2 = confinement total \n");
  for (int j=0; j<NombreCases*NombreCases;j++){

    printf("Nombre d'habitants dans la zone %d : \n", j+1);
    c=scanf("%d", &InfoCases[j][0]);
    if (c==0){
      printf("ERREUR, IL FAUT UN NOMBRE\n");
    }

    printf("Taux d'infection ?\n");
    a=scanf("%d", &InfoCases[j][1]);
    if (InfoCases[j][1]>2 || InfoCases[j][1]<0){
      printf("Veuillez repréciser le taux d'infection par un chiffre entre 0 et 2\n");
      a=scanf("%d", &InfoCases[j][1]);
    }
    else if(a==0){
      printf("ERREUR, IL FAUT UN NOMBRE\n");
    }


    printf("Qualité du confinement \n");
    b = scanf("%d", &InfoCases[j][2]);
    if (InfoCases[j][2]>2 ||InfoCases[j][2]<0){
      printf("Veuillez repréciser la qualité du confinement par un chiffre entre 0 et 2 \n");
      b=scanf("%d", &InfoCases[j][2]);
    }
    else if (b==0){
      printf("ERREUR, IL FAUT UN NOMBRE\n");
    }

    Np+=InfoCases[j][0];/*On somme pour avoir le nombre total de particules voulues*/
  }

  Np += 2*(NombreCases-1)*NombrePartFixesParCote;//En comptant les particules des murs
  printf("Nombre total de particules = %d \n", Np);


  return(Np);
}


///////////////////////////////////////////Créations des évènements individuels//////////////////////////////////////////////////////////////////////
void CreationEvenementChocMurDroite(Event *e, Particle *p, int i, int compteur, double Lmax_x){
  e[compteur].ia = i;
  e[compteur].ib=-1;
  e[compteur].type = right;

  if (p[i].vx > 0) {
    e[compteur].time=(Lmax_x-p[i].Rad-p[i].x)/p[i].vx;
  }
  /*calcul du temps avant la prochaine collision pour la particule i (on prend en compte 4 types de coll, d'où compt qui permet d'avancer)*/
  else {
    e[compteur].time = -1;/*Si la particule va vers la droite, on ignore la collision sur le mur de gauche*/
  }

}


void CreationEvenementChocMurGauche(Event *e, Particle *p, int i, int compteur, double Lmin_x){
  e[compteur].ia = i;
  e[compteur].ib=-1;
  e[compteur].type = left;

  if (p[i].vx > 0){
    e[compteur].time=-1;
  }
  else{
    e[compteur].time = (Lmin_x + p[i].Rad - p[i].x)/p[i].vx;
  }
}


void CreationEvenementChocMurHaut(Event *e, Particle *p, int i, int compteur, double Lmax_y){
  e[compteur].ia = i;
  e[compteur].ib=-1;
  e[compteur].type = top;

  if (p[i].vy > 0) {
    e[compteur].time = (Lmax_y-p[i].Rad-p[i].y)/p[i].vy;
  }
  else{
    e[compteur].time = -1;
  }
}


void CreationEvenementChocMurBas(Event *e, Particle *p, int i, int compteur, double Lmin_y){
  e[compteur].ia = i;
  e[compteur].ib=-1;
  e[compteur].type = bottom;

  if (p[i].vy > 0){
    e[compteur].time = -1;
  }
  else{
    e[compteur].time = (Lmin_y+p[i].Rad-p[i].y)/p[i].vy;
  }
}


void CreationEvenementsChocsParticules(Event *e, Particle *p, int i, int j, int compteur){

  float dx, dy, dvx, dvy, A, B, C, Omega;/*Pour les calculs de temps de collisions*/

  dx = p[i].x-p[j].x; /*Cf la fichier de Julian sur la resolution théorique*/
  dy= p[i].y-p[j].y;
  dvx= p[i].vx-p[j].vx;
  dvy= p[i].vy-p[j].vy;
  A=pow(dvy,2)+pow(dvx,2);/*pow permet de mettre à la puissance 2*/
  B=2*((dvx*dx)+(dvy*dy));
  C=pow(dx,2)+pow(dy,2)-pow(p[i].Rad+ p[j].Rad,2);
  Omega=pow(B,2)-(4*A*C);

  if (pow(B,2)<4*A*C|| B>0 ){/*Pour éviter de prendre en compte des évènements impossibles*/
    e[compteur].time=-1;
  }
  else {
    e[compteur].time=(-B-sqrt(Omega))/(2*A);
  }
  e[compteur].ia=i;
  e[compteur].ib=j;
  e[compteur].type = particle;
}


///////////////////////////////////////////Création du tableau des évènements/////////////////////////////////////////////////////////////////////
void CreationDesEvenements(Event *e, Particle* p, float Vraidt, int Np, double Lmax_x, double Lmax_y, double Lmin_x, double Lmin_y){


  int compteur = 0;/*va de 1 à 4*Np+1+Np*(Np-1)/2*/
  /////////////////Création de l'évenement fictif 'Animation' pour fluidifier l'affichage
  e[compteur].ia=-1;//Evenement 0 = Animation
  e[compteur].ib=-1;/*Pour les collisions non-interparticulaires, on met -1 pour ib*/
  e[compteur].time=Vraidt; /*Et non dt, afin d'avoir des images à temps dt régulier (même quand il y a des phenomènes physiques entre)*/
  e[compteur].type=animation;
  compteur ++;

  for (int i=0;i<Np;i++){
    CreationEvenementChocMurDroite(e, p, i, compteur, Lmax_x);
    compteur ++;
    CreationEvenementChocMurGauche(e, p, i, compteur, Lmin_x);
    compteur ++;
    CreationEvenementChocMurHaut(e, p, i, compteur, Lmax_y);
    compteur ++;
    CreationEvenementChocMurBas(e, p, i, compteur, Lmin_y);
    compteur ++;
    for(int j=i+1; j<Np; j++){
      CreationEvenementsChocsParticules(e, p, i, j, compteur);
      compteur ++;
    }
  }
}


///////////////////////////////////////////Mise à jours des paramètres associés aux sphères/////////////////////////////////////////////////
void MiseAJourPositionsParticules(Particle *p, int Np, float tmin){

  for (int i=0;i<Np;i++){
    p[i].x=p[i].x+(tmin*p[i].vx);
    p[i].y=p[i].y+(tmin*p[i].vy);
  }
}


void TestDeRemission(Particle*p, int Np, double dt, double delai_guerison){
  for (int i=0;i<Np;i++){
    if(p[i].State == 1){ /*Test de rémission*/
      p[i].Contamination_time += dt;
      if (p[i].Contamination_time > delai_guerison){
        p[i].State = 2;
      }
    }
  }
}



void MiseAJourEtatsParticules(Event *e, Particle *p, float nombre_aleatoire, int k){

  if(p[e[k].ia].State == 1 && p[e[k].ib].State == 0){
    if (nombre_aleatoire < 0.3){/*La probabilité qu'un contact soit infectueux est estimée à 0.3*/
      p[e[k].ib].State = 1;
      p[e[k].ib].Contamination_time = 0;
    }
  }else if(p[e[k].ia].State == 0 && p[e[k].ib].State == 1){
    if (nombre_aleatoire < 0.3){
      p[e[k].ia].State = 1;
      p[e[k].ia].Contamination_time = 0;
    }
  }
}

void AffichageEvenementAnimation(int a,double * Compteurtemps, double ** DimCases, Particle * p, double * Vraidt, Graphics * gw, double dt, double FPS, int NombreCases, double **TableauDeCouleurs, int **CompteurZones, int** InfoCases, FILE * Epidemie){

  gw->draw(p, FPS, NombreCases, TableauDeCouleurs, CompteurZones, DimCases, InfoCases, Compteurtemps);
  *Vraidt=dt;
  int NombreTotalDeMalades=0;
  *Compteurtemps+=dt;
  fprintf(Epidemie,"%g     ", *Compteurtemps);
  for (int m=0; m<(NombreCases*NombreCases); m++){
    NombreTotalDeMalades+=CompteurZones[m][2];
    if (a==0){
      fprintf(Epidemie, "     %d     ", (CompteurZones[m][2]*100)/CompteurZones[m][0]); /*On enregistre dans le fichier les valeurs du compteur en pourcentage,
    calculées dans draw*/
    }
  }
  if (a==1){
      fprintf(Epidemie, "     %d     ", NombreTotalDeMalades );
  }
  fprintf(Epidemie, "\n");
}



void MiseAJourVitessesParticules(int Np, double delai_guerison,int a,double * Compteurtemps, double ** DimCases, Event * e, Particle * p, double * Vraidt, Graphics * gw, double dt, double FPS, int NombreCases, double **TableauDeCouleurs, int **CompteurZones, int** InfoCases, FILE * Epidemie, double tmin, int k, double Lmax_x, double Lmax_y, double Lmin_x, double Lmin_y){

  float Dx, Dy, DxUni, DyUni, DVx, DVy, ProdScal, Norme;
  float nombre_aleatoire = drand48();

  if(e[k].type==right || e[k].type==left){
    p[e[k].ia].x=Recadrage(p[e[k].ia].x,0,p[e[k].ia].Rad,Lmax_x, Lmin_x);
    p[e[k].ia].vx *= -1;/*inversion du signe de la composante d'interet de la vitesse*/
    *Vraidt=*Vraidt-tmin;
  }

  else if (e[k].type==top || e[k].type==bottom){
    p[e[k].ia].y=Recadrage(p[e[k].ia].y,1,p[e[k].ia].Rad,Lmax_y, Lmin_y);
    p[e[k].ia].vy *= -1;
    *Vraidt=*Vraidt-tmin;/*On enlève le tmin afin de rester fluide dans l'affichage*/
  }

  else if(e[k].type==particle){/*Si c'est un choc entre particules*/

    if (p[e[k].ia].State != 3){
      p[e[k].ia].x=Recadrage(p[e[k].ia].x,0,p[e[k].ia].Rad,Lmax_x, Lmin_x);
      p[e[k].ia].y=Recadrage(p[e[k].ia].y,1,p[e[k].ib].Rad,Lmax_y, Lmin_y);
    }
    if (p[e[k].ib].State != 3){
      p[e[k].ib].x=Recadrage(p[e[k].ib].x,0,p[e[k].ib].Rad,Lmax_x, Lmin_x);
      p[e[k].ib].y=Recadrage(p[e[k].ib].y,1,p[e[k].ib].Rad,Lmax_x, Lmin_y);
    }

    Dx=p[e[k].ib].x-p[e[k].ia].x;
    Dy=p[e[k].ib].y-p[e[k].ia].y;
    Norme=sqrt(pow(Dx,2)+pow(Dy,2));
    DxUni=Dx/Norme;
    DyUni=Dy/Norme;
    DVx=p[e[k].ia].vx-p[e[k].ib].vx;
    DVy=p[e[k].ia].vy-p[e[k].ib].vy;
    ProdScal=((DxUni*DVx)+(DyUni*DVy));
    p[e[k].ia].vx=p[e[k].ia].vx-((2*p[e[k].ib].Mass)/(p[e[k].ia].Mass+p[e[k].ib].Mass))*ProdScal*DxUni; /*Cf fichier de Julian sur la mise à jour des vitesses !*/
    p[e[k].ia].vy=p[e[k].ia].vy-((2*p[e[k].ib].Mass)/(p[e[k].ia].Mass+p[e[k].ib].Mass))*ProdScal*DyUni;
    p[e[k].ib].vx=p[e[k].ib].vx+((2*p[e[k].ia].Mass)/(p[e[k].ia].Mass+p[e[k].ib].Mass))*ProdScal*DxUni;
    p[e[k].ib].vy=p[e[k].ib].vy+((2*p[e[k].ia].Mass)/(p[e[k].ia ].Mass+p[e[k].ib].Mass))*ProdScal*DyUni;
    *Vraidt=*Vraidt-tmin;

    MiseAJourEtatsParticules(e, p, nombre_aleatoire, k);


  }else{
      TestDeRemission(p,Np, dt,delai_guerison);
      AffichageEvenementAnimation(a,Compteurtemps, DimCases, p, Vraidt, gw, dt, FPS, NombreCases, TableauDeCouleurs, CompteurZones, InfoCases, Epidemie);

  }
}


//////////////////////////////////////Calcul du temps minimal////////////////////////////////////////////////////////////////////////////////
float IndiceProchainEvenement(Event * e, int N){
  double tmin=e[0].time;/*On commence avec dt, si on trouve une collision à un temps inférieur OU EGALE à dt,
  c'est la collision qui prime (sinon la boule sort de la boite)*/
  int k=0;
  for(int j=0;j<N;j++){
    if (e[j].time > 0 && e[j].time<=tmin){
      tmin = e[j].time;
      k=j;
    }
  }
  return(k);
}


////////////////////////////////////////Calcul du prochain évènement/////////////////////////////////////////////////////////////////////
void RechercheProchainEvenement(int a,double * Compteurtemps, Event * e, double ** DimCases, Particle * p, double * Vraidt, Graphics * gw, double dt, double FPS, int Np, double Lmax_x, double Lmax_y, double Lmin_x, double Lmin_y, double delai_guerison, int NombreCases, int TailleEv, double **TableauDeCouleurs, int **CompteurZones, int** InfoCases, FILE * Epidemie){

  int k = IndiceProchainEvenement(e, TailleEv);
  float tmin = e[k].time;
  MiseAJourPositionsParticules(p, Np, tmin);
  MiseAJourVitessesParticules(Np, delai_guerison,a,Compteurtemps, DimCases, e, p, Vraidt, gw, dt, FPS, NombreCases, TableauDeCouleurs, CompteurZones, InfoCases, Epidemie, tmin, k, Lmax_x, Lmax_y, Lmin_x, Lmin_y);

}


////////////////////////////////////////Generation des couleurs des particules en fonction des zones/////////////////////////////////////////////////
void GenerationDeCouleursAleatoires(double **tableauCouleurs, int NombreCases){
  int i,j;
  for(i=0; i<NombreCases*NombreCases; i++){
    for(j=0; j<3; j++){
      tableauCouleurs[i][j] = drand48();
    }
  }
}


int main(){
////////////////////////////////////Initialisation////////////////////////////////////////////////////////////////////
  double FPS=60., dt=0.1; /*dt = 1 heure*/
  double temps_guerison = 14*24*dt; /*temps de guerison estimé : 20 jours*/
  double Vraidt=dt;
  double diameter=0.2, Radius=diameter/2, ParticleMass=0.2;
  int Pix=750; /*Nombre de pixels pour la fenêtre*/
  double Lmax_x=10, Lmax_y = 10, Lmin_x=0, Lmin_y=0;



  int NombreCases=3; /*Nombre de cases par côté*/
  int NombreAretes = NombreCases*(2*NombreCases-2); /*Nombre d'arêtes (côtés des cases) hors parois limites*/
  int NombreTrous = 4*NombreAretes; /*Il faut qu'il y ait le même nombre de trous sur chaque arête,
  donc que ce nombre soit proportionnel au  nombre d'arêtes*/
  double LargeurTrous=3*Radius; /*largeur des trous de la parois poreuse*/
  double LongueurCarre = Lmax_x - Lmin_x;
  double LCases = LongueurCarre/NombreCases; /*Longueur d'une arête d'une case*/
  double RayonPartFixe = (((NombreAretes*LongueurCarre)/(NombreTrous*NombreCases)) - LargeurTrous)/2;
  /*Le nombre de trous est égal au nombre de particules fixes*/
  int NombreParticulesFixesParCote = (NombreTrous/NombreAretes)*NombreCases + 1;
  double DistInterPart = 2*RayonPartFixe + LargeurTrous;



  int **InfoCases;/*Tableau du nombre d'habitants par case, du taux d'infection, et de la qualité du conf*/
  InfoCases= (int **)malloc((NombreCases*NombreCases)* sizeof(int *));
  for (int i=0; i<(NombreCases*NombreCases); i++){
    InfoCases[i]=(int *)malloc(3 *sizeof(int));
  }


  double **DimCases;/*Tableau contenant les bornes de chacune des cases (lmin_x, lmax_x, lmin_y, lmax_y)*/
  DimCases= (double **)malloc((NombreCases*NombreCases)* sizeof(double *));
  for (int i=0; i<(NombreCases*NombreCases); i++){
    DimCases[i]=(double *)malloc(4 *sizeof(double));
  }


  double **TableauDeCouleurs; /*Tableau contenant les codes rgb des couleurs des particules dans chaque case*/
  TableauDeCouleurs = (double **)malloc((NombreCases*NombreCases)* sizeof(double *));
  for (int i=0; i<(NombreCases*NombreCases); i++){
    TableauDeCouleurs[i]=(double *)malloc(3 *sizeof(double));
  }
  GenerationDeCouleursAleatoires(TableauDeCouleurs, NombreCases);


  RemplissageTableauDimensionsCases(DimCases, NombreCases, LCases, Lmax_x, Lmin_x, Lmax_y, Lmin_y);
  int Np = RemplissageTableauInformationsCases(InfoCases, NombreCases, NombreParticulesFixesParCote);
  int TailleEv = 4*Np+1+ Np*(Np-1)/2; /*Taille du tableau des evenements*/


  Graphics gw(Np, Pix, Lmin_x, Lmin_y, Lmax_x, Lmax_y, diameter);/*Open a window to plot particles in*/
  srand48(time(NULL));/*inititalize random numbers*/
  Particle *p= (Particle *) malloc(Np *sizeof( Particle)); /*an array of particles*/
  InitialisationParticules(p, Radius, ParticleMass, Lmin_x, Lmin_y, LCases, DistInterPart, RayonPartFixe, NombreCases, NombreParticulesFixesParCote, InfoCases, DimCases, temps_guerison); /*place particles in box*/
  Event *e = (Event *) malloc( (TailleEv)* sizeof(Event)); /*Tableau des evenements*/


  int ** CompteurZones; /*Tableau qui tient le compte des individus malades, remis, sains et total dans chaque case*/
  CompteurZones = (int **)malloc((NombreCases*NombreCases)* sizeof(int *));
  for (int i=0; i<(NombreCases*NombreCases); i++){
    CompteurZones[i]=(int *)malloc(4 *sizeof(int));
  }


  double Compteurtemps=0;/*Décompte du temps pour le tracé*/
  FILE* Epidemie = NULL;/*Fichier pour tracer le nombre de malades dans chaque case*/
  Epidemie = fopen("Epidemie.txt", "w");/*w pour write*/
  int a=0, rep;
  printf("Si vous voulez tracer le nombre total de malades en fonction du temps, tapez 1, si vous voulez tracer le pourcentage de malades par zones, tapez 0\n");
  rep= scanf("%d", &a);
  if (rep==0 || a>1){
    printf("Veuillez mieux répondre à la question svp \n");
    rep=scanf("%d", &a);
  }
  if (Epidemie != NULL){



////////////////////////////////////  ANIMATION   //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    gw.draw(p,FPS, NombreCases, TableauDeCouleurs, CompteurZones, DimCases, InfoCases, &Compteurtemps); /*Position initiale*/
    int k, l;
    for (l=0; l<100000;l++){
      CreationDesEvenements(e, p, Vraidt, Np, Lmax_x, Lmax_y, Lmin_x, Lmin_y);
      RechercheProchainEvenement(a, &Compteurtemps, e, DimCases, p, &Vraidt, &gw, dt, FPS, Np, Lmax_x, Lmax_y, Lmin_x, Lmin_y, temps_guerison, NombreCases, TailleEv, TableauDeCouleurs, CompteurZones, InfoCases, Epidemie);
    }
    fclose(Epidemie);
  }

  free(InfoCases);
  free(TableauDeCouleurs);
  free(CompteurZones);
  free(DimCases);
  free(p);
  free(e);


  return(0);
}
