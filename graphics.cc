//https:www.cairographics.org/manual/cairo-Paths.html fix
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <time.h>
#include <assert.h>
#include "Graphics.h"


void
Graphics::draw(Particle *p, double FPS, int NombreDeCases, double **couleursZones, int** CompteurZones, double** DimCases, int** InfoCases, double *compteurTemps){
  assert(FPS > 1);
  //  const int FPS=50;//frames per second
  struct timespec tm={0,(long int)(1000000000/FPS)}; //Sleep in nanoseconds between frames
  XEvent event;//check if window closed and finish
  if( XCheckWindowEvent(dsp,da, DestroyNotify , &event)){ XCloseDisplay(dsp); exit(1);}

  double gamma = DimY/(lmax_y-lmin_y + diam);//scaling between physical units and pixels
  double alpha= gamma*diam/2.;//DO NOT Touch ME!


  char ***nbParticulesTxt;
  nbParticulesTxt = (char ***)malloc((NombreDeCases*NombreDeCases)* sizeof(char **));
  for (int i=0; i<(NombreDeCases*NombreDeCases); i++){
    nbParticulesTxt[i]=(char **)malloc(4 *sizeof(char*));
  }
  for(int i=0; i<NombreDeCases*NombreDeCases; i++){
    for(int j=0; j<4; j++){
      nbParticulesTxt[i][j]=(char *)malloc(20 *sizeof(char));
    }
  }


char *dtTxt; // contient le temps écoulé en heures
dtTxt = (char *)malloc(40*sizeof(char));


char **ZonesTxt;
ZonesTxt = (char **)malloc((NombreDeCases*NombreDeCases)*sizeof(char*));
for(int i=0; i<NombreDeCases*NombreDeCases; i++){
  ZonesTxt[i]=(char*)malloc(20*sizeof(char));
}


  cairo_push_group(cr); //start drawing
  cairo_set_source_rgb(cr, 0.0, 0.19, .19);//dark green background
  cairo_paint (cr); //clear screen with green

  if(1){// Concerne les particules solides



    for(int i=0; i<NombreDeCases*NombreDeCases; i++){// Remise à zéro des compteurs
      for(int j = 0; j < 4; j++){
        CompteurZones[i][j] = 0;
      }
    }


    for(int i=0;i<Np;i++){ // On parcourt l'intégralité des particules existantes
      for(int j = 0; j < NombreDeCases*NombreDeCases; j++){ //On parcourt chacune des zones
        if (p[i].x > DimCases[j][0] && p[i].x < DimCases[j][1] && p[i].y > DimCases[j][2] && p[i].y < DimCases[j][3] && p[i].Zone >=0 ){
          //La vitesse de la particule et sa zone sont modifiées si elle vient de changer de zones
          if (p[i].Zone!=j){
            p[i].Zone = j; // La zone de la particule est modifiée
            double Norme = sqrt(pow(p[i].vx,2)+pow(p[i].vy,2));
            if (InfoCases[j][2]==0){ /*Aucun confinement*/
              p[i].vx = 0.5*p[i].vx/Norme;/*Même direction mais norme différente*/
              p[i].vy = 0.5*p[i].vy/Norme;
            }
            else if (InfoCases[j][2]==1){/*Confinement peu respecté*/
              p[i].vx = 0.2*p[i].vx/Norme;
              p[i].vy = 0.2*p[i].vy/Norme;
            }
            else if (InfoCases[j][2]==2){/*Confinement stricte*/
              p[i].vx = 0.1*p[i].vx/Norme;
              p[i].vy = 0.1*p[i].vy/Norme;
            }
          }
          if (p[i].State == 0){
            cairo_set_source_rgb(cr, 0.49, 0.64, 0.0);//vert pour particule saine
            CompteurZones[j][1] += 1;
            CompteurZones[j][0] += 1; // à chaque particule parcourue on incrémente le nombre de particules totales présentes dans la zone j+1 considérée
          }
          if (p[i].State == 1){
            cairo_set_source_rgb(cr, 0.79, 0.29, 0.29);//rouge pour particule infectée
            CompteurZones[j][2] += 1;
            CompteurZones[j][0] += 1;
          }
          if (p[i].State == 2){
            cairo_set_source_rgb(cr, 1, 0.5, 0.72);//rose pour particule rétablie
            CompteurZones[j][3] += 1;
            CompteurZones[j][0] += 1;
          }

          cairo_new_sub_path(cr) ;
          cairo_arc(cr,   (alpha + gamma* (p[i].x -lmin_x)) ,  DimY-(alpha + gamma*(p[i].y - lmin_y)), gamma*p[i].Rad, 0, 2 * M_PI);// Après que la particule ait été caractérisée et comptée, on l'affiche (cercle)
          cairo_fill(cr);
          cairo_set_source_rgb(cr, couleursZones[j][0], couleursZones[j][1], couleursZones[j][2]);// on colorie le centre de la particule selon sa zone
          //printf("Zone %d", p[i].Zone);
          cairo_new_sub_path(cr) ;
          cairo_arc(cr,   (alpha + gamma* (p[i].x -lmin_x)) ,  DimY-(alpha + gamma*(p[i].y - lmin_y)), gamma*p[i].Rad - alpha/4, 0, 2 * M_PI);
          cairo_fill(cr);//draw all particles with solid color
        }
      }
      if(p[i].State == 3){
        cairo_set_source_rgb(cr, 0.34, 0.34, 0.34);//gris pour les particules fixes
        cairo_new_sub_path(cr) ;
        cairo_arc(cr,   (alpha + gamma* (p[i].x -lmin_x)) ,  DimY-(alpha + gamma*(p[i].y - lmin_y)), gamma*p[i].Rad, 0, 2 * M_PI);
        cairo_fill(cr);
      }
    }

//////////////////////////////////////////////////////////////////////////////////////////////  AFFICHAGE DU NOMBRE DE PARTICULES PAR ZONES   ////////////////////////////////////////////////////////////////////////////////////////////////////

    for (int k = 0; k <NombreDeCases*NombreDeCases; k++){
      for (int m = 0; m<4; m++){

        cairo_move_to(cr, gamma*lmax_x*(1.0*((k)%NombreDeCases)/NombreDeCases + 1.0/10), gamma*lmax_y*(1.0 - 1.0*((k)/NombreDeCases + 1)/NombreDeCases + 1.5/10 + 1.0/30*(m)));
        //On écrit dans le coin en haut à gauche de chaque zone le nombre de particules totales, saines infectées, rétablies respectives.

        switch (m){ // selon le type (total, saines, infectées , rétablies) on utilise une couleur différente
          case 0: // nombre total en blanc
          cairo_set_source_rgb(cr, 1, 1, 1);
          sprintf(nbParticulesTxt[k][m], "%d au total", CompteurZones[k][m]); // On remplit le tableau nbParticulesTxt qui contient ces valeurs sous forme de caractères
          break;
          case 1:// nombre de particules saines en vert
          cairo_set_source_rgb(cr, 0.49, 0.64, 0.0);
          sprintf(nbParticulesTxt[k][m], "%d sains", CompteurZones[k][m]); // On remplit le tableau nbParticulesTxt qui contient ces valeurs sous forme de caractères
          break;
          case 2: // nombre de particules infectées en rouge
          cairo_set_source_rgb(cr, 0.79, 0.29, 0.29);
          sprintf(nbParticulesTxt[k][m], "%d infectés", CompteurZones[k][m]); // On remplit le tableau nbParticulesTxt qui contient ces valeurs sous forme de caractères
          break;
          case 3: //nombre de particules retablies en rose
          cairo_set_source_rgb(cr, 1, 0.5, 0.72);
          sprintf(nbParticulesTxt[k][m], "%d rétablis", CompteurZones[k][m]); // On remplit le tableau nbParticulesTxt qui contient ces valeurs sous forme de caractères
          break;

        }


        cairo_set_font_size(cr, 15); // Police
        cairo_show_text (cr, nbParticulesTxt[k][m]);//affichage

      }

      cairo_set_source_rgb(cr, couleursZones[k][0], couleursZones[k][1], couleursZones[k][2]);
      cairo_move_to(cr, gamma*lmax_x*(1.0*((k)%NombreDeCases)/NombreDeCases + 1.0/10), gamma*lmax_y*(1.0 - 1.0*((k)/NombreDeCases + 1)/NombreDeCases + 1.5/20));
      sprintf(ZonesTxt[k], "Zone %d", k+1);
      cairo_set_font_size(cr, 20); // Police
      cairo_show_text (cr, ZonesTxt[k]);//affichage
    }



    if(0){//RING for outter circle
      cairo_set_source_rgb(cr, .2, 0.79, 0.79);//dark blue for particles
      for(int i=0;i<Np;i++){
        cairo_new_sub_path(cr) ;
        cairo_arc(cr,  (alpha + gamma* (p[i].x -lmin_x)) ,  DimY - (alpha + gamma*(p[i].y - lmin_y)), 1.4*gamma*p[i].Rad, 0, 2 * M_PI);
      }
      cairo_stroke (cr); // hollow particles
    }
    frame(alpha, (alpha + gamma* (lmax_x -lmin_x)), alpha, DimY-alpha);//draw square border
    cairo_pop_group_to_source(cr); //finished drawing operations for this set of positions
    cairo_paint(cr);//send to screen
    cairo_surface_flush(sfc); //send to x11
    nanosleep( &tm , NULL); //this sets the animations speed
    XFlush(dsp);//sync X11 to cairo


  //https://cairographics.org/Xlib/



  }


  cairo_set_source_rgb(cr, 1, 1, 1);//dark blue for particles
  cairo_move_to(cr, gamma*lmax_x*(0.55), gamma*lmax_y*(0.95));
  sprintf(dtTxt, "Nombre d'heures écoulées: %d", int(floor(10*(*compteurTemps))));
  cairo_set_font_size(cr, 20); // Police
  cairo_show_text (cr, dtTxt);//affichage

  free(nbParticulesTxt);
  free(ZonesTxt);




}

void Graphics::frame(double xmin, double xmax, double ymin, double ymax){
  cairo_set_source_rgb (cr, 1, 1, 0);//yellow border
  cairo_rectangle (cr,xmin, ymin, xmax-xmin, ymax-ymin);
  cairo_stroke (cr);
}
