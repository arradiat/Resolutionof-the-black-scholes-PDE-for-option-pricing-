/**
 * general.h
 *Classe contenant les methodes communes utilisant dans le cadre de la resolution
 *Classe de l'edp de black schole et de l'edp de black scholes réduite
 *
 *      Author: mame diarra
 */
#include <iostream>

#ifndef GENERAL_H_
#define GENERAL_H_
namespace edp{

/**Calcul de la matrice L de la decomposition LU d'un matrice A de taille n*n */
double* MatriceL(double* A,int n);

/**Calcul de la matrice U de la decomposition LU d'un matrice A de taille n*n */
double* MatriceU(double* A,int n);

/** fonction utlisant la rdecomposition LU d'une matrice
 *  A e taille n pour resoudre le systeme Ax=b*/
double * res_LU(double *A, double *b, int n);

/** Fonction calculant le produit d'une matrice de taille n*n M1 et
 *  d'un vecteur colonne de longueur n V */
double* Produitvecteur( double *M1, double* V,int n);
/** Fonction calculant la somme de deux vecteurs colonnes de longueur n */
double * SommeVecteur(double * V1,double*V2, int n);
/** Fonction affichant un vecteur colonne de taille n */
void affichage_vecteur(double *D,int n);
/** Fonction affichant une matrice  de taille n*n */
void affichage_matrice(double *D,int n,int p);

/**Fonction calculant le max entre deux rééls */
double max(double a , double b);
/**
	 * methode permettant de definir l'erreur relative entre la resolution
	 * de l'edp de black & sholes completes et le retour inverse de l'edp
	 * Resolution_reduite
	 */
	double *erreur_relative(double *a,double*b,int Nespace );


}

#endif /* GENERAL_H_ */
