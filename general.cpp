/*
 * general.cpp
 *
 *  Created on: 5 janv. 2020
 *      Author: mame diarra
 */
#include "general.h"
#include <iostream>
#include <cmath>
namespace edp{
/**Calcul de la matrice L de la decomposition LU d'un matrice A de taille n*n */
double* MatriceL(double* A,int n){
	int taille=n;
	double* S=new double[taille+1];
	S[0]=1.0;
	S[1]=A[taille*0+0];
	for(int k=2;k<=taille;k++){
		S[k]=A[taille*(k-1)+k-1]*S[k-1]-A[taille*(k-1)+(k-2)]*S[k-2]*A[taille*(k-2)+(k-1)];
	}
	double *L=new double[taille*taille];
	for(int i=0;i<taille;i++){
		for(int j=0;j<taille;j++){
			L[taille*i+j]=0.0;
		}
	}
	for(int i=0;i<taille;i++){
		L[taille*i+i]=1.0;
	}
	for(int i=1;i<taille;i++){
		L[i*taille+(i-1)]=A[i*taille+(i-1)]*(S[i-1]/S[i]);
	}
	delete []S;

	return L;
}
/**Calcul de la matrice U de la decomposition LU d'un matrice A de taille n*n */
double* MatriceU(double* A,int n){
	int taille=n;
	double* S=new double[taille+1];
	S[0]=1.0;
	S[1]=A[taille*0+0];
	for(int k=2;k<=taille;k++){
		S[k]=A[taille*(k-1)+k-1]*S[k-1]-A[taille*(k-1)+(k-2)]*S[k-2]*A[taille*(k-2)+(k-1)];
	}
	double *U=new double[taille*taille];
	for(int i=0;i<taille;i++){
		for(int j=0;j<taille;j++){
			U[(i*taille)+j]=0.0;
		}
	}
	for(int i=0;i<taille;i++){
		U[taille*i+i]=S[i+1]/S[i];
	}
	for(int i=0;i<(taille-1);i++){
		U[i*taille+(i+1)]=A[i*taille+(i+1)];
	}
	delete[]S;
	return U;

}
/** fonction utlisant la rdecomposition LU d'une matrice
 *  A e taille n pour resoudre le systeme Ax=b*/
double * res_LU(double *A, double *b,int n ){
	/** fait la decomppsition LU de A*/
	double *U= MatriceU(A,n);
	double *L= MatriceL(A,n);
	/** on rÈsoud Ly=b */
	double* y = new double [n];
	y[0]=b[0]/L[0*n+0];
	for (int i=1;i<n;i++){
		double s=0;
		for (int j=0;j<i;j++){
			s=s+L[i*n+j]*y[j];
		}
		y[i]=(b[i]-s)/L[i*n+i];
	}
	delete []L;
	/** on rÈsoud Ux=y */
	double * x= new double [n];
	x[n-1]=y[n-1]/U[(n-1)*n+n-1];
	for (int i=n-2;i>=0;i--){
		double s=0;
		for (int j=i+1;j<n;j++){
			s=s+U[i*n+j]*x[j];
		}
		x[i]=(y[i]-s)/U[i*n+i];
	}
	delete []U;
	delete[] y;
	return (x);
}
/** Fonction calculant le produit d'une matrice de taille n*n M1 et
 *  d'un vecteur colonne de longueur n V */
double* Produitvecteur( double * M1,  double* V,int n){
	double *res=new double[n];
	for(int i=0;i<n;i++){
		res[i]=0.0;
		for(int j=0;j<n;j++){
			res[i]=res[i]+M1[i*n+j]*V[j];
		}
	}
	return res;
}
/** Fonction calculant la somme de deux vecteurs colonnes de longueur n */
double* SommeVecteur(double * V1, double *V2,int n){
	double *sum=new double[n];
	for (int i=0;i<n;i++){
		sum[i]=V1[i]+V2[i];
	}
	return(sum);
}
/** Fonction affichant un vecteur colonne de taille n */
void affichage_matrice(double *D,int n,int p){
	for(int i=0;i<n;i++){
		for(int j=0;j<p;j++){
			std::cout<<D[i*p+j]<< "  ||  ";
		}
		std::cout <<"\n";
	}
}
/** Fonction affichant une matrice  de taille n*n */
void affichage_vecteur(double *D,int n){
	for(int i=0;i<n;i++){
		std::cout<<D[i]<<"\n";
	}
}
/**Fonction calculant le max entre deux rééls */
double max(double a,double  b){
	if(a>=b)
		return a;
	else
		return b;
}
/**
 * methode permettant de definir l'erreur relative entre la resolution
 * de l'edp de black & sholes completes et le retour inverse de l'edp
 * Resolution_reduite
 */
double *erreur_relative(double *a,double*b,int Nespace ){
	double *res=new double[Nespace];
	for(int i=0;i<Nespace;i++){
		res[i]=(abs(a[i]-b[i])/b[i]);
	}
	return res;
}
}//namespace
