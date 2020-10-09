/*
 * Resolution_complete.cpp
 *
 *  Created on: 29 d»c. 2019
 *      Author: imane
 */
#include "Put.h"
#include "Call.h"
#include <cmath>
#include <stdlib.h>
#include "Resolution_complete.h"
#include "general.h"
#include "black-scholes.h"

namespace edp{
/** Constructeur valué de la classe*/
Resolution_complete::Resolution_complete(double T, double L,double r,double sigma,double K, int N, int M){
	dT=T/M;
	dS=L/(N+1);
	/* Allocation de m√àmoire pour la matrice A*/
	A=new double[N*N];
	/*On remplit la matrice A de 0*/
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			A[N*i+j]=0.0;
		}
	}/**remplissage de la diagonale  */
	for(int i=0;i<N;i++){
		A[N*i+i]= 1+0.5*r*dT*(1+i)+0.5*pow(sigma,2)*pow(i+1,2)*dT;
	}
	/**remplissage de la diagonale inferieur */
	for(int i=1;i<N;i++){
		A[N*i+(i-1)]= -0.25*pow(sigma,2)*dT*pow(i+1,2);
	}
	/**remplissage de la diagonale sup√©rieur */
	for(int i=0;i<N-1;i++){
		A[N*i+(i+1)]= -0.5*r*dT*(i+1) - 0.25*pow(sigma,2)*dT*pow(i+1,2);
	}
	/* Allocation de m»moire pour la matrice B*/
	B=new double[N*N];
	/*On remplit la matrice B de 0*/
	for(int i=0;i<=N-1;i++){
		for(int j=0;j<=N-1;j++){
			B[(i*N)+j]=0.0;
		}
	}
	/**remplissage de la diagonale  */
	for(int i=0;i<N;i++){
		B[(i*N)+i]= 1-0.5*r*dT*(1+i)-0.5*pow(sigma,2)*dT*pow(i+1,2)-r*dT;
	}
	/**remplissage de la diagonale infÈrieur */
	for(int i=1;i<N;i++){
		B[i*N+(i-1)]= 0.25*pow(sigma,2)*dT*pow(i+1,2);
	}
	/**remplissage de la diagonale supÈrieur */
	for(int i=0;i<N-1;i++){
		B[N*i+(i+1)]= 0.5*r*dT*(i+1) + 0.25*pow(sigma,2)*dT*pow(i+1,2);
	}

}
/** Destrueteur de la classe */
Resolution_complete::~Resolution_complete(){
	delete[]A;/** on supprime A*/
	delete[]B;/** on supprime B*/
}
/**
 * methode permettant de calculer le veceteur Fn de la résolution
 *  dans le cas ou le pay_off est un PUt*/
double* Resolution_complete::VecteurF_put(double K,double sigma,double r,int n,int N){
	/*Allocation de la m√©moire pour le vecteur F*/
	double * F=new double[N];
	/*On remplit le vecteur  F*/
	F[0]=K*exp(r*n*dT);
	for(int i=1;i<N;i++){
		F[i]=0.0;
	}

	return F;

}
/**
 * methode permettant de calculer le veceteur Fn de la résolution
 *  dans le cas ou le pay_off est un call*/
double* Resolution_complete::VecteurF_call(double K,double sigma,double r,int n,int N){
	/*Allocation de la mémoire pour le vecteur F*/
	double * F=new double[N];
	/*On remplit la matrice F*/
	F[N-1]=(K*exp(r*n*dT));
	for(int i=0;i<(N-1);i++){
		F[i]=0.0;
	}
	return F;
}
/**methode permettant de recupere le membre de donnes A*/
double* Resolution_complete::get_A(){
	return A;
}
/**methode permettant de recupere le membre de donnes B*/
double* Resolution_complete::get_B(){
	return B;
}
/**methode permettant de recupere le membre de donnes dT*/
double  Resolution_complete::get_dt(){
	return dT;
}
/**methode permettant de recupere le membre de donnes dS*/
double  Resolution_complete::get_ds(){
	return dS;
}
/** methode permettant de resoudre une edp de blackschole sachant que le payoff est un put
 *  ainsi que la discretisation en temps et en espace */
double * Resolution_complete::res_edp_put(BlackScholes Bs,int Mtemps,int Nespace ,double k){
	/**Recuperation de parametre de l'edp Bs*/
	double r=Bs.get_r();
	double L=Bs.get_L();
	double T=Bs.get_T();
	double sigma=Bs.get_sigma();
	/** création de l'objet calcul C */
	Resolution_complete C( T,  L, r, sigma,k, Nespace,  Mtemps);
	/** recuperation de la matrice B et de dS de calcul*/
	double *B=C.get_B();
	double ds=C.get_ds();
	/**initialisation du vecteur solution*/
	double *sol= new double[Nespace];
	for ( int i=0;i<Nespace;i++){
		double s_i=(i+1)*ds;
		sol[i]=max(0.0,k-s_i);//premiere colonne ligne i
	}
	for (int m=0;m<Mtemps;m++){
		/**Resolution_complete du vecteur Fn */
		double*Cn=Resolution_complete::VecteurF_put(k,sigma,r,m+1, Nespace);
		/**Resolution_complete du vecteur B*Un */
		double * Bn=Produitvecteur(B,sol,Nespace);
		double *b=SommeVecteur(Bn,Cn,Nespace);
		/** resolution du systeme AUn+1=BUn+ Fn */
		double* res=res_LU(C.get_A(),b,Nespace);
		/** recuperation du vecteur solution*/
		for(int i=0;i<Nespace;i++){
			sol[i]=res[i];
		}
		/** suppresion des tableaux*/
		delete[]Cn;
		delete[]Bn;
		delete[]b;
		delete[]res;
	}
	std::cout<<"Solution complete pour un put \n";
	std::cout<<"--------------------\n";
	return sol;
}
/** methode permettant de resoudre une edp de blackschole sachant que le payoff est un call
 *  ainsi que la discretisation en temps et en espace */
double *  Resolution_complete::res_edp_call(BlackScholes Bs,int Mtemps,int Nespace ,double k){
	/**Recuperation de parametre de l'edp Bs*/
	double r=Bs.get_r();
	double L=Bs.get_L();//borne espace
	double T=Bs.get_T();
	double sigma=Bs.get_sigma();
	/** création de l'objet calcul C */
	Resolution_complete C( T,  L, r, sigma,k, Nespace,  Mtemps);
	/** recuperation de la matrice B et de dS de calcul*/
	double *B=C.get_B();
	double ds=C.get_ds();
	/**initialisation du vecteur solution*/
	double *sol= new double[Nespace];
	for ( int i=0;i<Nespace;i++){
		double s_i=i*ds;
		sol[i]=max(0.0,s_i-k);//premiere colonne ligne i
	}
	for (int m=0;m<Mtemps;m++){
		/**Resolution_complete du vecteur Fn */
		double*Cn=Resolution_complete::VecteurF_call(k,sigma,r,m+1, Nespace);
		/**Resolution_complete du vecteur B*Un */
		double * Bn=Produitvecteur(B,sol,Nespace);
		double *b=SommeVecteur(Bn,Cn,Nespace);
		/** resolution du systeme AUn+1=BUn+ Fn */
		double* res=res_LU(C.get_A(),b,Nespace);
		/** recuperation du vecteur solution*/
		for(int i=0;i<Nespace;i++){
			sol[i]=res[i];
		}

		delete[]Cn;
		delete[]Bn;
		delete[]b;
		delete[]res;
	}
	std::cout<<"Solution complete pour un call\n";
	std::cout<<"--------------------\n";
	/** suppresion des tableaux*/
	return sol;
}
/**
 * methode permettant de resoudre une edp de blackschole connaissant le
 * strike du payoff ainsi que la discretisation en temps et en espace,
 * LA methode affiche le resultat  */
double * Resolution_complete::res_edp(BlackScholes Bs,int Mtemps,int Nespace ,double k){
	/**Recuperation de parametre de l'edp Bs*/
	char a=(*(Bs.get_payoff())).Pay_Off::get_type();
	double *sol;
	if(a=='p'){
		/** le payoff est un put
		 * utilisation de la fonction resolution dans le cas d'un put*/
		sol=Resolution_complete::res_edp_put(Bs, Mtemps, Nespace , k);
	}
	else{
		/** le payoff est un call
		 * utilisation de la fonction resolution dans le cas d'un call*/
		sol=Resolution_complete::res_edp_call(Bs, Mtemps, Nespace , k);
	}
	/** Affichage de la solution*/

	affichage_vecteur(sol,Nespace);
	return(sol);

}

}//namespace
