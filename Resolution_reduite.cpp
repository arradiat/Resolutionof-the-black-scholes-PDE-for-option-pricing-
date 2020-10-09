/*
 * Resolution_reduite.cpp
 *
 *  Created on: 3 janv. 2020
 *      Author: mame diarra
 */
#include "Put.h"
#include "Call.h"
#include <cmath>
#include <stdlib.h>
#include "Resolution_reduite.h"
#include "Pay_off.h"
#include "general.h"
#include "black-scholes.h"

namespace edp{
/** Constructeur valuÈ de la classe*/
reduite::reduite(int M,int N, double T1, double  L1,int K,double sigma , double r){
	dS=(log(L1/K)-log(0.001/K))/((N+1));/*Discretisation spatiale: pas en espace de la resolution*/
	double Tprim=0.5*(sigma*sigma)*T1;

	dT=Tprim/M;/*Discretisation temporelle: pas en temps de la resolution*/
	T=T1;
	L=L1;
	A=new double[N*N];
	/*On remplit la matrice A*/
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++){
			A[N*i+j]=0.0;
		}
	}
	/**remplissage de la diagonale  */
	for(int i=0;i<N;i++){
		A[N*i+i]= 1+(2*dT)/(dS*dS);
	}
	/**remplissage de la diagonale inferieur */
	for(int i=1;i<N;i++){
		A[N*i+(i-1)]= -dT/(dS*dS);
	}
	/**remplissage de la diagonale supÈrieur */
	for(int i=0;i<N-1;i++){
		A[N*i+(i+1)]= -dT/(dS*dS);

	}
}
/** Destructeur de la classe */
reduite::~reduite(){

	delete[]A;/** on supprime A*/
}
/**
 * methode permettant de calculer le veceteur Fn de la r√©solution
 *  dans le cas ou le pay_off est un call*/
double* reduite::VecteurF_call(double K,double sigma,double r,int n,int N,int L){
	/*Allocation de la m√©moire pour le vecteur F*/
	double * F=new double[N];
	double k =2*r/(sigma*sigma);
	double w=exp(0.5*(k-1)*log(L/K)+0.25*(k+1)*(k+1)*n*dT)*exp((2*r*n*dT)/(sigma*sigma));
	/*On remplit le vecteur  F*/
	F[N-1]=(dT/(dS*dS))*w;
	for(int i=0;i<(N-1);i++){
		F[i]=0.0;
	}
	return F;
}
double* reduite::VecteurF_put(double K,double sigma,double r,int n,int N,int L){
	/*Allocation de la m√©moire pour le vecteur F*/
	double * F=new double[N];
	double k =2*r/(sigma*sigma);
	double w=(dT/(dS*dS))*(exp(0.5*(k-1)*log(0.001/K)+0.25*(k+1)*(k+1)*n*dT)*exp(-(2*r*n*dT)/(sigma*sigma)));
	/*On remplit la matrice F*/
	F[0]=w;
	for(int i=1;i<=(N-1);i++){
		F[i]=0.0;
	}
	return F;
}

/**methode permettant de recupere le membre de donnes A*/
double* reduite::get_A(){
	return A;
}
/**methode permettant de recupere le membre de donnes dT*/
double  reduite::get_dt(){
	return dT;
}
/**methode permettant de recupere le membre de donnes dS*/
double  reduite::get_ds(){
	return dS;

}

/** methode permettant de resoudre une edp de blackschole sachant que le payoff est un put
 *  ainsi que la discretisation en temps et en espace */

double * reduite::res_reduite_put(BlackScholes Bs, int Mtemps, int Nespace, double K){
	/**Recuperation de parametre de l'edp Bs*/
	double r=Bs.get_r();
	double L=Bs.get_L();//borne espace
	double T=Bs.get_T();
	double sigma=Bs.get_sigma();
	double k=2*r/(sigma*sigma);
	/** cr√©ation de l'objet Resolution_reduite red */
	reduite red(Mtemps, Nespace, T, L,K,sigma ,r);
	/** recuperation de la matrice A et de dS de calcul*/
	double * A_=red.get_A();
	double ds=red.get_ds();
	/**initialisation du vecteur solution*/
	double* sol=new double[Nespace];
	for ( int i=0;i<Nespace;i++){
		double s_i=(i+1)*ds;
		sol[i]=max(0.0, (exp(0.5*(k-1)*s_i)-exp(0.5*(k+1)*s_i)));//condition initiale put
	}
	for (int m=0;m<Mtemps;m++){
		/**Calcul du vecteur Fn */
		double*Cn=reduite::VecteurF_put(k,sigma,r,m+1, Nespace,L);
		double *b=SommeVecteur(sol,Cn,Nespace);
		/** resolution du systeme AUn+1=Un+ Fn */
		double* res=res_LU(A_,b,Nespace);
		/** recuperation du vecteur solution*/
		for ( int i=0;i<Nespace;i++){
			sol[i]=res[i];// colonne m+1 ligne i


		}
		/** suppresion des tableaux*/
		delete[]Cn;
		delete[]b;
		delete[]res;
	}
	return sol;
}
/** methode permettant de resoudre une edp de blackschole sachant que le payoff est un call
 *  ainsi que la discretisation en temps et en espace */

double * reduite::res_reduite_call(BlackScholes Bs, int Mtemps, int Nespace, double K){
	/**Recuperation de parametre de l'edp Bs*/
	double r=Bs.get_r();
	double L=Bs.get_L();//borne espace
	double T=Bs.get_T();
	double sigma=Bs.get_sigma();
	double k=2*r/(sigma*sigma);
	/** cr√©ation de l'objet Resolution_reduite red */
	reduite red(Mtemps, Nespace, T, L,K,sigma ,r);
	/** recuperation de la matrice A et de dS de calcul*/
	double * A_=red.get_A();
	double ds=red.get_ds();

	/**initialisation du vecteur solution*/
	double* sol=new double[Nespace];
	for ( int i=0;i<Nespace;i++){
		double s_i=(i+1)*ds;
		sol[i]=max(0.0, exp(0.5*(k+1)*s_i)-exp(0.5*(k-1)*s_i));//condition initiale call
	}
	for (int m=0;m<Mtemps;m++){
		/**Calcul du vecteur Fn */
		double*Cn=reduite::VecteurF_call(k,sigma,r,m+1, Nespace,L);
		double *b=SommeVecteur(sol,Cn,Nespace);
		/** resolution du systeme AUn+1=Un+ Fn */
		double* res=res_LU(A_,b,Nespace);
		/** recuperation du vecteur solution*/
		for ( int i=0;i<Nespace;i++){
			sol[i]=res[i];// colonne m+1 ligne i


		}
		/** suppresion des tableaux*/
		delete[]Cn;
		delete[]b;
		delete[]res;
	}
	std::cout<<"Solution réduite pour un call \n";
	std::cout<<"--------------------\n";
	return sol;
}

double * reduite::res_reduite(BlackScholes Bs, int Mtemps, int Nespace, double K){
	char a=(*(Bs.get_payoff())).Pay_Off::get_type();


	double* sol=new double[Nespace];
	//double alpha=(1-k)/2;
	if(a=='p'){
		sol=reduite::res_reduite_put( Bs,  Mtemps,  Nespace,  K);

	}
	else{
		sol=reduite::res_reduite_call( Bs,  Mtemps,  Nespace,  K);

	}

	return sol;
}

/**
 * methode permettant de definir les changements de variables inverses
 * pour passe de l'edp Resolution_reduite ‡ l'edp de black and scholes
 */
double * reduite::red_compl(double* sol, int Nespace, double K, double sigma , double r , int T){
	double * C= new double[Nespace];
	double k=2*r/(sigma*sigma);
	double alpha=(1-k)/2;
	double beta=-0.25*(k+1)*(k+1);
	double tau=(sigma*sigma)*T*0.5;
	for (int i =0;i<Nespace;i++){
		double x_i=(i+1)*dS;
		C[i]=K*exp(alpha*x_i+beta*tau)*sol[i];
	}
	//transformatio pour passer de la solution de la reduite a celle de la c
	affichage_vecteur(C,Nespace);
	return C;

}




}
