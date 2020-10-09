/*
 * main.cpp
 *
 *  Created on: 30 déc. 2019
 *      Author: mamediarra
 */
#include <fstream>
#include "Pay_off.h"
#include "Put.h"
#include "Call.h"
#include "Resolution_complete.h"
#include "black-scholes.h"
#include<stdlib.h>
#include<ostream>
#include<iostream>
#include <cmath>
#include "general.h"
#include "Resolution_reduite.h"
using namespace edp;
using namespace std;
int main(){
	/**initialisation des paramétres de l'edp pour les tests*/
	int T=1;
	double r=0.1;
	double sigma=0.1;
	double K=100;
	int L=300;
	int M=960;
	int N=960;
	/**
	 * CAS OU LE PAYOFF EST UN PUT */


	/**creation d'un Put avec comme strike K et de l'edp  blackScholes */
	Put put(K);
	BlackScholes Bs(&put,r,T,L,sigma);

	/**creation d'un objet calcul ave les parametres de notre EDP*/
	Resolution_complete C(T, L,r,sigma, K, M, N);

	/**Resolution_complete du vecteur solution de l'equation compléte et affichage des etapes */
	double * sol1=C.res_edp(Bs,N,M,100);


	/** creation d'un objet Resolution_reduite correspondant a la resolution de l'équation black scholes*/
	reduite red(N,M,T,L,K,sigma ,r );

	/**Resolution_complete de la solution de l'equation Resolution_reduite et  affichage des etapes */
	double* sol2=red.res_reduite(Bs,M,N,K);


	/**Passage de la solution de l'equation réduite a la solution de l'equation compléte*/
	double* sol=red.red_compl(sol2,N,K,sigma,r,T);

	/** ecriture du resultat de l'equation réduite  dans le fichier solution_reduite.txt*/
	ofstream objetfichier;
	objetfichier.open("solution_reduite_put.txt", ios::out); //on ouvrre le fichier en ecriture
	if (objetfichier.bad()) //permet de tester si le fichier s'est ouvert sans probleme
		return 1;
	/** ecriture */
	for(int i =0;i<N;i++){
		objetfichier << sol[i] << endl;//*
	}
	objetfichier.close(); //on ferme le fichier pour liberer la mémoire

	/** ecriture du resultat de l'equation complete  dans le fichier solution_complete.txt*/
	ofstream objetfichier1;
	objetfichier1.open("solution_complete_put.txt", ios::out); //on ouvrre le fichier en ecriture
	if (objetfichier1.bad()) //permet de tester si le fichier s'est ouvert sans probleme
		return 1;
	for(int i =0;i<N;i++){
		objetfichier1 << sol1[i] << endl;//*
	}

	objetfichier1.close(); //on ferme le fichier pour liberer la mémoire

	/**
	 * CAS OU LE PAYOFF EST UN CALL */

	/**creation d'un Call avec comme strike K et de l'edp  blackScholes */

	Call call(K);
	BlackScholes Bs1(&call,r,T,L,sigma);

	/**creation d'un objet calcul ave les parametres de notre EDP*/
	Resolution_complete C1(T, L,r,sigma, K, M, N);

	/**Resolution_complete du vecteur solution de l'equation compléte et affichage des etapes */
	double * solC=C1.res_edp(Bs1,N,M,K);


	/** creation d'un objet Resolution_reduite correspondant a la resolution de l'équation black scholes*/
	reduite red1(N,M,T,L,K,sigma ,r );

	/**Resolution_complete de la solution de l'equation Resolution_reduite et  affichage des etapes */
	double* solC1=red1.res_reduite(Bs1,M,N,K);
	/**Passage de la solution de l'equation réduite a la solution de l'equation compléte*/
	double* solC2=red.red_compl(solC1,N,K,sigma,r,T);

	/** ecriture du resultat de l'equation réduite  dans le fichier solution_reduite.txt*/
	ofstream fichier;
	fichier.open("solution_reduite_call.txt", ios::out); //on ouvrre le fichier en ecriture
	if (objetfichier.bad()) //permet de tester si le fichier s'est ouvert sans probleme
		return 1;
	/** ecriture */
	for(int i =0;i<N;i++){
		fichier << solC2[i] << endl;//*
	}
	fichier.close(); //on ferme le fichier pour liberer la mémoire

	/** ecriture du resultat de l'equation complete  dans le fichier solution_complete.txt*/
	ofstream fichier1;
	fichier1.open("solution_complete_call.txt", ios::out); //on ouvrre le fichier en ecriture
	if (fichier1.bad()) //permet de tester si le fichier s'est ouvert sans probleme
		return 1;
	for(int i =0;i<N;i++){
		fichier1 << solC[i] << endl;//*
	}

	fichier1.close(); //on ferme le fichier pour libérer la mémoire
	/**
	 * Récuperation et ecriture de la matrice A de la résolution de l'edp compléte
	 *  dans un fichier txt */
	double * A= C.get_A();
	ofstream fichier_A;
	fichier_A.open("Matrice_A_comp.txt", ios::out); //on ouvrre le fichier en ecriture
	if (fichier_A.bad()) //permet de tester si le fichier s'est ouvert sans probleme
		return 1;
	/** ecriture */
	for(int i =0;i<N;i++){
		for (int j=0; j<N; j++){
			fichier_A << A[i*N+j]<<";";//*
		}
		fichier_A << endl;
	}
	fichier_A.close(); //on ferme le fichier pour liberer la mémoire
	/**
	 * Récuperation et ecriture de la matrice A de la résolution de l'edp reduite
	 *  dans un fichier txt */
	double * A1=red.get_A();
	ofstream fichier_A1;
	fichier_A1.open("Matrice_A_reduite.txt", ios::out); //on ouvrre le fichier en ecriture
	if (fichier_A1.bad()) //permet de tester si le fichier s'est ouvert sans probleme
		return 1;
	/** ecriture */
	for(int i =0;i<N;i++){
		for (int j=0; j<N; j++){
			fichier_A1 << A1[i*N+j]<<";";//*
		}
		fichier_A1 << endl;
	}
	fichier_A1.close(); //on ferme le fichier pour liberer la mémoire
	/**
	 * Récuperation et ecriture de l'erreur relative  entre les deux résolutions */
	ofstream erreur;
	erreur.open("erreur_Put.txt", ios::out); //on ouvrre le fichier en ecriture
	if (erreur.bad()) //permet de tester si le fichier s'est ouvert sans probleme
		return 1;
	/** ecriture */
	for(int i =0;i<N;i++){
			erreur << erreur_relative(sol,sol1,N)[i] <<endl;//*


	}
	erreur.close(); //on ferme le fichier pour liberer la mémoire

	/**
		 * Récuperation et ecriture de l'erreur relative  entre les deux résolutions */
		ofstream erreur1;
		erreur1.open("erreur_call.txt", ios::out); //on ouvrre le fichier en ecriture
		if (erreur1.bad()) //permet de tester si le fichier s'est ouvert sans probleme
			return 1;
		/** ecriture */
		for(int i =0;i<N;i++){
				erreur1 << erreur_relative(solC,solC2,N)[i] <<endl ;//*

		}
		erreur.close(); //on ferme le fichier pour liberer la mémoire

	return 0;
}



