/*
 * reduite.h
 *
 *  Created on: 3 janv. 2020
 *      Author: mame diarra
 */

#ifndef REDUITE_H_
#define REDUITE_H_
#include "black-scholes.h"
#include <cmath>
#include <iostream>
#include "Pay_off.h"
#include "Put.h"
#include "Call.h"
#include <stdlib.h>
#include "general.h"
namespace edp{
class reduite{
	double dT; /** pas en temps de la resolution */
	double dS; /** pas en espace de la resolution */
	int T; /** borne temporelle*/
	int L;/** borne spatiale*/
	double * A;/** matrice A de la resolution  */
public:
	/**methode permettant de recupere le membre de donnes A*/
	double* get_A();
	/**methode permettant de recupere le membre de donnes dT*/
	double get_dt();
	/**methode permettant de recupere le membre de donnes dS*/
	double get_ds();
	/** Destrueteur de la classe */
	~reduite();
	/** Constructeur valuÈ de la classe*/
	reduite(int M, int N, double  T, double L, int K,double sigma, double r);
	/**
	 * methode permettant de calculer le vecteur Fn de la r√©solution
	 *  dans le cas ou le pay_off est un PUt*/
	double* VecteurF_put(double K,double sigma,double r,int n,int N,int L);
	/**
	 * methode permettant de calculer le vecteur Fn de la r√©solution
	 *  dans le cas ou le pay_off est un call*/
	double* VecteurF_call(double K,double sigma,double r,int n,int N,int L);
	/** methode permettant de resoudre une edp Resolution_reduite  sachant que le payoff est un put
	 *  ainsi que la discretisation en temps et en espace */
	double * res_reduite_put(BlackScholes Bs, int Mtemps, int Nespace, double K);
	/** methode permettant de resoudre une edp Resolution_reduite  sachant que le payoff est un call
	 *  ainsi que la discretisation en temps et en espace */
	double * res_reduite_call(BlackScholes Bs, int Mtemps, int Nespace, double K);
	/**
	 * methode permettant de resoudre une edp Resolution_reduite connaissant le
	 * strike du payoff ainsi que la discretisation en temps et en espace,
	 * LA methode affiche le resultat ainsi que les etapes de la resolution  */
	double * res_reduite(BlackScholes Bs,int M, int N,double K);
	/**
	 * methode permettant de definir les changements de variables inverses
	 * pour passe de l'edp Resolution_reduite ‡ l'edp de black and scholes
	 */
	double *red_compl(double * sol, int Nespace, double K, double sigma, double r , int T );


};
}





#endif /* REDUITE_H_ */
