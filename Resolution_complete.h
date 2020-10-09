
/**
 *
 *Classe Calcul.h implémentant la resolution de l'edp de blackscholes .
 *Classe. Elle a comme mebre de données le pas en temps dT, le pas
 *Classe. en espace dS et les matrice A, B de la resolution
 *
 *      Author: imane
 */

#ifndef CALCUL_H_
#define CALCUL_H_
#include "black-scholes.h"
#include <cmath>
#include <iostream>
#include "general.h"

namespace edp{

class Resolution_complete{

private:
	double dT;/** pas en temps de la resolution */

	double dS;/** pas en espace  de la resolution */

	double * A; /** matrice A de la resolution  */

	double * B;/** matrice B de la resolution  */

public:

	/**methode permettant de recupere le membre de donnes A*/
	double* get_A();
	/**methode permettant de recupere le membre de donnes B*/
	double* get_B();
	/**methode permettant de recupere le membre de donnes dT*/
	double get_dt();
	/**methode permettant de recupere le membre de donnes dT*/
	double get_ds();

	/** Constructeur valué de la classe*/
	Resolution_complete(double T, double S,double r,double sigma,double K, int N, int M);

	/** Destrueteur de la classe */
	~Resolution_complete();

	/**
	 * methode permettant de calculer le veceteur Fn de la résolution
	 *  dans le cas ou le pay_off est un PUt*/

	double* VecteurF_put(double K,double sigma,double r,int n,int N);

	/**
	 * methode permettant de calculer le veceteur Fn de la résolution
	 * dans le cas ou le pay_off est un call */
	double* VecteurF_call(double K,double sigma,double r,int n,int N);

	/** methode permettant de resoudre une edp de blackschole sachant que le payoff est un put
	 *  ainsi que la discretisation en temps et en espace */
	double* res_edp_put(BlackScholes Bs,int M,int N,double K);

	/** methode permettant de resoudre une edp de blackschole sachant que le payoff est un call
	 *  ainsi que la discretisation en temps et en espace */
	double* res_edp_call(BlackScholes Bs,int M,int N,double K);

	/**
	 * methode permettant de resoudre une edp de blackschole connaissant le
	 * strike du payoff ainsi que la discretisation en temps et en espace,
	 * LA methode affiche le resultat ainsi que les etapes de la resolution  */
	double * res_edp(BlackScholes Bs, int M, int N,double K);


};

}



#endif /* CALCUL_H_ */
