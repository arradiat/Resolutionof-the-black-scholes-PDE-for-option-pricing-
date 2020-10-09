/**
 * Classe black-scholes de l'Èquation au derivÈ partielle de Black & Scholes
 * Les variables sont formÈes d'un pointeur vers le payoff en question ( call ou put)
 * et les paramËtres de l'Èquation differentielle ainsi que les bornes L et T
 *  Les mÈthodes sont constituÈes de l'ensemble des getters et du constructeur */

#include "Pay_off.h"
#include<cmath>
#include "general.h"
#ifndef BLACK_SCHOLES_H_
#define BLACK_SCHOLES_H_
namespace edp{

class BlackScholes {

	Pay_Off* pay_off_; /* Pointeur vers le payoff souhaitÈ (Call/Put)*/
	double r_;  /*taux d'intÈret*/
	double T_;  /*borne temporelle*/
	double L_;  /*borne spatiale*/
	double sigma_; /* VolatilitÈ*/


public:
	/* constructeur surchargÈ crÈant
	 * une Èquation avec les paramËtres*/
	BlackScholes(Pay_Off* pay_off,double r,double T,double L,double sigma);
	/**methode permettant de recupere le membre de donnees r_*/
	double get_r();
	/**methode permettant de recupere le membre de donnees T_*/
	double get_T();
	/**methode permettant de recupere le membre de donnees L_*/
	double get_L();
	/**methode permettant de recupere le membre de donnees sigma_*/
	double get_sigma();
	/**methode permettant de recupere le membre de donnees pay_off_*/
	Pay_Off* get_payoff();
};


}/// namespace
#endif /* BLACK_SCHOLES_H_ */
