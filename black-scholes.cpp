/*
 * black-scholes.cpp
 *
 */
#include "general.h"
#include "Pay_off.h"
#include "Put.h"
#include "Call.h"
#include "black-scholes.h"
#include<stdlib.h>
#include <cmath>
namespace edp{

/** Constructeur valuÈ de la classe*/
BlackScholes::BlackScholes(Pay_Off* pay_off,double r,double T,double L,double sigma){

	pay_off_=pay_off;
	r_=r;
	T_=T;
	L_=L;
	sigma_=sigma;

}

/**methode permettant de recupere le membre de donnees r_*/

double BlackScholes::get_r(){
	return(r_);
}
/**methode permettant de recupere le membre de donnees T_*/
double BlackScholes::get_T(){
	return(T_);
}
/**methode permettant de recupere le membre de donnees L_*/
double BlackScholes::get_L(){
	return(L_);
}
/**methode permettant de recupere le membre de donnees sigma_*/
double BlackScholes::get_sigma(){
	return(sigma_);
}
/**methode permettant de recupere le membre de donnees payoff_*/
Pay_Off* BlackScholes::get_payoff(){
	return(pay_off_);
}


}


