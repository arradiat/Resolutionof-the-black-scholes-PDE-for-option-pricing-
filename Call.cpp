/*
 * Call.cpp
 *
 */
#include "general.h"
#include "black-scholes.h"
#include "Pay_off.h"
#include "Call.h"
#include<stdlib.h>
#include <cmath>
namespace edp{

/** Constructeur valuÈ de la classe*/
Call::Call( double K1):Pay_Off('c'){
	K=K1;
}

/**methode permettant de recupere le membre de donnees K_*/
double Call::get_K(){
	return K;
}
/*
 * Methode virtuelle pure retournant le payoff(call)
 */

double Call::operator() (const double& S) const {
	return max(S-K,0.0); // Standard European call pay-off
}

}
