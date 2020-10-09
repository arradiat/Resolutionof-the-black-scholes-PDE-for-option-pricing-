/*
 * Put.cpp
 *
 *  Created on: 30 d√©c. 2019
 *      Author: mamediarra
 */
#include"general.h"
#include "black-scholes.h"
#include "Pay_off.h"
#include "Put.h"
#include<stdlib.h>
#include <cmath>


namespace edp{

/** Constructeur valuÈ de la classe*/
Put::Put( double K1):Pay_Off('p'){
	K=K1;
}

/**methode permettant de recupere le membre de donnees K_*/
double Put::get_K(){
	return K;
}

/*
 * Methode virtuelle pure retournant le payoff(put)
 */
double Put::operator() (const double& S) const {
	double k=K;
	return max(k-S, 0.0);
}


}


