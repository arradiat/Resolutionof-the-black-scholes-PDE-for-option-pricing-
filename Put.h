/*
 * Put.h
 *
 *Classe Put(classe fille) hÈritant de la class mËre Payoff.
 *   ComposÈe de la variable K(strike)
 *   DÈfinie un constructeur, destructeur par dÈfaut,
 *   un getter et la methode virtuelle pure () de la classe mËre.

 */
#include "black-scholes.h"
#include "Pay_off.h"
#include "general.h"
#ifndef PUT_H_
#define PUT_H_
namespace edp{

class Put : public Pay_Off {

private:

	double K; /* Strike */

public:

	/**constructeur surchargÈ crÈant un put*/
	Put( double K);

	/** destructeur par dÈfaut */
	virtual ~Put() {};

	/**methode permettant de recupere le membre de donnees K_*/
	double get_K();

	/**
	 * Implementation de la mÈthode virtuelle pure () de la classe mËre */
	virtual double operator() (const double& S) const;
};
} ///namespace
#endif /* PUT_H_ */
