/*
 * Call.h
 *
 *Classe Call(classe fille) hÈritant de la class mËre Payoff.
 *   ComposÈe de la variable K(strike)
 *   DÈfinie un constructeur, destructeur par dÈfaut,
 *   un getter et la methode virtuelle pure () de la classe mËre.

 */
#include "general.h"
#include "Pay_off.h"
#include "black-scholes.h"
#ifndef CALL_H_
#define CALL_H_
namespace edp{

class Call : public Pay_Off {

private:

	double K; /* Strike */

public:

	/**constructeur surchargÈ crÈant un call*/
	Call( double K);

	/** destructeur par dÈfaut */
	virtual ~Call() {};

	/**methode permettant de recupere le membre de donnees K_*/
	double  get_K();

	/**
	 * Implementation de la mÈthode virtuelle pure () de la classe mËre */
	virtual double operator() (const double& S ) const;
};

} ///namespace

#endif /* CALL_H_ */
