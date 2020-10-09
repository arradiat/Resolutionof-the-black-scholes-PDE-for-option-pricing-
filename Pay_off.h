/**
 * Pay_off.h
 *
 * classe abstraite Payoff composÈe d'un membre
 * de donnÈes type_ qui  reprÈsente le type de l'option (put ou call) , d'un
 * constructeur ,d'un desructeur par dÈfaut ,d'un getter qui retourne le type de l'option et enfin
 * de  l'operateur () qui est un functeur
 * (c'est a dire fonction d'objet) ce qui veut dire qu'on peut appeler l'objet comme on
 * on peut appeler une  fonction. L'appel du foncteur calcule la valeur du Payoff et la renvoie
 */


#ifndef PAY_OFF_H_
#define PAY_OFF_H_
namespace edp{

class Pay_Off {

	char type_; /* Type de l'option (payoff): Call/Put*/

public:

	/*Constructeur valuÈ de la classe*/
	Pay_Off(char  a);


	/*Destructeur par dÈfaut*/
	virtual ~Pay_Off() {}; /** destructeur */

	/**methode permettant de recupere le membre de donnees type_*/

	char get_type();



	/** Operateur () prend comme argument
	 * la valeur de l'actif sous jacent S.
	 * Methode virtuelle pure*/

	virtual double operator() (const double& S) const = 0;


};
} ///namespace
#endif /* PAY_OFF_H_ */

