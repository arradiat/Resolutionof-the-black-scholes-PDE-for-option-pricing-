/*
 * Pay_off.cpp
 *
 *  Created on: 29 d√©c. 2019
 *      Author: mamediarra
 */

#include "Pay_off.h"
#include<stdlib.h>
#include <cmath>
namespace edp{

/*Constructeur surchargÈ */

Pay_Off::Pay_Off(char  a){
	type_=a;

}

/**methode permettant de recupere le membre de donnees type_*/
char Pay_Off::get_type(){
	return (type_);
}



}


