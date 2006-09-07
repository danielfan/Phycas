#ifndef CIPRES_HELPER_H__
#define CIPRES_HELPER_H__
#include <iostream>
#include "ConfigDependentHeaders.h"
#include <CipresIDL/api1/CipresS.h>

std::ostream & operator << (std::ostream &outs, const CipresIDL_api1::DataMatrix &matrix);
std::ostream & operator << (std::ostream &outs, const CipresIDL_api1::RawMatrix &matrix);
std::ostream & operator << (std::ostream &outs, const CipresIDL_api1::Characters &characters);
std::ostream & operator << (std::ostream &outs, const CipresIDL_api1::Tree  &tree);
std::ostream & operator << (std::ostream &outs, const CipresIDL_api1::TaxonSeqSeq  &taxonSubsets);
std::ostream & operator << (std::ostream &outs, const CipresIDL_api1::TaxonSeq  &taxonSeq);
std::ostream & operator << (std::ostream &outs, const CipresIDL_api1::CharStateLookup &matrix);
std::ostream & operator << (std::ostream &outs, const CipresIDL_api1::StateSet &stateSet);




#endif // CIPRES_HELPER_H__
