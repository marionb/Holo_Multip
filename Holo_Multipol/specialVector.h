#ifndef SPECIALVECTOR_H
#define SPECIALVECTOR_H

#include <vector>

template<class Type>
struct SpecialVector
{
    std::vector<Type> positive;
    std::vector<Type> negative;
};
#endif
