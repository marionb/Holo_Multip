#ifndef SPECIALVECTOR_H
#define SPECIALVECTOR_H

#include <vector>
#include <cmath>
#include <assert.h>

template<class Type>
struct SpecialVector
{
    std::vector<Type> positive;
    std::vector<Type> negative;
    int size;

    SpecialVector()
    {
        size=0;
    }

    /**
    *resize the positive and negative vector and set the new size
    *@param max int size of expansion
    */
    void expand(int max)
    {
        positive.resize(max+1);
        negative.resize(max);
        size=negative.size()+positive.size();
    }

    /**
    *add an elemenent at a certain position
    *@param elem    the element to be added
    *@param pos     int the position where to add the element; can be positive as well as negative but has to be within the bounds of the vector.
    */
    void addElement(Type elem, int pos)
    {
        if(size/2<std::abs(pos))
        {
                std::cout<<"can not place element; vector to small!\n";
                return;
        }
        if(pos<0)
        {
            negative[-pos-1]=elem; //the position element at -1 is placed in the position 0, -2 at -1 ... because there is no 0
        }else //pos is positive
        {
            positive[pos]=elem;
        }

    }

    /**
    *get an element at a certain position
    *@param pos     integer the position of the wanted element; can be positive or negative
    */
    Type getElement(int pos)
    {
        assert(size/2<std::abs(pos));

        if (pos<0)
        {
            return negative[-pos-1]; //the position element at -1 is placed in the position 0, -2 at -1 ... because there is no 0
        }else //pos is positive
        {
            return positive[pos];
        }
    }
};
#endif
