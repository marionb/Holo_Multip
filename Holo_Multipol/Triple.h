#ifndef DOMEGA_H
#define DOMEGA_H


template <class T1, class T2, class T3>
struct Triple
{    //Structure of three accociated elements of different types -> making life easier
    T1 first;
    T2 second;
    T3 third;
    /**
    * Constructor with given element
    */
    Triple(T1 e1, T2 e2, T3 e3):first(e1),second(e2),third(e3)
    {}

    /**
    * Default consstructor
    */
    Triple()
    {}
};

#endif
