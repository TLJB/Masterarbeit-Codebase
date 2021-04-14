#ifndef KRONECKER_H
#define KRONECKER_H


// deal does not appear to have a native function for the kronecker 
// delta
double kronecker (int A,int B)
{
    if (A==B)
    { 
        return 1;
    }
    else
    {
        return 0;
    }
}

#endif 