#ifndef SGT_H
#define SGT_H

namespace sgt
{

typedef long long int LongInt;

template<class T>
void releaseMemory(std::vector<T> *v)
{
    if(v)
    {
       v->clear();
       v->shrink_to_fit();
       delete v;
    }
    v = nullptr;
}
/**
 * @brief diff: Calculate the n-th discrete difference along the given axis.
 * @param S: signal
 * @return : n-th discrete difference
 */
template<class T>
std::vector<T> * diff(typename std::vector<T>::iterator vBegin, typename std::vector<T>::iterator vEnd, int dfRng = 1, bool absValues = false)
{
    Q_UNUSED(dfRng)
    std::vector<T> * diffRslt = new std::vector<T>();

    for(typename std::vector<T>::iterator sIt = vBegin; sIt != (vEnd - 1); ++sIt)
        diffRslt->push_back(absValues ? std::abs((*(sIt+1)) - (*(sIt))) : (*(sIt+1)) - (*(sIt)));


    return diffRslt;
}

template<class T>
double diffOfPoint(typename std::vector<T>::iterator itStart, typename std::vector<T>::iterator point, typename std::vector<T>::iterator itEnd, int dfRange = 1, bool absValues = false)
{
    double diffOfPointRslt = 0;

    if(dfRange == 0)
        return *point;

    typename std::vector<T>::iterator it = point;

    for(++it; it != itEnd && it < point + dfRange; ++it)
        diffOfPointRslt += absValues ? std::abs(*it) : *it;

    if(point - itStart >= dfRange)
    {
        for(it = point - dfRange; it < point - 1; ++it)
            diffOfPointRslt -= absValues ? std::abs(*it) : *it;
    }

    return diffOfPointRslt / dfRange;
}

template<class T>
std::vector<T> * diffRang(typename std::vector<T>::iterator itStart, typename std::vector<T>::iterator itEnd, int dfRange = 1, bool absValues = false)
{
    std::vector<T> * diffRslt = new std::vector<T>();

    for(std::vector<T>::iterator currIt = itStart; currIt != itEnd - 1; ++currIt)
        diffRslt->push_back(static_cast<T>(diffOfPoint<T>(itStart, currIt, itEnd, dfRange, absValues)));


    return diffRslt;
}

template<class T>
std::vector<T> * diff(typename std::vector<T>::reverse_iterator vBegin, typename std::vector<T>::reverse_iterator vEnd, int dfRng = 1, bool absValues = false)
{
    Q_UNUSED(dfRng)
    std::vector<T> * diffRslt = new std::vector<T>();

    for(typename std::vector<T>::reverse_iterator sIt = vBegin; sIt != (vEnd + 1); ++sIt)
        diffRslt->push_back(absValues ? std::abs((*(sIt+1)) - (*(sIt))) : (*(sIt+1)) - (*(sIt)));


    return diffRslt;
}

template<class T>
std::vector<T> * diff(std::vector<T> *S, bool absValues = false)
{
    std::vector<T> * diffRslt = new std::vector<T>();

    for(std::vector<T>::iterator sIt = S->begin(); sIt != (S->end() - 1); ++sIt)
        diffRslt->push_back(absValues ? std::abs((*(sIt+1)) - (*(sIt))) :(*(sIt+1)) - (*(sIt)));

    return diffRslt;
}

template<class T>
double diff(typename std::vector<T>::iterator vBegin, typename std::vector<T>::iterator vEnd, quint64 point, int dfRng = 1)
{
    double output = 0;
    int size = vEnd - vBegin;

    if(dfRng == 0)
        return (*vBegin);

    for(size_t i = point + 1; i <= point + dfRng && i < size; ++i)
        output += *(vBegin + i);

    for(size_t i2 = point - dfRng ; i2 <= point - 1; ++i2)
        output -= *(vBegin - i2);

    return output / dfRng;
}

template<class T>
std::vector<LongInt> * indexes_of_min(typename std::vector<T>::iterator vBegin, typename std::vector<T>::iterator vEnd, bool absValues = false)
{
    std::vector<LongInt> * indexes = nullptr;
    std::vector<T> *d = sgt::diff<T>(vBegin, vEnd, absValues);

    if(d && d->size() > 0)
    {
        std::vector<T> d1(d->begin(), d->end() - 1);
        std::vector<T> d2(d->begin() + 1, d->end());

        std::vector<T>::iterator d1It = d1.begin();
        std::vector<T>::iterator d2It = d2.begin();
        std::vector<T>::iterator itEnd = d1.end();

        LongInt index = 0;
        indexes = new std::vector<LongInt>();

        for(; d1It != itEnd; ++d1It, ++d2It)
        {
            if((*d1It)*(*d2It) < 0 && (*d1It) < 0)
                indexes->push_back(index + 1);
            ++index;
        }

    }

    delete d;

    if(indexes && indexes->size() == 0)
    {
        delete indexes;
        indexes = nullptr;
    }

    return indexes;
}

template<class T>
std::vector<LongInt> * indexes_of_min(typename std::vector<T>::reverse_iterator vBegin, typename std::vector<T>::reverse_iterator vEnd, bool absValues = false)
{
     std::vector<LongInt> * indexes = nullptr;
    std::vector<T> *d = sgt::diff<T>(vBegin, vEnd, absValues);
    if(d && d->size() > 0)
    {
        std::vector<T> d1(d->begin(), d->end() - 1);
        std::vector<T> d2(d->begin() + 1, d->end());

        std::vector<T>::iterator d1It = d1.begin();
        std::vector<T>::iterator d2It = d2.begin();
        std::vector<T>::iterator itEnd = d1.end();

        LongInt index = 0;
       indexes = new std::vector<LongInt>();

        for(; d1It != itEnd; ++d1It, ++d2It)
        {
            if((*d1It)*(*d2It) < 0 && (*d1It) < 0)
                indexes->push_back(index + 1);
            ++index;
        }
    }

    delete d;

    if(indexes && indexes->size() == 0)
    {
        delete indexes;
        indexes =nullptr;
    }

    return indexes;
}

template<class T>
std::vector<LongInt> * indexes_of_min(std::vector<T> *S, bool absValues = false)
{
    std::vector<LongInt> * indexes = nullptr;
    std::vector<T> *d = sgt::diff<T>(S, absValues);
    if(d && d->size() > 0)
    {
        std::vector<T> d1(d->begin(), d->end() - 1);
        std::vector<T> d2(d->begin() + 1, d->end());

        std::vector<T>::iterator d1It = d1.begin();
        std::vector<T>::iterator d2It = d2.begin();
        std::vector<T>::iterator itEnd = d1.end();

        LongInt index = 0;
        indexes = new std::vector<LongInt>();

        for(; d1It != itEnd; ++d1It, ++d2It)
        {
            if((*d1It)*(*d2It) < 0 && (*d1It) < 0)
                indexes->push_back(index + 1);
            ++index;
        }
    }

    delete d;

    if(indexes && indexes->size() == 0)
    {
        delete indexes;
        indexes =nullptr;
    }

    return indexes;
}

template<class T>
std::vector<LongInt> * indexes_of_max(typename std::vector<T>::iterator vBegin, typename std::vector<T>::iterator vEnd, bool absValues = false)
{
    std::vector<LongInt> * indexes = nullptr;
    std::vector<T> *d = sgt::diff<T>(vBegin, vEnd, absValues);

    if(d && d->size() > 0)
    {
        std::vector<T> d1(d->begin(), d->end() - 1);
        std::vector<T> d2(d->begin() + 1, d->end());

        std::vector<T>::iterator d1It = d1.begin();
        std::vector<T>::iterator d2It = d2.begin();
        std::vector<T>::iterator itEnd = d1.end();

        LongInt index = 0;
        indexes = new std::vector<LongInt>();

        for(; d1It != itEnd; ++d1It, ++d2It)
        {
            if((*d1It)*(*d2It) < 0 && (*d1It) > 0)
                indexes->push_back(index + 1);
            ++index;
        }
    }

    delete d;

    if(indexes && indexes->size() == 0)
    {
        delete indexes;
        indexes =nullptr;
    }

    return indexes;
}

template<class T>
std::vector<LongInt> * indexes_of_max(typename std::vector<T>::reverse_iterator vBegin, typename std::vector<T>::reverse_iterator vEnd, bool absValues = false)
{
    std::vector<LongInt> * indexes = nullptr;
    std::vector<T> *d = sgt::diff<T>(vBegin, vEnd, absValues);
    if(d && d->size() > 0)
    {
        std::vector<T> d1(d->begin(), d->end() - 1);
        std::vector<T> d2(d->begin() + 1, d->end());

        std::vector<T>::iterator d1It = d1.begin();
        std::vector<T>::iterator d2It = d2.begin();
        std::vector<T>::iterator itEnd = d1.end();

        LongInt index = 0;
        indexes = new std::vector<LongInt>();

        for(; d1It != itEnd; ++d1It, ++d2It)
        {
            if((*d1It)*(*d2It) < 0 && (*d1It) > 0)
                indexes->push_back(index + 1);
            ++index;
        }
    }

    delete d;

    if(indexes && indexes->size() == 0)
    {
        delete indexes;
        indexes =nullptr;
    }

    return indexes;
}

template<class T>
std::vector<LongInt> * indexes_of_max(std::vector<T> *S, bool absValues = false)
{
    std::vector<LongInt> * indexes = nullptr;
    std::vector<T> *d = sgt::diff<T>(S, absValues);
    if(d && d->size())
    {
        std::vector<T> d1(d->begin(), d->end() - 1);
        std::vector<T> d2(d->begin() + 1, d->end());

        std::vector<T>::iterator d1It = d1.begin();
        std::vector<T>::iterator d2It = d2.begin();
        std::vector<T>::iterator itEnd = d1.end();

        LongInt index = 0;
        indexes = new std::vector<LongInt>();

        for(; d1It != itEnd; ++d1It, ++d2It)
        {
            if((*d1It)*(*d2It) < 0 && (*d1It) > 0)
                indexes->push_back(index + 1);
            ++index;
        }
    }

    delete d;

    if(indexes && indexes->size() == 0)
    {
        delete indexes;
        indexes =nullptr;
    }

    return indexes;
}

template<class T>
std::vector<LongInt> * min_max_indexes(typename std::vector<T>::iterator vBegin, typename std::vector<T>::iterator vEnd)
{
    std::vector<LongInt> * indexes = nullptr;
    std::vector<T> *d = sgt::diff<T>(vBegin, vEnd);

    if(d && d->size() > 0)
    {
        std::vector<T>::iterator dIt = d->begin();
        std::vector<T>::iterator itEnd = d->end();

        LongInt index = 0;
        LongInt lidx_min = 0, lidx_max = 0;
        indexes = new std::vector<LongInt>();

        T lmin = *dIt, lmax = *dIt;
        ++dIt;

        for(; dIt != itEnd; ++dIt)
        {
            ++index;
            if(*dIt < lmin)
            {
                lmin = *dIt;
                lidx_min = index;
            }

            if(*dIt > lmax)
            {
                lmax = *dIt;
                lidx_max = index;
            }
        }

        indexes->push_back(lidx_min);
        indexes->push_back(lidx_max);
    }

    delete d;

    return indexes;
}

template<class T>
std::vector<LongInt> *nonzeroIndexes(std::vector<T> *S)
{
    std::vector<T> S1(S->begin(), S->end() - 1);
    std::vector<T> S2(S->begin() + 1, S->end());

    std::vector<T>::iterator s1It = S1.begin();
    std::vector<T>::iterator s2It = S2.begin();
    std::vector<T>::iterator itEnd = S1.end();

    LongInt index = 0;
    std::vector<LongInt> * indexes = new std::vector<LongInt>();

    for(;s1It != itEnd; ++s1It, ++s2It)
    {
        if((*s1It) * (*s2It) < 0)
            indexes->push_back(index);
        ++index;
    }

    if(indexes->size() == 0)
    {
        delete indexes;
        indexes =nullptr;
    }

    return indexes;
}

template<class T>
std::vector<LongInt> *nonzeroIndexes(typename std::vector<T>::iterator vBegin, typename std::vector<T>::iterator vEnd)
{
    std::vector<T> S1(vBegin, vEnd - 1);
    std::vector<T> S2(vBegin + 1, vEnd);

    std::vector<T>::iterator s1It = S1.begin();
    std::vector<T>::iterator s2It = S2.begin();
    std::vector<T>::iterator itEnd = S1.end();

    LongInt index = 0;
    std::vector<LongInt> * indexes = new std::vector<LongInt>();

    for(;s1It != itEnd; ++s1It, ++s2It)
    {
        if((*s1It) * (*s2It) < 0)
            indexes->push_back(index);
        ++index;
    }

    if(indexes->size() == 0)
    {
        delete indexes;
        indexes =nullptr;
    }

    return indexes;
}

template<class T>
std::vector<LongInt> *nonzeroIndexes(typename std::vector<T>::reverse_iterator vBegin, typename std::vector<T>::reverse_iterator vEnd)
{
    std::vector<T> S1;//(vBegin, vEnd - 1);
    std::vector<T> S2;//(vBegin + 1, vEnd);

    copy<T>(vBegin, vEnd-1, &S1);
    copy<T>(vBegin+1, vEnd, &S2);

    std::vector<T>::iterator s1It = S1.begin();
    std::vector<T>::iterator s2It = S2.begin();
    std::vector<T>::iterator itEnd = S1.end();

    LongInt index = 0;
    std::vector<LongInt> * indexes = new std::vector<LongInt>();

    for(;s1It != itEnd; ++s1It, ++s2It)
    {
        if((*s1It) * (*s2It) < 0)
            indexes->push_back(index);
        ++index;
    }

    if(indexes->size() == 0)
    {
        delete indexes;
        indexes =nullptr;
    }

    return indexes;
}

template<class T>
std::vector<sgt::LongInt> * not_duplicate(std::vector<T> *S)
{
    sgt::LongInt N = S->size();
    std::vector<sgt::LongInt> *idx = new std::vector<sgt::LongInt>();
    std::vector<float> S1(N-2), S2(N-2), S3(N-2), dup(N-2);

    sgt::copy(S->begin()+1, S->end()-1, &S1);
    sgt::copy(S->begin(), S->end()-2, &S2);
    sgt::copy(S->begin()+2, S->end(), &S3);

    for(sgt::LongInt i = 0; i < N-2; ++i)
        dup[i] = S1[i] == S2[i] && S1[i] == S3[i] ? 1 : 0;

    idx->push_back(0);

    sgt::LongInt index = 0;
    for(sgt::LongInt i = 1; i < N-1; ++i)
    {
        ++index;
        if(dup[i-1] == 0)
            idx->push_back(index);
    }

    idx->push_back(N-1);

    return idx;
}

template<class T>
std::vector<sgt::LongInt> * not_duplicate(typename std::vector<T>::iterator itStart, typename std::vector<T>::iterator itEnd)
{
    sgt::LongInt N = itEnd - itStart;
    std::vector<sgt::LongInt> *idx = new std::vector<sgt::LongInt>();
    std::vector<T> S1(N-2), S2(N-2), S3(N-2), dup(N-2);

    sgt::copy(itStart+1, itEnd-1, &S1);
    sgt::copy(itStart, itEnd-2, &S2);
    sgt::copy(itStart+2, itEnd, &S3);

    for(sgt::LongInt i = 0; i < N-2; ++i)
        dup[i] = S1[i] == S2[i] && S1[i] == S3[i] ? 1 : 0;

    idx->push_back(0);

    sgt::LongInt index = 0;
    for(sgt::LongInt i = 1; i < N-1; ++i)
    {
        ++index;
        if(dup[i-1] == 0)
            idx->push_back(index);
    }

    idx->push_back(N-1);

    return idx;
}

template<class T>
LongInt zero_crossing(typename std::vector<T>::iterator vBegin, typename std::vector<T>::iterator vEnd)
{
    LongInt index = -1;
    typename std::vector<T>::iterator it = vBegin + 1;
    for(;it != vEnd; ++it)
    {
        if(((*(it-1)) > 0 && (*it <=0)) || ((*(it-1) < 0) && (*it >= 0)))
        {
            index = it - vBegin;
            break;
        }
    }
    return index;
}

template<class T>
LongInt zero_crossing(typename std::vector<T>::reverse_iterator rBegin, typename std::vector<T>::reverse_iterator rEnd)
{
    LongInt index = -1;
    typename std::vector<T>::reverse_iterator it = rBegin + 1;
    for(;it != rEnd; ++it)
    {
        if(((*(it-1)) > 0 && (*it <=0)) || ((*(it-1) < 0) && (*it >= 0)))
        {
            index = it - rBegin;
            break;
        }
    }
    return index;
}

template<class T>
std::vector<LongInt> * zero_crossing_all(typename std::vector<T>::iterator vBegin, typename std::vector<T>::iterator vEnd)
{
    std::vector<LongInt> * indexes = new std::vector<LongInt>();
    LongInt index = -1;
    typename std::vector<T>::iterator it = vBegin + 1 < vEnd ? vBegin + 1 : vBegin;
    for(;it != vEnd; ++it)
    {
//        qDebug() << *it;
        if(((*(it-1)) > 0 && (*it <=0)) || ((*(it-1) < 0) && (*it >= 0)))
        {
            index = it - vBegin;
            indexes->push_back(index);
        }
    }

    if(indexes->size() == 0)
    {
        delete indexes;
        indexes =nullptr;
    }

    return indexes;
}

/**
 * @brief sum: compute the sume of two vector they have to have the at least n elements z[i] = x[i] +(fact) * y[i]
 * @param xStart: iterator to the begining of vector x
 * @param yStart: iterator to the begining of vector y
 * @param n: the amount of elements that has to be sum
 * @param zStart: the vector where the result will be stored
 * @param fact: this is a cosntant value that could be 1 or -1 or even other necesary value
 */
template<class T>
void sum(typename std::vector<T>::iterator xStart, typename std::vector<T>::iterator yStart, sgt::LongInt n, typename std::vector<T>::iterator zStart, int fact=1)
{
    std::vector<T>::iterator end = xStart + n;
    for(;xStart != end; ++xStart, ++yStart, ++zStart)
    {
        (*zStart) = (*xStart) + fact * (*yStart);
    }
}

/**
 * @brief indexes_of: find all the indexes of the values that satisfy the predicate
 * @param xStart: iterator to the begining of vector x
 * @param xEnd: iterator to the end of vector x
 * @param pred: unitary predicate
 */
template<class T, class UnaryPredicate>
std::vector<LongInt> * indexes_of(typename std::vector<T>::iterator itStart, typename std::vector<T>::iterator itEnd, UnaryPredicate pred, int offset = 0)
{
    sgt::LongInt index = 0;
    std::vector<sgt::LongInt> *indexes = new std::vector<sgt::LongInt>();

    for(;itStart != itEnd; ++itStart)
    {
        if(pred(*itStart))
            indexes->push_back(index + offset);
        ++index;
    }

    if(indexes->size() == 0)
    {
        delete indexes;
        indexes = nullptr;
    }
    return indexes;
}

template<class T, class UnaryPredicate>
std::vector<LongInt> * indexes_of(typename std::vector<T>::reverse_iterator itStart, typename std::vector<T>::reverse_iterator itEnd, UnaryPredicate pred, int offset = 0)
{
    sgt::LongInt index = 0;
    std::vector<sgt::LongInt> *indexes = new std::vector<sgt::LongInt>();

    for(;itStart != itEnd; ++itStart)
    {
        if(pred(*itStart))
            indexes->push_back(index + offset);
        ++index;
    }

    if(indexes->size() == 0)
    {
        delete indexes;
        indexes =nullptr;
    }
    return indexes;
}

template<class T>
std::vector<LongInt> * indexes_of_first_show_up(std::vector<T> *S)
{
    std::vector<sgt::LongInt> *indexes = new std::vector<sgt::LongInt>();

    std::vector<T>::iterator itStart = S->begin();
    std::vector<T>::iterator itEnd = S->end();
    std::vector<T>::reverse_iterator init = S->rend();
    bool found;
    sgt::LongInt index = 0;
    indexes->push_back(index);
    ++itStart;

    for(;itStart != itEnd; ++itStart)
    {
        ++index;
        found = false;
        for(std::vector<T>::reverse_iterator itBackward = S->rend() - index; itBackward != init; ++itBackward)
        {
            if(*itBackward == *itStart)
            {
                found = true;
                break;
            }
        }
        if(!found)
        indexes->push_back(index);
    }

    if(indexes->size() == 0)
    {
        delete indexes;
        indexes =nullptr;
    }

    return indexes;
}

/**
 * @brief who_satisfy_this: return a vector where each position have 1 if the element satisfy the predicate or 0 in other case
 * @param xStart: iterator to the begining of the vector
 * @param xEnd: iterator to the end of the vector
 * @param pred: predicate that have to be satisfy
 * @return
 */
template<class T, class UnaryPredicate>
std::vector<LongInt> * who_satisfy_this(typename std::vector<T>::iterator xStart, typename std::vector<T>::iterator xEnd, UnaryPredicate pred)
{
    std::vector<sgt::LongInt> * result = new std::vector<sgt::LongInt>();
    for(;xStart != xEnd; ++xStart)
        result->push_back(pred(*xStart) ? 1 : 0);

    return result;
}

template<class T, class UnaryPredicate>
std::vector<LongInt> * who_satisfy_this(typename std::vector<T>::reverse_iterator xStart, typename std::vector<T>::reverse_iterator xEnd, UnaryPredicate pred)
{
    std::vector<sgt::LongInt> * result = new std::vector<sgt::LongInt>();
    for(;xStart != xEnd; ++xStart)
        result->push_back(pred(*xStart) ? 1 : 0);

    return result;
}

template<class T, class UnaryPredicate>
bool any_of(typename std::vector<T>::reverse_iterator xStart, typename std::vector<T>::reverse_iterator xEnd, UnaryPredicate pred)
{
    for(;xStart != xEnd; ++xStart)
    {
        if(pred(*xStart))
            return true;
    }
    return false;
}
//===========================----- filters -----=========================================//

template <class T>
std::vector<T> * lowPassF(std::vector<T> *S, float dt, float RC)
{
    std::vector<T> * y = new std::vector<T>(S->size());
    float alf = dt / (RC + dt);
    (*y)[0] = alf * (*S)[0];
    for(LongInt i = 1; i < S->size(); ++i)
    {
        (*y)[i] = alf * (*S)[i] + (1 - alf)* (*y)[i-1];
    }

    return y;
}

//===========================----- Operations on Vectors -----===========================//
template<class X, class Y>
std::vector<Y> * convertTo(std::vector<X> *v)
{
    std::vector<Y> *z = new std::vector<Y>(v->size());
    for(sgt::LongInt i = 0; i < (sgt::LongInt)z->size(); ++i)
        (*z)[i] = static_cast<Y>((*v)[i]);

    return z;
}

template<class T>
/**
 * @brief arange: Array of evenly spaced values.
 * @param init
 * @param end
 * @param step
 * @return
 */
std::vector<T> * arange(T init, T end, T step = 1)
{
    std::vector<T> *result = new std::vector<T>();

    while(ceil((end - init)/step) < end)
    {
        result->push_back(init);
        init += step;
    }

    return result;
}

template<class T>
std::vector<T> * append(std::vector<T> *S, std::vector<T> *Z)
{
    for(std::vector<T>::iterator zIt = Z->begin(); zIt != Z->end(); ++zIt)
        S->push_back(*(zIt));

    return S;
}

template<class T>
std::vector<T> * copy(typename std::vector<T>::iterator itStart, typename std::vector<T>::iterator itEnd, std::vector<T> *Z)
{
    Z->resize(itEnd - itStart);
    typename std::vector<T>::iterator itStartResult = Z->begin();
    for(;itStart != itEnd; ++itStart, ++itStartResult)
        *itStartResult = *itStart;

    return Z;
}

template<class T>
std::vector<T> * copy(typename std::vector<T>::reverse_iterator itStart, typename std::vector<T>::reverse_iterator itEnd, std::vector<T> *Z)
{
    if(!Z)
        Z = new std::vector<T>();

    for(T t ;itStart != itEnd; ++itStart, t =*itStart)
        Z->push_back(*itStart);

    return Z;
}

template<class X, class Z>
void copy(typename std::vector<X>::iterator itStart, typename std::vector<X>::iterator itEnd, typename std::vector<Z>::iterator itStartResult)
{
    for(;itStart != itEnd; ++itStart, ++itStartResult)
        *itStartResult = static_cast<Z>((*itStart));
}

template<class T>
std::vector<T> * copy_and_overwrite(typename std::vector<T>::iterator itStart, typename std::vector<T>::iterator itEnd, typename std::vector<T>::iterator startOverwrite, typename std::vector<T>::iterator endOverwrite, T x, std::vector<T> *&outVect)
{
    if(!outVect)
        outVect = new std::vector<T>();

    outVect->resize(itEnd - itStart);
    typename std::vector<T>::iterator itStartResult = outVect->begin();

    for(;itStart != itEnd; ++itStart, ++itStartResult)
        *itStartResult = itStart >= startOverwrite && itStart <= endOverwrite ? x : *itStart;

    return outVect;
}

template<class T>
void reverse(typename std::vector<T>::iterator itStart, typename std::vector<T>::iterator itEnd)
{
    while((itStart != itEnd) && (itStart != --itEnd))
            std::iter_swap(itStart++, itEnd);
}
//===========================----- Extrema -----=========================================//
enum ExtremaType
{
    UNSET_EXTREMA_TYPE  = 0,
    EXT_MAXIMUM_VALUE   = 1,
    EXT_MINIMUM_VALUE   = 2,
    EXT_MAXIMUM_TIME    = 3,
    EXT_MINIMUM_TIME    = 4,
    EXT_ZERO            = 5
};

template<class T>
class Extrema
{
public:
    const LongInt NON_ZERO = NAN;
    std::vector<T>      *local_max_values;
    std::vector<T>    *local_max_time;
    std::vector<T>      *local_min_values;
    std::vector<T>    *local_min_time;
    std::vector<LongInt>    *indzer;
    Extrema()
    {
        local_max_values = nullptr;
        local_max_time = nullptr;
        local_min_values = nullptr;
        local_min_time = nullptr;
        indzer = nullptr;
    }
    Extrema(std::vector<T> *lMaxV,
                 std::vector<T> *lMaxT,
                 std::vector<T> *lMinV,
                 std::vector<T> *lMinT,
                 std::vector<LongInt> *indZeros)
    {
        local_max_values = lMaxV;
        local_max_time = lMaxT;
        local_min_values = lMinV;
        local_min_time = lMinT;
        indzer = indZeros;
    }

    T invalidValue(){return NON_ZERO;}

    bool hasMaximum()
    {
        return !local_max_values || (local_max_values && local_max_values->size() == 0) ? false : true;
    }
    bool hasMinimum()
    {
        return !local_min_values || (local_min_values && local_min_values->size() == 0) ? false : true;
    }
    bool hasZeros()
    {
        return !indzer || (indzer && indzer->size() == 0) ? false : true;
    }

    T getMaxAt(long long int p)
    {
        T maxAt = invalidValue();
        if(p == -1 && maxCount() > 0)
            maxAt = (*local_max_values)[maxCount() - 1];
        else if(maxCount() > 0 && p < maxCount())
            maxAt = (*local_max_values)[p];
        return maxAt;
    }
    T getTimeOfMax(T max)
    {
        typename std::vector<T>::iterator itMaxPos = std::find(maxBegin(), maxEnd(), max);
        T timeOfMax = itMaxPos == maxEnd() ? invalidValue() : getTimeOfMaxAt(itMaxPos - maxBegin());
        return timeOfMax;
    }
    typename std::vector<T> getMaxInRange(T startTimePos, T endTimePos)
    {
        typename std::vector<T> maxVect;

        if(maxCount() == 0)
            return maxVect;

        int p = 0;
        T pivotTime = getTimeOfMaxAt(p);
        while ((p < maxCount() && pivotTime <= endTimePos) || (endTimePos == static_cast<T>(-1) && p < maxCount()))
        {
            if((endTimePos == static_cast<T>(-1) && pivotTime >= startTimePos) ||
               (pivotTime >= startTimePos && pivotTime <= endTimePos))
            {
                maxVect.push_back(getMaxAt(p));
            }
            ++p;
            pivotTime = getTimeOfMaxAt(p);
        }

        return maxVect;
    }
    T getTimeOfMaxAt(long long int p)
    {
        T timeAt = invalidValue();
        if(p == -1 && maxCount() > 0)
            timeAt = (*local_max_time)[maxCount() - 1];
        else if(maxCount() > 0 && p < maxCount())
            timeAt = (*local_max_time)[p];
        return timeAt;
    }
    void removeMax(T max)
    {
        T timeMax = getTimeOfMax(max);
        if(timeMax != invalidValue())
        {
            typename std::vector<T>::iterator itToTimeMax = std::find(local_max_time->begin(), local_max_time->end(), timeMax);
            if(itToTimeMax != local_max_time->end())
                local_max_time->erase(itToTimeMax);
        }
        typename std::vector<T>::iterator itToMax = std::find(maxBegin(), maxEnd(), max);
        if(itToMax != maxEnd())
            local_max_values->erase(itToMax);
    }

    int maxCount(){return !local_max_values || (local_max_values && local_max_values->size() == 0) ? 0 : static_cast<int>(local_max_values->size());}

    typename std::vector<T>::iterator maxBegin(){return local_max_values->begin();}
    typename std::vector<T>::iterator maxEnd(){return local_max_values->end();}
    typename std::vector<T>::reverse_iterator rMaxBegin(){return local_max_values->rbegin();}
    typename std::vector<T>::reverse_iterator rMaxEnd(){return local_max_values->rend();}

    T getMinAt(long long int p)
    {
        T minAt = invalidValue();
        if(p == -1 && minCount() > 0)
            minAt = (*local_min_values)[minCount() - 1];
        else if(minCount() > 0 && p < minCount())
            minAt = (*local_min_values)[p];
        return minAt;
    }
    T getTimeOfMinAt(long long int p)
    {
        T timeAt = invalidValue();
        if(p == -1 && minCount() > 0)
            timeAt = (*local_min_time)[minCount() - 1];
        else if(minCount() > 0 && p < minCount())
            timeAt = (*local_min_time)[p];
        return timeAt;
    }

    int minCount(){return !local_min_values || (local_min_values && local_min_values->size() == 0) ? 0 : static_cast<int>(local_min_values->size());}

    T getTimeOfMin(T min)
    {
        typename std::vector<T>::iterator itMinPos = std::find(minBegin(), minEnd(), min);
        T timeOfMin = itMinPos == minEnd() ? invalidValue() : getTimeOfMinAt(itMinPos - minBegin());
        return timeOfMin;
    }
    typename std::vector<T> getMinInRange(T startTimePos, T endTimePos)
    {
        typename std::vector<T> minVect;

        if(minCount() == 0)
            return minVect;

        int p = 0;
        T pivotTime = getTimeOfMinAt(p);
        while ((p < minCount() && pivotTime <= endTimePos)  || (endTimePos == static_cast<T>(-1) && p < maxCount()))
        {
            //If endTimePos == -1 looking for between the start time position and the end
            if((endTimePos == static_cast<T>(-1) && pivotTime >= startTimePos) ||
                    (pivotTime >= startTimePos && pivotTime <= endTimePos))
            {
                minVect.push_back(getMinAt(p));
            }

            ++p;
            pivotTime = getTimeOfMinAt(p);
        }

        return minVect;
    }
    void removeMin(T min)
    {
        T timeMin = getTimeOfMin(min);
        if(timeMin != invalidValue())
        {
            typename std::vector<T>::iterator itToTimeMin = std::find(local_min_time->begin(), local_min_time->end(), timeMin);
            if(itToTimeMin != local_min_time->end())
                local_min_time->erase(itToTimeMin);
        }
        typename std::vector<T>::iterator itToMin = std::find(minBegin(), minEnd(), min);
        if(itToMin != minEnd())
            local_min_values->erase(itToMin);
    }

    T getZeroAt(long long int p)
    {
        T zeroAt = invalidValue();
        if(p == -1 && zeroCount() > 0)
            zeroAt = (*indzer)[zeroCount() - 1];
        else if(zeroCount() > 0 && p < zeroCount())
            zeroAt = (*indzer)[p];
        return zeroAt;
    }
    int zeroCount(){return !indzer || (indzer && indzer->size() == 0) ? 0 : static_cast<int>(indzer->size());}
    void removeExtremeWithoutZeroCrossing()
    {
        if(zeroCount() == 0)
            return;

        T endTime;
        T startTime = static_cast<T>(0);

        typename std::vector<T> maxToRemove;
        typename std::vector<T> minToRemove;

        for(int i0 = 0; i0 <= zeroCount(); ++i0)
        {
            endTime = i0 == zeroCount() ? static_cast<T>(-1) :getZeroAt(i0);
            //look for the biggest of maximums
            if(maxCount() > 0)
            {
                typename std::vector<T> maxVect = getMaxInRange(startTime, endTime);
                if(maxVect.size() > 1)
                {
                    typename std::vector<T>::iterator largestMax = std::max_element(maxVect.begin(),maxVect.end());
                    for (int i1 = 0; i1 < maxVect.size(); ++i1)
                    {
                        if(maxVect[i1] != *largestMax)
                            maxToRemove.push_back(maxVect[i1]);
                    }
                }
            }
            //look for the smallest of minimums
            if(minCount() > 0)
            {
                typename std::vector<T> minVect = getMinInRange(startTime, endTime);
                if(minVect.size() > 1)
                {
                    typename std::vector<T>::iterator smallestMin = std::min_element(minVect.begin(),minVect.end());
                    for (int i2 = 0; i2 < minVect.size(); ++i2)
                    {
                        if(minVect[i2] != *smallestMin)
                            minToRemove.push_back(minVect[i2]);
                    }
                }
            }
            startTime = endTime;
        }
        //remove the redundant maximums
        for (int i3 = 0; i3 < maxToRemove.size(); ++i3)
            removeMax(maxToRemove[i3]);

        //remove the redundant minimums
        for (int i4 = 0; i4 < minToRemove.size(); ++i4)
            removeMin(minToRemove[i4]);
    }
    T rangeBetweenExtrema(long long p)
    {
        T result = invalidValue();
        if(p == -1)
        {
            if(maxCount() > 0 && maxCount() == minCount())
                result = std::abs((*local_max_values)[maxCount() - 1] - (*local_min_values)[maxCount() - 1]);
            else if(maxCount() > minCount())
            {
                result = std::abs((*local_max_values)[maxCount() -1]);
            }
            else if(maxCount() < minCount())
            {
                result = std::abs((*local_min_values)[minCount() -1]);
            }
        }
        else if((maxCount() > 0 && p < maxCount()) || (minCount() > 0 && p < minCount()))
        {
            if(maxCount() == minCount() || (maxCount() > minCount() && p < minCount()) || (maxCount() < minCount() && p < maxCount()))
            {
                result = std::abs((*local_max_values)[p] - (*local_min_values)[p]);
            }
            else if(maxCount() > minCount() && p >= minCount())
            {
                result = std::abs((*local_max_values)[p]);
            }
            else if(maxCount() < minCount() && p >= maxCount())
            {
                result = std::abs((*local_min_values)[p]);
            }
        }
        return result;
    }
    typename std::vector<T>::iterator minBegin(){return local_min_values->begin();}
    typename std::vector<T>::iterator minEnd(){return local_min_values->end();}
    typename std::vector<T>::reverse_iterator rMinBegin(){return local_min_values->rbegin();}
    typename std::vector<T>::reverse_iterator rMinEnd(){return local_min_values->rend();}

    T operator()(long long int i, ExtremaType extT)
    {
        switch (extT) {
        case EXT_MAXIMUM_VALUE:
            return getMaxAt(i);
            break;
        case EXT_MINIMUM_VALUE:
            return getMinAt(i);
            break;
        case EXT_MINIMUM_TIME:
            return getTimeOfMinAt(i);
            break;
        case EXT_MAXIMUM_TIME:
            return getTimeOfMaxAt(i);
            break;
        case EXT_ZERO:
            return static_cast<T>(getZeroAt(i));
            break;
        default:
            break;
        }
        return NAN;
    }

    const T& operator()(long long int i, ExtremaType extT) const
    {
        switch (extT) {
        case EXT_MAXIMUM_VALUE:
            return getMaxAt(i);
            break;
        case EXT_MINIMUM_VALUE:
            return getMinAt(i);
            break;
        case EXT_MINIMUM_TIME:
            return getTimeOfMinAt(i);
            break;
        case EXT_MAXIMUM_TIME:
            return getTimeOfMaxAt(i);
            break;
        default:
            break;
        }
    }

    ~Extrema()
    {
        releaseMemory<T>(local_max_values);
        releaseMemory<T>(local_max_time);
        releaseMemory<T>(local_min_values);
        releaseMemory<T>(local_min_time);
        releaseMemory<LongInt>(indzer);
    }
};

template<class T>
class MaxMinSpline
{
public:
    std::vector<T> *max_spline;
    std::vector<T> *min_spline;
    Extrema<T>            *ext;

    MaxMinSpline()
    {
        max_spline  = nullptr;
        min_spline  = nullptr;
        ext         = nullptr;
    }

    MaxMinSpline(std::vector<T> *maxSpline, std::vector<T> *minSpline, Extrema<T> *_ext)
    {
        max_spline  = maxSpline;
        min_spline  = minSpline;
        ext = _ext;
    }

    ~MaxMinSpline()
    {
        releaseMemory(max_spline);
        releaseMemory(min_spline);
        if(ext)
            delete ext;
    }
};

/**
 * @brief find_extrema_simple: Performs extrema detection, where extremum is defined as a point, that is above/below its neighbours.
 * @param S: signal's values
 */
template<class T>
sgt::Extrema<T> * find_extrema_simple(typename std::vector<T>::iterator itStart, typename std::vector<T>::iterator itEnd, bool absValues = false)
{
    int size = itEnd - itStart;
    std::vector<sgt::LongInt>   *d = nullptr;
    std::vector<sgt::LongInt> *indzer = sgt::nonzeroIndexes<T>(itStart, itEnd);

    //if any S[i] == 0
    if(std::any_of(itStart, itEnd, [](T x){/*qDebug() << x;*/ return x == 0.0;}))
    {
        std::vector<sgt::LongInt> *iz =sgt::indexes_of<T>(itStart, itEnd, [](T x){return x == 0.0;});

        d = sgt::diff<sgt::LongInt>(iz, absValues);

        std::vector<sgt::LongInt> *indz = nullptr;

        if(std::any_of(d->begin(), d->end(), [](T x){return x == 1;}))
        {
//            std::vector<sgt::LongInt> *zer =sgt::indexes_of<float>(S->begin(), S->end(), [](float x){return x == 0.0;});
            std::vector<sgt::LongInt> *dz = sgt::who_satisfy_this<T>(itStart, itEnd, [](T x){return x == 0.0;});
            dz->push_back(0);
            dz->insert(dz->begin(), 0);

            std::vector<sgt::LongInt> *debz = sgt::indexes_of<sgt::LongInt>(dz->begin(), dz->end(), [](sgt::LongInt x){return x == 1;});
            std::vector<sgt::LongInt> *finz = sgt::indexes_of<sgt::LongInt>(dz->begin(), dz->end(), [](sgt::LongInt x){return x == -1;}, -1);

            if(!indz)
                indz = new std::vector<sgt::LongInt>();

            //TODO: check this for file203 rpeakIndex = 0
            if(debz && finz)
                for(sgt::LongInt i = 0; i < (sgt::LongInt)debz->size(); ++i)
                    indz->push_back(std::round(((*debz)[i] + (*finz)[i])/2.0));

            sgt::releaseMemory<sgt::LongInt>(debz);
            sgt::releaseMemory<sgt::LongInt>(finz);
            sgt::releaseMemory<sgt::LongInt>(dz);
            sgt::releaseMemory<sgt::LongInt>(iz);
        }
        else
        {
            sgt::releaseMemory<sgt::LongInt>(indz);
            indz = iz;
        }

        //TODO: verify this....
        if(indzer && indz && indz->size() > 0)
            sgt::append<sgt::LongInt>(indzer, indz);
        if(indzer && indzer->size() > 0)
            std::sort(indzer->begin(), indzer->end());

        sgt::releaseMemory<sgt::LongInt>(indz);
        sgt::releaseMemory<sgt::LongInt>(d);
    }

    //Finds local extrema
    std::vector<T> *sDiff = sgt::diff<T>(itStart, itEnd, absValues);
    std::vector<sgt::LongInt> *indmin = sgt::indexes_of_min<T>(itStart, itEnd, absValues);
    std::vector<sgt::LongInt> *indmax = sgt::indexes_of_max<T>(itStart, itEnd, absValues);

    //When two or more points have the same value
    if(std::any_of(sDiff->begin(), sDiff->end(),[](T x){return x == 0.0;}))
    {
        std::vector<sgt::LongInt> imax;
        std::vector<sgt::LongInt> imin;

        std::vector<sgt::LongInt> *bad = sgt::who_satisfy_this<T>(sDiff->begin(), sDiff->end(), [](T x){return x == 0.0;});
        bad->push_back(0);
        bad->insert(bad->begin(),0);
        std::vector<sgt::LongInt> *dd = sgt::diff<sgt::LongInt>(bad, absValues);

        sgt::releaseMemory<sgt::LongInt>(bad);

        std::vector<sgt::LongInt> *debs =sgt::indexes_of<sgt::LongInt>(dd->begin(), dd->end(), [](sgt::LongInt x){return x == 1;});
        std::vector<sgt::LongInt> *fins =sgt::indexes_of<sgt::LongInt>(dd->begin(), dd->end(), [](sgt::LongInt x){return x == -1;});

        sgt::releaseMemory<sgt::LongInt>(dd);

        if((*debs)[0] == 1)
        {
            if(debs->size() >1)
            {
                //debs, fins = debs[1:], fins[1:];
                debs->erase(debs->begin());
                fins->erase(fins->begin());
            }
            else
            {
                debs->clear();
                fins->clear();
            }
        }

        if(debs->size() > 0)
        {
            if((*fins)[fins->size()-1] == size-1)
                if(debs->size() > 1)
                {
                    //debs, fins = debs[:-1], fins[:-1]
                    debs->pop_back();
                    fins->pop_back();
                }
                else
                {
                    debs->clear();
                    fins->clear();
                }
        }

        sgt::LongInt lc = debs->size();

        if(lc > 0)
        {
            for(sgt::LongInt k = 0; k < lc; ++k)
            {
                if((*debs)[k]-1 >=0 && (*sDiff)[(*debs)[k]-1] > 0)
                {
                    if((*sDiff)[(*fins)[k]] < 0)
                        imax.push_back(floor(((*fins)[k]+(*debs)[k])/2.0)); //the original code uses round

                }
                else
                {
                    if((*sDiff)[(*fins)[k]] > 0)
                        imin.push_back(floor(((*fins)[k]+ (*debs)[k])/2.0)); //the original code uses round
                }
            }
        }

        sgt::releaseMemory<sgt::LongInt>(debs);
        sgt::releaseMemory<sgt::LongInt>(fins);

        if(indmax && imax.size() > 0)
        {
            //indmax = indmax.tolist()
            foreach (sgt::LongInt x, imax)
            {
                indmax->push_back(x);
            }
            std::sort(indmax->begin(), indmax->end());
        }

        if(indmin && imin.size() > 0)
        {
            //indmin = indmin.tolist()
            foreach (sgt::LongInt x, imin)
            {
                indmin->push_back(x);
            }
            std::sort(indmin->begin(), indmin->end());
        }
    }

    std::vector<T> *lMaxValues = new std::vector<T>();
    std::vector<T> *lMinValues = new std::vector<T>();
    std::vector<T> *lMaxTimes = new std::vector<T>();
    std::vector<T> *lMinTimes = new std::vector<T>();

    for(sgt::LongInt i = 0; indmax && i < (sgt::LongInt)indmax->size(); ++i)
    {
        //lMaxValues->push_back((*S)[(*indmax)[i]]);
        lMaxValues->push_back(*(itStart+(*indmax)[i]));
        //lMaxTimes->push_back((*signalTime)[(*indmax)[i]]);
        lMaxTimes->push_back((*indmax)[i]);
    }
    for(sgt::LongInt i = 0; indmin && i < (sgt::LongInt)indmin->size(); ++i)
    {
        //lMinValues->push_back((*S)[(*indmin)[i]]);
        lMinValues->push_back(*(itStart+(*indmin)[i]));
        //lMinTimes->push_back((*signalTime)[(*indmin)[i]]);
        lMinTimes->push_back((*indmin)[i]);
    }

    sgt::releaseMemory<T>(sDiff);
    sgt::releaseMemory<sgt::LongInt>(indmax);
    sgt::releaseMemory<sgt::LongInt>(indmin);

    return new sgt::Extrema<T>(lMaxValues, lMaxTimes, lMinValues, lMinTimes, indzer);
}

/**
 * @brief find_extrema_simple: Performs extrema detection, where extremum is defined as a point, that is above/below its neighbours.
 * @param S: signal's values
 */
template<class T>
sgt::Extrema<T> * find_extrema_simple(typename std::vector<T>::reverse_iterator itStart, typename std::vector<T>::reverse_iterator itEnd, bool absValues = false)
{
    int size = itEnd-itStart;
    std::vector<sgt::LongInt>   *d = nullptr;
    std::vector<sgt::LongInt> *indzer = sgt::nonzeroIndexes<T>(itStart, itEnd);

    //if any S[i] == 0
    if(sgt::any_of<T>(itStart, itEnd, [](T x){return x == 0.0;}))
    {
        std::vector<sgt::LongInt> *iz =sgt::indexes_of<T>(itStart, itEnd, [](T x){return x == 0.0;});

        d = sgt::diff<sgt::LongInt>(iz, absValues);

        std::vector<sgt::LongInt> *indz = nullptr;

        if(std::any_of(d->begin(), d->end(), [](T x){return x == 1;}))
        {
            std::vector<sgt::LongInt> *dz = sgt::who_satisfy_this<T>(itStart, itEnd, [](T x){return x == 0.0;});
            dz->push_back(0);
            dz->insert(dz->begin(), 0);

            std::vector<sgt::LongInt> *debz =sgt::indexes_of<sgt::LongInt>(dz->begin(), dz->end(), [](sgt::LongInt x){return x == 1;});
            std::vector<sgt::LongInt> *finz =sgt::indexes_of<sgt::LongInt>(dz->begin(), dz->end(), [](sgt::LongInt x){return x == -1;}, -1);

            if(!indz)
                indz = new std::vector<sgt::LongInt>();

            //TODO: check this for file104 rpeakIndex = 0
            if(debz && finz)
                for(sgt::LongInt i = 0; i < (sgt::LongInt)debz->size(); ++i)
                    indz->push_back(std::round(((*debz)[i] + (*finz)[i])/2.0));

            sgt::releaseMemory<sgt::LongInt>(debz);
            sgt::releaseMemory<sgt::LongInt>(finz);
            sgt::releaseMemory<sgt::LongInt>(dz);
            sgt::releaseMemory<sgt::LongInt>(iz);
        }
        else
        {
            sgt::releaseMemory<sgt::LongInt>(indz);
            indz = iz;
        }

        //TODO: verify this....
        if(indz && indz->size() > 0)
            sgt::append<sgt::LongInt>(indzer, indz);
        if(indzer && indzer->size() > 0)
            std::sort(indzer->begin(), indzer->end());

        sgt::releaseMemory<sgt::LongInt>(indz);
        sgt::releaseMemory<sgt::LongInt>(d);
    }

    //Finds local extrema
    std::vector<T> *sDiff = sgt::diff<T>(itStart, itEnd, absValues);
    if(sDiff->size() == size +1)
        sDiff->pop_back();

    std::vector<sgt::LongInt> *indmin = sgt::indexes_of_min<T>(itStart, itEnd, absValues);
    std::vector<sgt::LongInt> *indmax = sgt::indexes_of_max<T>(itStart, itEnd, absValues);

    //When two or more points have the same value
    if(std::any_of(sDiff->begin(), sDiff->end(),[](T x){return x == 0.0;}))
    {
        std::vector<sgt::LongInt> imax;
        std::vector<sgt::LongInt> imin;

        std::vector<sgt::LongInt> *bad = sgt::who_satisfy_this<T>(sDiff->begin(), sDiff->end(), [](T x){return x == 0.0;});
        bad->push_back(0);
        bad->insert(bad->begin(),0);
        std::vector<sgt::LongInt> *dd = sgt::diff<sgt::LongInt>(bad, absValues);

        sgt::releaseMemory<sgt::LongInt>(bad);

        std::vector<sgt::LongInt> *debs =sgt::indexes_of<sgt::LongInt>(dd->begin(), dd->end(), [](sgt::LongInt x){return x == 1;});
        std::vector<sgt::LongInt> *fins =sgt::indexes_of<sgt::LongInt>(dd->begin(), dd->end(), [](sgt::LongInt x){return x == -1;});

        sgt::releaseMemory<sgt::LongInt>(dd);

        if((*debs)[0] == 1)
        {
            if(debs->size() >1)
            {
                //debs, fins = debs[1:], fins[1:];
                debs->erase(debs->begin());
                fins->erase(fins->begin());
            }
            else
            {
                debs->clear();
                fins->clear();
            }
        }

        if(debs->size() > 0)
        {
            if((*fins)[fins->size()-1] == size-1 || (*fins)[fins->size()-1] == size)
                if(debs->size() > 1)
                {
                    //debs, fins = debs[:-1], fins[:-1]
                    debs->pop_back();
                    fins->pop_back();
                }
                else
                {
                    debs->clear();
                    fins->clear();
                }
        }

        sgt::LongInt lc = debs->size();

        if(lc > 0)
        {
            for(sgt::LongInt k = 0; k < lc; ++k)
            {
                if((*debs)[k]-1 >=0 && (*sDiff)[(*debs)[k]-1] > 0)
                {
                    if((*sDiff)[(*fins)[k]] < 0)
                        imax.push_back(floor(((*fins)[k]+(*debs)[k])/2.0)); //the original code uses round

                }
                else
                {
                    if((*sDiff)[(*fins)[k]] > 0)
                        imin.push_back(floor(((*fins)[k]+ (*debs)[k])/2.0)); //the original code uses round
                }
            }
        }

        sgt::releaseMemory<sgt::LongInt>(debs);
        sgt::releaseMemory<sgt::LongInt>(fins);

        if(indmax && imax.size() > 0)
        {
            //indmax = indmax.tolist()
            foreach (sgt::LongInt x, imax)
            {
                indmax->push_back(x);
            }
            std::sort(indmax->begin(), indmax->end());
        }

        if(indmin && imin.size() > 0)
        {
            //indmin = indmin.tolist()
            foreach (sgt::LongInt x, imin)
            {
                indmin->push_back(x);
            }
            std::sort(indmin->begin(), indmin->end());
        }
    }

    std::vector<T> *lMaxValues = new std::vector<T>();
    std::vector<T> *lMinValues = new std::vector<T>();
    std::vector<T> *lMaxTimes = new std::vector<T>();
    std::vector<T> *lMinTimes = new std::vector<T>();

    for(sgt::LongInt i = 0; indmax && i < (sgt::LongInt)indmax->size(); ++i)
    {
        //lMaxValues->push_back((*S)[(*indmax)[i]]);
        lMaxValues->push_back(*(itStart+(*indmax)[i]));
        //lMaxTimes->push_back((*signalTime)[(*indmax)[i]]);
        lMaxTimes->push_back((*indmax)[i]);
    }
    for(sgt::LongInt i = 0; indmin && i < (sgt::LongInt)indmin->size(); ++i)
    {
        //lMinValues->push_back((*S)[(*indmin)[i]]);
        lMinValues->push_back(*(itStart+(*indmin)[i]));
        //lMinTimes->push_back((*signalTime)[(*indmin)[i]]);
        lMinTimes->push_back((*indmin)[i]);
    }

    sgt::releaseMemory<T>(sDiff);
    sgt::releaseMemory<sgt::LongInt>(indmax);
    sgt::releaseMemory<sgt::LongInt>(indmin);

    return new sgt::Extrema<T>(lMaxValues, lMaxTimes, lMinValues, lMinTimes, indzer);
}

template<class T>
sgt::Extrema<T> *find_extrema_parabol(typename std::vector<T>::iterator itStart, typename std::vector<T>::iterator itEnd, bool absValues = false)
{
    std::vector<sgt::LongInt>   *d = nullptr;
    std::vector<sgt::LongInt> *indzer = sgt::nonzeroIndexes<T>(itStart, itEnd);

    //if any S[i] == 0
    if(std::any_of(itStart, itEnd, [](T x){return x == 0.0;}))
    {
        std::vector<sgt::LongInt> *iz =sgt::indexes_of<T>(itStart, itEnd, [](T x){return x == 0.0;});

        d = sgt::diff<sgt::LongInt>(iz);

        std::vector<sgt::LongInt> *indz = nullptr;

        if(std::any_of(d->begin(), d->end(), [](T x){return x == 1;}))
        {
            std::vector<sgt::LongInt> *dz = sgt::who_satisfy_this<T>(itStart, itEnd, [](T x){return x == 0.0;});
            dz->push_back(0);
            dz->insert(dz->begin(), 0);

            std::vector<sgt::LongInt> *debz =sgt::indexes_of<sgt::LongInt>(dz->begin(), dz->end(), [](sgt::LongInt x){return x == 1;});
            std::vector<sgt::LongInt> *finz =sgt::indexes_of<sgt::LongInt>(dz->begin(), dz->end(), [](sgt::LongInt x){return x == -1;}, -1);

            if(!indz)
                indz = new std::vector<sgt::LongInt>();

            for(sgt::LongInt i = 0; i < (sgt::LongInt)debz->size(); ++i)
                indz->push_back(std::round(((*debz)[i] + (*finz)[i])/2.0));

            sgt::releaseMemory<sgt::LongInt>(debz);
            sgt::releaseMemory<sgt::LongInt>(finz);
            sgt::releaseMemory<sgt::LongInt>(dz);
        }
        else
        {
            sgt::releaseMemory<sgt::LongInt>(indz);
            indz = iz;
        }

        sgt::append<sgt::LongInt>(indzer, indz);
        std::sort(indzer->begin(), indzer->end());

        sgt::releaseMemory<sgt::LongInt>(indz);
        sgt::releaseMemory<sgt::LongInt>(d);
    }


    //dt = float(T[1]-T[0])
    T dt = 1;
    T scale = 2.0*dt*dt;
//    float scale = 2.0f*dt*dt;
    //idx = self._not_duplicate(S)
    std::vector<sgt::LongInt> *idx = not_duplicate<T>(itStart, itEnd);

    //T = T[idx]
    //S = S[idx]
    std::vector<T> *cT = new std::vector<T>();
    std::vector<T> *cS = new std::vector<T>();
    for(sgt::LongInt i = 0; i < (sgt::LongInt)idx->size(); ++i)
    {
        cT->push_back((*idx)[i]);
        cS->push_back(*(itStart + (*idx)[i]));
    }

    //p - previous
    //0 - current
    //n - next
    std::vector<T> Tp, T0, Tn, Sp, S0, Sn;

    //Tp, T0, Tn = T[:-2], T[1:-1], T[2:]
    sgt::copy(cT->begin(), cT->end() - 2, &Tp);
    sgt::copy(cT->begin()+1, cT->end() - 1, &T0);
    sgt::copy(cT->begin()+2, cT->end(), &Tn);

    //Sp, S0, Sn = S[:-2], S[1:-1], S[2:]
    sgt::copy(cS->begin(), cS->end() - 2, &Sp);
    sgt::copy(cS->begin()+1, cS->end() - 1, &S0);
    sgt::copy(cS->begin()+2, cS->end(), &Sn);

    //a = Sn + Sp - 2*S0
    //b = 2*(Tn+Tp)*S0 - ((Tn+T0)*Sp+(T0+Tp)*Sn)
    //c = Sp*T0*Tn -2*Tp*S0*Tn + Tp*T0*Sn
    //TnTp, T0Tn, TpT0 = Tn-Tp, T0-Tn, Tp-T0
    //scale = Tp*Tn*Tn + Tp*Tp*T0 + T0*T0*Tn - Tp*Tp*Tn - Tp*T0*T0 - T0*Tn*Tn
    //a = a/scale
    //b = b/scale
    //c = c/scale
    //a[a==0] = 1e-14
    //tVertex = -0.5*b/a
    T a = 0.0f, b = 0.0f, c = 0.0f, tVertex = 0.0f;
    T TnTp = 0.0f, T0Tn = 0.0f, TpT0 = 0.0f;
    idx->clear();

    std::vector<T> *local_max_pos = new std::vector<T>();
    std::vector<T> *local_max_val = new std::vector<T>();
    std::vector<T> *local_min_pos = new std::vector<T>();
    std::vector<T> *local_min_val = new std::vector<T>();

    for(sgt::LongInt i = 0; i < (sgt::LongInt)Tp.size(); ++i)
    {

        scale = Tp[i]*Tn[i]*Tn[i] + Tp[i]*Tp[i]*T0[i] + T0[i]*T0[i]*Tn[i] - Tp[i]*Tp[i]*Tn[i] - Tp[i]*T0[i]*T0[i] - T0[i]*Tn[i]*Tn[i];
        TnTp = Tn[i]-Tp[i];
        T0Tn = T0[i]-Tn[i];
        TpT0 = Tp[i]-T0[i];

        a = (T0Tn*Sp[i] + TnTp*S0[i] + TpT0*Sn[i])/scale;
        a == a == 0.0f ? 1e-14 : a;
        b = ((S0[i]-Sn[i])*Tp[i]*Tp[i] + (Sn[i]-Sp[i])*T0[i]*T0[i] + (Sp[i]-S0[i])*Tn[i]*Tn[i])/scale;
        c = (T0[i]*Tn[i]*T0Tn*Sp[i] + Tn[i]*Tp[i]*TnTp*S0[i] + Tp[i]*T0[i]*TpT0*Sn[i])/scale;

        tVertex = -0.5*b/a;
        if((tVertex<T0[i]+0.5*(Tn[i]-T0[i])) && (tVertex>=T0[i]-0.5*(T0[i]-Tp[i])))
        {
            if(a < 0)
            {
                local_max_pos->push_back(tVertex);
                local_max_val->push_back(a*tVertex*tVertex + b*tVertex + c);
            }
            else
            {
                local_min_pos->push_back(tVertex);
                local_min_val->push_back(a*tVertex*tVertex +b*tVertex + c);
            }
        }

    }

    return new sgt::Extrema<T>(local_max_val, local_max_pos, local_min_val, local_min_pos, indzer);
}

template <class T>
class Utils
{
private:
     Utils();
public:
     /**
      * Swap two positions.
      */
     static void swap(std::vector<T> &arr, int i, int j)
     {
         double a = arr[i];
         arr[i] = arr[j];
         arr[j] = a;
     }

    /**
     * To restore the max-heap condition when a node's priority is increased. We move up the heap, exchaning the node at position k with its parent
     * (at postion k/2) if necessary, continuing as long as a[k/2] &lt; a[k] or until we reach the top of the heap.
     */
    static void shiftUp(std::vector<T> &arr, int k)
    {
        while(k > 1 && arr[k/2] < arr[k])
        {
            swap(arr, k, k/2);
            k = k/2;
        }
    }
    /**
     * To restore the max-heap condition when a node's priority is decreased. We move down the heap, exchanging the node at position k with the larger
     * of that node's two children if necessary and stopping when the node at k is not smaller than either child or the bottom is reached. Note that
     * if n is even and k is n/2, then the node at k has only one child -- this case must be treated properly.
     */
    static void shiftDown(std::vector<T> &arr, long long k, long long n)
    {
        while(2*k <= n)
        {
            long long j = 2 * k;
            if(j < n && arr[j] < arr[j + 1])
            {
                j++;
            }
            if(arr[k] >= arr[j])
            {
                break;
            }
            swap(arr, k, j);
            k = j;
        }
    }
    /**
     * The Euclidean distance.
     */
    static double euclideanDistance(std::vector<T> *x, std::vector<T>  *y)
    {
        if(x->size() != y->size())
            return std::numeric_limits<T>::max();


        double sum = 0.0;
        for(int i = 0; i < x->size(); i++)
        {
            double componentSum = (*x)[i] - (*y)[i];
            sum += componentSum * componentSum;
        }

        return std::sqrt(sum);
    }
    static double squaredEuclideanDistance(std::vector<T> *x, std::vector<T>  *y)
    {
        if(x->size() != y->size())
            return 0.0; //TOV: coul this be a valid output //throw new IllegalArgumentException("Input vector sizes are different.");


        double sum = 0.0;
        for(int i = 0; i < x->size(); i++)
        {
            double componentSum = (*x)[i] - (*y)[i];
            sum += componentSum * componentSum;
        }

        return sum;
    }

    static double centroid(std::vector<T> *x)
    {
        double _centroid = 0.0;
        size_t N = x->size();

        for(size_t i = 0; i < N; ++i)
            _centroid += (*x)[i];

        _centroid /= N;

        return _centroid;
    }
    static double radius(std::vector<T> *x)
    {
        double _radius = 0.0;
        size_t N = x->size();

        for(size_t i = 1; i < N; ++i)
        {
            double tmp = (*x)[i] - (*x)[0];
            _radius += tmp * tmp;
        }
        _radius /= N;
        _radius = std::sqrt(_radius);

        return _radius;
    }
    static double diameter(std::vector<T> *x)
    {
        double _diameter = 0.0;
        size_t N = x->size();

        for(size_t i = 1; i < N; ++i)
            for(size_t j = 1; j < N; ++j)
            {
                 double tmp = (*x)[i] - (*x)[j];
                 _diameter += tmp * tmp;
            }
        _diameter /= N * (N-1);
        _diameter = std::sqrt(_diameter);

        return _diameter;
    }
};

}

#endif // SGT_H
