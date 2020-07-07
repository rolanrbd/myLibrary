#include "Clustering.h"
#include <QQueue>
#include <QtAlgorithms>

//=========--------Utils--------========//
template <class T>
void clt::Utils<T>::swap(std::vector<T> &arr, int i, int j)
{
    double a = arr[i];
    arr[i] = arr[j];
    arr[j] = a;
}

template <class T>
void clt::Utils<T>::siftUp(std::vector<T> &arr, int k)
{
    while (k > 1 && arr[k/2] < arr[k])
    {
        swap(arr, k, k/2);
        k = k/2;
    }
}

template <class T>
void clt::Utils<T>::siftDown(std::vector<T> &arr, int k, int n)
{
    while (2*k <= n)
    {
        int j = 2 * k;
        if (j < n && arr[j] < arr[j + 1])
        {
            j++;
        }
        if (arr[k] >= arr[j])
        {
            break;
        }
        swap(arr, k, j);
        k = j;
    }
}

template <class T>
double clt::Utils<T>::euclideanDistance(std::vector<T> *x, std::vector<T> *y)
{
    if (x->size() != y->size())
        return 0.0;


    double sum = 0.0;
    for (int i = 0; i < x->size(); i++)
    {
        double componentSum = (*x)[i] - (*y)[i];
        sum += componentSum * componentSum;
    }

    return std::sqrt(sum);
}

template<class T>
double clt::Utils<T>::squaredEuclideanDistance(std::vector<T> *x, std::vector<T> *y)
{
    if (x->size() != y->size())
        return 0.0;


    double sum = 0.0;
    for (int i = 0; i < x->size(); i++)
    {
        double componentSum = (*x)[i] - (*y)[i];
        sum += componentSum * componentSum;
    }

    return sum;
}

template<class T>
double clt::Utils<T>::centroid(std::vector<T> *x)
{
    double _centroid = 0.0;
    size_t N = x->size();

    for(size_t i = 0; i < N; ++i)
        _centroid += (*x)[i];

    _centroid /= N;

    return _centroid;
}

template<class T>
double clt::Utils<T>::radius(std::vector<T> *x)
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

template<class T>
double clt::Utils<T>::diameter(std::vector<T> *x)
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

//=========--------IntHeapSelect--------========//

clt::IntHeapSelect::IntHeapSelect(int k)
{
    heap.resize(k);
    *this = clt::IntHeapSelect::IntHeapSelect(heap);
}

clt::IntHeapSelect::IntHeapSelect(std::vector<int> &hp)
{
    heap = hp;

    k = (int) hp.size();
    n = 0;
    sorted = false;
}

clt::IntHeapSelect::~IntHeapSelect()
{
    heap.clear();
}

void clt::IntHeapSelect::add(int datum)
{
    sorted = false;
    if (n < k)
    {
        heap[n++] = datum;
        if (n == k)
        {
            heapify(heap);
        }
    }
    else
    {
        n++;
        if (datum < heap[0])
        {
            heap[0] = datum;
            clt::Utils<int>::siftDown(heap, 0, k-1);
        }
    }
}

int clt::IntHeapSelect::peek()
{
    return heap[0];
}

int clt::IntHeapSelect::get(int i)
{
    if (i > std::min(k, n) - 1)
        return -1;


    if (i == k-1)
        return heap[0];


    if (!sorted)
    {
        sort(heap, std::min(k,n));
        sorted = true;
    }

    return heap[k-1-i];
}

void clt::IntHeapSelect::sort()
{
    if (!sorted)
    {
        sort(heap, std::min(k,n));
        sorted = true;
    }
}

void clt::IntHeapSelect::heapify(std::vector<int> &arr)
{
    int n = (int) arr.size();
    for (int i = n / 2 - 1; i >= 0; i--)
        clt::Utils<int>::siftDown(arr, i, n - 1);
}

void clt::IntHeapSelect::sort(std::vector<int> a, int n)
{
    int inc = 1;
    do
    {
        inc *= 3;
        inc++;
    } while (inc <= n);

    do
    {
        inc /= 3;
        for (int i = inc; i < n; i++)
        {
            int v = a[i];
            int j = i;
            while (a[j - inc] < v)
            {
                a[j] = a[j - inc];
                j -= inc;
                if (j < inc)
                {
                    break;
                }
            }
            a[j] = v;
        }
    } while (inc > 1);
}

//=========--------Linkage--------========//

void clt::Linkage::init(std::vector<std::vector<double> > * _proximity, size_t size)
{
    Q_UNUSED(_proximity)

    linkageSize = size; //(int)_proximity->size();
    clt::longInt proximitySize = linkageSize * ( (linkageSize+1) % 2 == 0 ? (linkageSize+1) / 2 : (linkageSize+1) / 2 + 1);
    proximity.resize(proximitySize);

    // row wise
    /*
    for (int i = 0, k = 0; i < size; i++) {
        double[] pi = proximity[i];
        for (int j = 0; j <= i; j++, k++) {
            this.proximity[k] = (float) pi[j];
        }
    }
    */

    // column wise
    /*
    for (int j = 0, k = 0; j < linkageSize; j++)
    {
        for (int i = j; i < linkageSize; i++, k++)
        {
            proximity[k] = (*_proximity)[i][j];
        }
    }*/
}

clt::longInt clt::Linkage::index(clt::longInt i, clt::longInt j)
{
    // row wise
    // return i > j ? i*(i+1)/2 + j : j*(j+1)/2 + i;
    // column wise
    return i > j ?
                proximity.size() - (linkageSize - j) * (linkageSize - j + 1)/2 + i - j :
                proximity.size() - (linkageSize - i) * (linkageSize-i+1)/2 + j - i;
}

clt::Linkage::~Linkage()
{
    proximity.clear();
}

clt::longInt clt::Linkage::size()
{
    return linkageSize;
}

float clt::Linkage::d(clt::longInt i, clt::longInt j)
{
    clt::longInt _index = index(i, j);
    return proximity[_index];
}

//=========--------WardLinkage--------========//

clt::WardLinkage::~WardLinkage()
{
    if(n)
    {
        n->clear();
        delete n;
    }
}

clt::WardLinkage::WardLinkage(std::vector<std::vector<double> > *_proximity, DistanceType distType, std::vector<std::vector<double> > *_centers)
{
    centers = _centers;
    init(_proximity, _centers->size());
    n = new std::vector<int>(linkageSize, 1);

    if(distType == EUCLIDEAN_DISTANCE)
    {
        for (int i = 0; i < proximity.size(); i++)
        {
            proximity[i] *= proximity[i];
        }
    }
    else if (distType == SQUARED_EUCLIDEAN_DISTANCE)
    {
        ;
    }
}

void clt::WardLinkage::merge(clt::longInt i, clt::longInt j)
{
    float nij = (*n)[i] + (*n)[j];

//    for (int k = 0; k < i; k++)
//    {
//        proximity[index(i, k)] = (d(i, k) * ((*n)[i] + (*n)[k]) + d(j, k) * ((*n)[j] + (*n)[k]) - d(j, i) * (*n)[k]) / (nij + (*n)[k]);
//    }

//    for (int k = i+1; k < j; k++)
//    {
//        proximity[index(k, i)] = (d(k, i) * ((*n)[i] + (*n)[k]) + d(j, k) * ((*n)[j] + (*n)[k]) - d(j, i) * (*n)[k]) / (nij + (*n)[k]);
//    }

//    for (int k = j+1; k <linkageSize; k++)
//    {
//        proximity[index(k, i)] = (d(k, i) * ((*n)[i] + (*n)[k]) + d(k, j) * ((*n)[j] + (*n)[k]) - d(j, i) * (*n)[k]) / (nij + (*n)[k]);
//    }

    for (clt::longInt k = 0; k < i; k++)
    {
        proximity[index(i, k)] = (clt::Utils<double>::squaredEuclideanDistance(&(*centers)[i], &(*centers)[k]) * ((*n)[i] + (*n)[k]) + clt::Utils<double>::squaredEuclideanDistance(&(*centers)[j], &(*centers)[k]) * ((*n)[j] + (*n)[k]) - clt::Utils<double>::squaredEuclideanDistance(&(*centers)[j], &(*centers)[i]) * (*n)[k]) / (nij + (*n)[k]);
    }

    for (clt::longInt k = i+1; k < j; k++)
    {
        proximity[index(k, i)] = (clt::Utils<double>::squaredEuclideanDistance(&(*centers)[k], &(*centers)[i]) * ((*n)[i] + (*n)[k]) + clt::Utils<double>::squaredEuclideanDistance(&(*centers)[j], &(*centers)[k]) * ((*n)[j] + (*n)[k]) - clt::Utils<double>::squaredEuclideanDistance(&(*centers)[j], &(*centers)[i]) * (*n)[k]) / (nij + (*n)[k]);
    }

    for (clt::longInt k = j+1; k <linkageSize; k++)
    {
        proximity[index(k, i)] = (clt::Utils<double>::squaredEuclideanDistance(&(*centers)[k], &(*centers)[i]) * ((*n)[i] + (*n)[k]) + clt::Utils<double>::squaredEuclideanDistance(&(*centers)[k], &(*centers)[j]) * ((*n)[j] + (*n)[k]) - clt::Utils<double>::squaredEuclideanDistance(&(*centers)[j], &(*centers)[i]) * (*n)[k]) / (nij + (*n)[k]);
    }

    (*n)[i] += (*n)[j];
}

//=========--------UPGMCLinkage--------========//

clt::UPGMCLinkage::UPGMCLinkage(std::vector<std::vector<double> > *_proximity)
{
    init(_proximity,_proximity->size());
    n.resize(_proximity->size());
    for (int i = 0; i < n.size(); i++)
    {
        n[i] = 1;
    }

    for (int i = 0; i < proximity.size(); i++)
    {
        proximity[i] *= proximity[i];
    }
}

void clt::UPGMCLinkage::merge(clt::longInt i, clt::longInt j)
{
    float nij = n[i] + n[j];

    for (clt::longInt k = 0; k < i; k++)
    {
        proximity[index(i, k)] = (d(i, k) * n[i] + d(j, k) * n[j] - d(j, i) * n[i] * n[j] / nij) / nij;
    }

    for (clt::longInt k = i+1; k < j; k++)
    {
        proximity[index(k, i)] = (d(k, i) * n[i] + d(j, k) * n[j] - d(j, i) * n[i] * n[j] / nij) / nij;
    }

    for (clt::longInt k = j+1; k <linkageSize; k++)
    {
        proximity[index(k, i)] = (d(k, i) * n[i] + d(k, j) * n[j] - d(j, i) * n[i] * n[j] / nij) / nij;
    }

    n[i] += n[j];
}

//=========--------WPGMCLinkage--------========//

clt::WPGMCLinkage::WPGMCLinkage(std::vector<std::vector<double> > *_proximity)
{
    init(_proximity, _proximity->size());
    for (int i = 0; i < proximity.size(); i++)
    {
        proximity[i] *= proximity[i];
    }
}

void clt::WPGMCLinkage::merge(clt::longInt i, clt::longInt j)
{
    for (clt::longInt k = 0; k < i; k++)
    {
        proximity[index(i, k)] = (d(i, k) + d(j, k)) / 2 - d(j, i) / 4;
    }

    for (clt::longInt k = i+1; k < j; k++)
    {
        proximity[index(k, i)] = (d(k, i) + d(j, k)) / 2 - d(j, i) / 4;
    }

    for (clt::longInt k = j+1; k < linkageSize; k++)
    {
        proximity[index(k, i)] = (d(k, i) + d(k, j)) / 2 - d(j, i) / 4;
    }
}

//=========--------FastPair--------========//

clt::FastPair::FastPair(std::vector<longInt> *_points, clt::Linkage *_linkage)
{
    points = _points;

    linkage = _linkage;

    npoints = _points->size();
    neighbor.resize(npoints);
    index.resize(npoints);
    distance.resize(npoints);

    // Find all neighbors. We use a conga line rather than calling getNeighbor.
    for (clt::longInt i = 0; i < npoints - 1; i++)
    {
        // find neighbor to p[0]
        clt::longInt nbr = i + 1;
        float nbd = FLOAT_MAX_VALUE;

        for (clt::longInt j = i + 1; j < npoints; j++)
        {
            float d = linkage->d((*points)[i], (*points)[j]);
            if (d < nbd)
            {
                nbr = j;
                nbd = d;
            }
        }

        // add that edge, move nbr to points[i+1]
        distance[(*points)[i]] = nbd;
        neighbor[(*points)[i]] = (*points)[nbr];
        (*points)[nbr] = (*points)[i + 1];
        (*points)[i + 1] = neighbor[(*points)[i]];
    }

    // No more neighbors, terminate conga line
    neighbor[(*points)[npoints - 1]] = (*points)[npoints - 1];
    distance[(*points)[npoints - 1]] = FLOAT_MAX_VALUE;

    // set where_are...
    for (int i = 0; i < npoints; i++)
    {
        index[(*points)[i]] = i;
    }
}

void clt::FastPair::findNeighbor(clt::longInt p)
{
    // if no neighbors available, set flag for UpdatePoint to find
    if (npoints == 1)
    {
        neighbor[p] = p;
        distance[p] = FLOAT_MAX_VALUE;
        return;
    }

    // find first point unequal to p itself
    int first = 0;
    if (p == (*points)[first])
    {
        first = 1;
    }

    neighbor[p] = (*points)[first];
    distance[p] = linkage->d(p, neighbor[p]);

    // now test whether each other point is closer
    for (clt::longInt i = first + 1; i < npoints; i++)
    {
        clt::longInt q = (*points)[i];
        if (q != p)
        {
            float d = linkage->d(p, q);
            if (d < distance[p])
            {
                distance[p] = d;
                neighbor[p] = q;
            }
        }
    }
}

void clt::FastPair::add(int p)
{
    findNeighbor(p);
    (*points)[index[p] = npoints++] = p;
}

void clt::FastPair::remove(int p)
{
    npoints--;
    clt::longInt q = index[p];
    index[(*points)[q] = (*points)[npoints]] = q;

    for (clt::longInt i = 0; i < npoints; i++)
    {
        if (neighbor[(*points)[i]] == p)
        {
            findNeighbor((*points)[i]);
        }
    }
}

double clt::FastPair::getNearestPair(std::vector<long long int> &pair)
{
    if (npoints < 2)
        return 0.0;

    double d = distance[(*points)[0]];
    int r = 0;
    for (clt::longInt i = 1; i < npoints; i++) {
        if (distance[(*points)[i]] < d) {
            d = distance[(*points)[i]];
            r = i;
        }
    }

    pair[0] = (*points)[r];
    pair[1] = neighbor[pair[0]];

    if (pair[0] > pair[1]) {
        int t = pair[0];
        pair[0] = pair[1];
        pair[1] = t;
    }

    return d;
}

void clt::FastPair::updatePoint(int p)
{
    neighbor[p] = p;    // flag for not yet found any
    distance[p] = FLOAT_MAX_VALUE;
    for (clt::longInt i = 0; i < npoints; i++)
    {
        clt::longInt q = (*points)[i];
        if (q != p)
        {
            float d = linkage->d(p, q);
            if (d < distance[p])
            {
                distance[p] = d;
                neighbor[p] = q;
            }
            if (neighbor[q] == p)
            {
                if (d > distance[q])
                {
                    findNeighbor(q);
                }
                else
                {
                    distance[q] = d;
                }
            }
        }
    }
}

void clt::FastPair::updateDistance(int p, int q)
{
    float d = linkage->d(p, q);

    if (d < distance[p])
    {
        distance[p] = q;
        neighbor[p] = q;
    }
    else if (neighbor[p] == q && d > distance[p])
    {
        findNeighbor(p);
    }

    if (d < distance[q])
    {
        distance[q] = p;
        neighbor[q] = p;
    }
    else if (neighbor[q] == p && d > distance[q])
    {
        findNeighbor(q);
    }
}

//=========--------HierarchicalClustering--------========//

clt::HierarchicalClustering::HierarchicalClustering(Linkage *linkage)
{
    clt::longInt n = linkage->size();

    merge.resize(n - 1);
    for(clt::longInt i = 0; i < n-1; ++i)
        merge[i].resize(2);

   std::vector<clt::longInt> id(n);
   std::vector<clt::longInt> points(n);

   //initializing from 0 to n-1
   std::iota(id.begin(), id.end(), 0);
   std::iota(points.begin(), points.end(), 0);

   height.resize(n - 1);

   //for test and eit signals
//   hlt::longInt timeToInform = (n-1) / 3;

    FastPair *fp = new FastPair(&points, linkage);

    for (clt::longInt i = 0; i < n - 1; i++)
    {
        height[i] = fp->getNearestPair(merge[i]);
        linkage->merge(merge[i][0], merge[i][1]);      // merge clusters into one
        fp->remove(merge[i][1]);                       // drop b
        fp->updatePoint(merge[i][0]);                  // and tell closest pairs about merger

        int p = merge[i][0];
        int q = merge[i][1];
        merge[i][0] = std::min(id[p], id[q]);
        merge[i][1] = std::max(id[p], id[q]);
        id[p] = n + i;
    }

    if (dynamic_cast<UPGMCLinkage*>(linkage) != nullptr || dynamic_cast<WPGMCLinkage*>(linkage) != nullptr   ||  dynamic_cast<WardLinkage*>(linkage) != nullptr ) {
        for (clt::longInt i = 0; i < height.size(); i++)
        {
            height[i] = std::sqrt(height[i]);
        }
    }

    delete fp;
}

clt::HierarchicalClustering::~HierarchicalClustering()
{
    merge.clear();
    height.clear();
}

std::vector<std::vector<long long int> > * clt::HierarchicalClustering::getTree()
{
    return &merge;
}

std::vector<double> *clt::HierarchicalClustering::getHeight()
{
    return & height;
}

std::vector<int> clt::HierarchicalClustering::partition(int k)
{
    clt::longInt n = merge.size() + 1;
    std::vector<int> membership(n);

    IntHeapSelect *heap = new IntHeapSelect(k);
    for (int i = 2; i <= k; i++) {
        heap->add(merge[n - i][0]);
        heap->add(merge[n - i][1]);
    }

    for (int i = 0; i < k; i++) {
        bfs(membership, heap->get(i), i);
    }

    return membership;
}

std::vector<int> clt::HierarchicalClustering::partition(double h)
{
    for (int i = 0; i < height.size() - 1; i++) {
        if (height[i] > height[i + 1])
            return std::vector<int>();

    }

    int n = (int) merge.size() + 1;
    int k = 2;
    for (; k <= n; k++) {
        if (height[n - k] < h) {
            break;
        }
    }

    if (k <= 2)
        return std::vector<int>();


    return partition(k - 1);
}

void clt::HierarchicalClustering::bfs(std::vector<int> &membership, int cluster, int id)
{
    int n = (int) merge.size() + 1;
   QQueue<int> queue;
    queue.enqueue(cluster);

    while (! queue.empty())
    {
        int i = queue.head();
        queue.dequeue();

        if (i < n)
        {
            membership[i] = id;
            continue;
        }

        i -= n;

        int m1 = merge[i][0];

        if (m1 >= n)
        {
            queue.enqueue(m1);
        }
        else
        {
            membership[m1] = id;
        }

        int m2 = merge[i][1];

        if (m2 >= n)
        {
            queue.enqueue(m2);
        }
        else
        {
            membership[m2] = id;
        }
    }
}

//=========--------Node--------========//

void clt::Node::releaseMemory()
{
    sum.clear();
    for(int i = 0; i < numChildren && i < children.size(); ++i)
    {
        if(children[i])
        {
            children[i]->releaseMemory();
            if(children[i])
                delete children[i];
            children[i] = nullptr;
        }
    }
    children.clear();
}

clt::Node::Node(int *_B, int *_d, float *_T)
{
    n = 0;
    refTo_B = _B;
    refTo_d = _d;
    refTo_T = _T;

    sum.resize(*refTo_d);

    parent = nullptr;
    numChildren = 0;
    children.resize(*refTo_B);
    for(int i = 0; i < *refTo_B; ++i)
        children[i] = nullptr;
}

clt::Node::~Node()
{
    releaseMemory();
}

int clt::Node::getNumChildren()
{
    return numChildren;
}

int clt::Node::getNumObservations()
{
    return n;
}

void clt::Node::setNumObservations(int _n)
{
    n = _n;
}

std::vector<float> *clt::Node::getSumOfObservations()
{
    return &sum;
}

std::vector<clt::Node *> * clt::Node::getChildren()
{
    return &children;
}

double clt::Node::distance(std::vector<float> x)
{
    double dist = 0.0;
    for (int i = 0; i < *refTo_d; i++)
    {
        double _d = sum[i] / n - x[i];
        dist += _d * _d;
    }
    return std::sqrt(dist);
}

double clt::Node::distance(clt::Node *node)
{
    double dist = 0.0;
    for (int i = 0; i < *refTo_d; i++)
    {
        double _d = ( sum[i] / n ) - ( /*(*node->getSumOfObservations() )[i] */node->sum[i] / node->n);
        dist += _d * _d;
    }
    return std::sqrt(dist);
}

clt::Leaf * clt::Node::search(std::vector<float> x)
{
    int index = 0;
    double smallest = children[0]->distance(x);

    // find the closest child node to this data point
    for (int i = 1; i < numChildren; i++)
    {
        double dist = children[i]->distance(x);
        if (dist < smallest)
        {
            index = i;
            smallest = dist;
        }
    }

    if (children[index]->numChildren == 0)
    {
        return (Leaf*) children[index];
    }
    else
    {
        return children[index]->search(x);
    }
}

void clt::Node::update(std::vector<float> x)
{
    n++;
    for (int i = 0; i < *refTo_d; i++) {
        sum[i] += x[i];
    }
}

clt::Node * clt::Node::add(std::vector<float> x, size_t id)
{    
    update(x);

    clt::Node * newRoot = nullptr;

    if(numChildren == 0)
    {
        if( dynamic_cast<Leaf*>(this) != nullptr )
            {
                Leaf *newLeaf = new Leaf(x, refTo_B, refTo_d, refTo_T, id);
                newRoot = add(newLeaf);
            }
    }
    else
    {
        int index = 0;
        double smallest = children[0]->distance(x);

        // find the closest child node to this data point
        for (int i = 1; i < numChildren; i++)
        {
            double dist = children[i]->distance(x);
            if (dist < smallest)
            {
                index = i;
                smallest = dist;
            }
        }

        if( dynamic_cast<Leaf*>(children[index]) != nullptr )
        {
            if (smallest > *refTo_T)
            {
                Leaf *newLeaf = new Leaf(x, refTo_B, refTo_d, refTo_T, id);
                newRoot = add(newLeaf);
            }
            else
                newRoot = children[index]->add(x, id);
        }
        else
            newRoot = children[index]->add(x, id);

    }

    return newRoot;
}

clt::Node * clt::Node::add(Node *node)
{
    clt::Node *newRoot = nullptr;
    if (numChildren < *refTo_B)
    {
        children[numChildren++] = node;
        node->parent = this;
    }
    else
    {
        bool newRootToBeUpdated = false;
        if (parent == nullptr)
        {
            parent = new Node(refTo_B, refTo_d, refTo_T);
            parent->add(this);
            newRoot = parent;
            newRootToBeUpdated = true;
        }
        else
        {
            parent->n = 0;
            int sumSize = (int)parent->sum.size();
            parent->sum.clear();
            parent->sum.resize(sumSize);

        }

        clt::Node *splitedNode = split(node);

        clt::Node *possibleRoot = parent->add(splitedNode);
        if(possibleRoot && !newRootToBeUpdated)
            newRoot = possibleRoot;

        for (int i = 0; i < parent->numChildren; i++)
        {
            parent->n += parent->children[i]->n;
            for (int j = 0; j < *refTo_d; j++)
            {
                parent->sum[j] += parent->children[i]->sum[j];
            }
        }
    }

    return newRoot;
}

clt::Node *clt::Node::split(Node *node)
{
    double farest = 0.0;
    int c1 = 0, c2 = 0;

    std::vector<std::vector<double>> dist(numChildren + 1);
    for(int i = 0; i < numChildren + 1; ++i)
        dist[i].resize(numChildren + 1);

    for (int i = 0; i < numChildren; i++)
    {
        for (int j = i + 1; j < numChildren; j++)
        {
            dist[i][j] = children[i]->distance(children[j]);
            dist[j][i] = dist[i][j];
            if (farest < dist[i][j])
            {
                c1 = i;
                c2 = j;
                farest = dist[i][j];
            }
        }

        dist[i][numChildren] = children[i]->distance(node);
        dist[numChildren][i] = dist[i][numChildren];
        if (farest < dist[i][numChildren])
        {
            c1 = i;
            c2 = numChildren;
            farest = dist[i][numChildren];
        }
    }

    int nc = numChildren;
    std::vector<Node*> child = children;

    // clean up this node.
    for(int i = 0; i < numChildren; ++i)
        children[i] = nullptr;

    numChildren = 0;
    n = 0;
    int sumSize = (int)sum.size();
    sum.clear();
    sum.resize(sumSize);

    Node *brother = new Node(refTo_B, refTo_d, refTo_T);

    for (int i = 0; i < nc; i++)
    {
        if (dist[i][c1] < dist[i][c2])
        {
            add(child[i]);
        }
        else
        {
            brother->add(child[i]);
        }
    }

    if (dist[nc][c1] < dist[nc][c2])
    {
        add(node);
    }
    else
    {
        brother->add(node);
    }

    for (int i = 0; i < numChildren; i++)
    {
        n += children[i]->n;
        for (int j = 0; j < *refTo_d; j++)
        {
            sum[j] += children[i]->sum[j];
        }
    }

    for (int i = 0; i < brother->numChildren; i++)
    {
        brother->n += brother->children[i]->n;
        for (int j = 0; j < *refTo_d; j++)
        {
            brother->sum[j] += brother->children[i]->sum[j];
        }
    }

    return brother;
}

bool clt::Node::leaf()
{
    return numChildren == 0 ? true : false;
}

bool clt::Node::parentOfLeaves()
{
    if(leaf())
        return false;

    for(int i = 0; i < numChildren; ++i)
        if(children[i]->leaf())
            return true;

    return false;
}

//=========--------Leaf--------========//

clt::Leaf::Leaf(std::vector<float> x, int *_B, int *_d, float *_T, size_t _id): Node(_B, _d, _T)
{
    id = _id;
    n = 1;
    for(int i = 0; i < *refTo_d; ++i)
        sum[i] = x[i];
}

clt::Leaf::~Leaf()
{

}

void clt::Leaf::add(std::vector<float> x)
{
    n++;
    for (int i = 0; i < *refTo_d; i++)
    {
        sum[i] += x[i];
    }
}

void clt::Leaf::setClusterLabel(int clLbl)
{
    clusterLabel = clLbl;
}

int clt::Leaf::getClusterLabel()
{
    return clusterLabel;
}

size_t clt::Leaf::getID()
{
    return id;
}

//=========--------BIRCH--------========//

clt::BIRCH::BIRCH(int _d, int _B, float _T)
{
    d = _d;
    B = _B;
    T = _T;
    OUTLIER = INT_MAX_VALUE;
    root = nullptr;        
}

clt::BIRCH::~BIRCH()
{
    centroids.clear();
    if(root)
        delete root;
    root = nullptr;
}

void clt::BIRCH::add(std::vector<float> x, size_t id)
{
    if (root == nullptr)
    {
        root = new Node( &B, &d, &T);
        Leaf *newLeaf = new Leaf(x, &B, &d, &T, id);
        root->add(newLeaf);
        root->update(x);
    }
    else
    {
        clt::Node *newRoot = root->add(x, id);
        if(newRoot)
            root = newRoot;
    }
}

int clt::BIRCH::getBrachingFactor()
{
    return B;
}

float clt::BIRCH::getMaxRadius()
{
    return T;
}

int clt::BIRCH::dimension()
{
    return d;
}

int clt::BIRCH::partition(int k, std::vector<std::vector<size_t> > *&clusters)
{
    return partition(k, 0, clusters);
}

int clt::BIRCH::partition(int k, int minPts, std::vector<std::vector<size_t> > *&clusters)
{

    std::vector<Leaf*> leaves;
    std::vector<Leaf*> OutlierLeaves;
    std::vector<std::vector<double>> centers;
    std::vector<std::vector<double>> OutlierCenters;
    QQueue<Node*> queue;
    queue.enqueue(root);

    while ( !queue.empty())
    {
        Node *node = queue.head();
        queue.dequeue();

        if (node->getNumChildren() == 0)
        {
            if (node->getNumObservations() >= minPts)
            {
                std::vector<double> x(d);
                for (int i0 = 0; i0 < d; i0++)
                {
                    x[i0] = (*node->getSumOfObservations())[i0] / node->getNumObservations();
                }
                centers.push_back(x);
                leaves.push_back((Leaf*)node);
            }
            else
            {
                std::vector<double> x(d);
                for (int i0 = 0; i0 < d; i0++)
                {
                    x[i0] = (*node->getSumOfObservations())[i0] / node->getNumObservations();
                }
                OutlierCenters.push_back(x);

                Leaf *leaf = (Leaf*) node;
                leaf->setClusterLabel(OUTLIER);                
                OutlierLeaves.push_back(leaf);
            }
        }
        else
        {
            for (int i1 = 0; i1 < node->getNumChildren(); i1++)
            {
                queue.enqueue( (*node->getChildren())[i1] );
            }
        }
    }

    clt::longInt totalNodes = leaves.size() + OutlierLeaves.size();
    return totalNodes;
    clt::longInt n = centers.size();

    if (n > k)
    {
        clusters = new  std::vector<std::vector<size_t> > (k,std::vector<size_t>());
        std::vector<std::vector<double>> proximity(n);

        Linkage *linkage = new WardLinkage(&proximity, SQUARED_EUCLIDEAN_DISTANCE, &centers);

        HierarchicalClustering *hc = new HierarchicalClustering(linkage);

        std::vector<int> y = hc->partition(k);

        // recalculate the centroids of each cluster
        centroids.clear();
        centroids.resize(k);
        for(int i2 = 0; i2 < k; ++i2)
            centroids[i2].resize(d);

        std::vector<int> nc(k,0);

        for (int i = 0; i < n; i++)
        {
            Leaf *leaf = leaves[i];
            int yi = y[i];
            //Identifying the cluster
            leaf->setClusterLabel(yi);
            (*clusters)[yi].push_back(leaf->getID());
            nc[yi] += leaf->getNumObservations();
            for (int j = 0; j < d; j++)
            {
                centroids[yi][j] += (*leaf->getSumOfObservations())[j];
            }
        }

        for (int i = 0; i < k; i++)
        {
            for (int j = 0; j < d; j++)
            {
                centroids[i][j] /= nc[i];
            }
        }

        delete linkage;
        delete hc;
    }
    else
    {
        clusters = new  std::vector<std::vector<size_t> > (n,std::vector<size_t>());
        //Identifying the cluster
        for (int i = 0; i < n; i++)
        {
            leaves[i]->setClusterLabel(i);
            (*clusters)[i].push_back(leaves[i]->getID());
        }
    }

    return n;
}

std::vector<std::vector<size_t> > *clt::BIRCH::getBirchClusters()
{
    std::vector<std::vector<size_t>> *clusters = new  std::vector<std::vector<size_t> > ();

    std::vector<Node*> parentOfleaves;
    QQueue<Node*> queue;

    queue.enqueue(root);
    emit notifyIncrementProgressBarValue();
    emit notifyProgressText(tr("Extracting clusters: Looking for clusters on the tree..."));

    while ( !queue.empty())
    {
        Node *node = queue.head();
        queue.dequeue();

        if(node->parentOfLeaves())
        {
            parentOfleaves.push_back(node);
        }
        else
        {
            for (int i1 = 0; i1 < node->getNumChildren(); i1++)
                queue.enqueue( (*node->getChildren())[i1] );
        }
    }


    QString textMsg = tr("Extracting clusters: ") + QString::number(parentOfleaves.size()) + tr(" clusters found.") + tr(" Looking for subclusters...");
    emit notifyIncrementProgressBarValue();
    emit notifyProgressText(textMsg);
    for(quint16 i = 0; i < parentOfleaves.size(); ++i)
    {
        std::vector<size_t> subCluster;
        Node* node = parentOfleaves[i];
        if(node->getNumObservations() == node->getNumChildren())
        {
            for(int j = 0; j < node->getNumChildren(); ++j)
            {
                Leaf *leaf= (Leaf*)(*node->getChildren())[j];
                subCluster.push_back( leaf->getID() );
            }
        }
        else
        {
            queue.clear();
            queue.enqueue(node);
            while ( !queue.empty())
            {
                Node *childNode = queue.head();
                queue.dequeue();

                Leaf *leaf;
                if(dynamic_cast<Leaf*>(childNode) != nullptr)
                {
                    leaf= (Leaf*) childNode;
                    subCluster.push_back( leaf->getID() );
                }
                if(!childNode->leaf())
                {
                    for (int i1 = 0; i1 < childNode->getNumChildren(); i1++)
                        queue.enqueue( (*childNode->getChildren())[i1] );
                }
            }
        }

        clusters->push_back(subCluster);
    }

    textMsg = tr("Extracting clusters: ") + QString::number(clusters->size()) + tr(" clusters found successfully.");
    emit notifyIncrementProgressBarValue();
    emit notifyProgressText(textMsg);

    return clusters;
}

int clt::BIRCH::predict(std::vector<float> x)
{
    if (centroids.empty())
        return -1; //TODO --> verify if this could be a rigth value

    Leaf *leaf = root->search(x);
    return leaf->getClusterLabel();
}

std::vector<std::vector<float>> * clt::BIRCH::getCentroids()
{
    return &centroids;
}

//=========--------Clustering--------========//

clt::Clustering::~Clustering()
{
    resetMemberVariables();
}

clt::Clustering::Clustering() :
    numFeatures(13),
    featureMax(nullptr),
    featureMin(nullptr),
    kClusters(nullptr),
    pointsVector(),
    clusterSize(100)
{

}

void clt::Clustering::resetMemberVariables()
{
    pointsVector.clear();
    if(kClusters)
    {
        delete kClusters;
        kClusters = nullptr;
    }

    if(featureMax)
    {
        featureMax->clear();
        delete featureMax;
        featureMax = nullptr;
    }

    if(featureMin)
    {
        featureMin->clear();
        delete featureMin;
        featureMin = nullptr;
    }

    signalValuesMV.clear();
    signalValuesMV.shrink_to_fit();
}

void clt::Clustering::_classifyDraws()
{
    //beat clustering

    readItems();
    normalizeVector(&pointsVector);

}

void clt::Clustering::_calculateModelDraw()
{

    return ;
}

void clt::Clustering::readItems()
{

    size_t beatSegmentsVectorSize = beatSegments->size();
    int index;

    float x0 = 0, xi = 0, x = 0;
    float y0_mV = 0, yi_mV = 0, y_mV = 0;
    float prevRR = 0, nextRR = 0;
    float cOffset_mV = -5.0f;
    float xOffset = 1000.0;

    std::vector<float> *max = new std::vector<float>(numFeatures, (-1)*FLT_MAX);
    std::vector<float> *min = new std::vector<float>(numFeatures, FLT_MAX);

    for(size_t beatIndex = 0; beatIndex < beatSegmentsVectorSize; ++beatIndex)
    {
        index = 0;
        std::vector<float> vPoints;

        //R-Peak in time, QRS-Onset in amplitude --> point of reference
        x0 = beat->getRPeak(); //Time
        reader.seek(beat->getQRSOnset());
        sample = reader.read();
        y0_mV = static_cast<float>(sample.voltage_mV);  //Amplitude

        //P-Onset
        xi = beat->getPOnset();
        x = x0 - xi + xOffset;
        vPoints.push_back(x);
        updateMinMax(min, max, index++, x);

        //P-Peak
        reader.seek(beat->getPPeak());
        sample = reader.read();
        xi = beat->getPPeak();
        yi_mV = static_cast<float>(sample.voltage_mV);
        x = x0 - xi + xOffset;
        y_mV = y0_mV - yi_mV + cOffset_mV;
        vPoints.push_back(x);
        vPoints.push_back(y_mV);
        updateMinMax(min, max, index++, x);
        updateMinMax(min, max, index++, y_mV);

        //P-Offset
        xi = beat->getPOffset();
        x = x0 - xi + xOffset;
        vPoints.push_back(x);
        updateMinMax(min, max, index++, x);

        //QRS-Onset
        xi = beat->getQRSOnset();
        x = x0 - xi + xOffset;
        vPoints.push_back(x);
        updateMinMax(min, max, index++, x);

        //R-Peak
        reader.seek(beat->getRPeak());
        sample = reader.read();
        yi_mV = static_cast<float>(sample.voltage_mV);
        y_mV = y0_mV - yi_mV + cOffset_mV;
        vPoints.push_back(y_mV);
        updateMinMax(min, max, index++, y_mV);

        //QRS-Offset
        xi = beat->getQRSOffset();
        x = x0 - xi + xOffset;
        vPoints.push_back(x);
        updateMinMax(min, max, index++, x);

        //T-Onset
        reader.seek(beat->getTOnset());
        sample = reader.read();
        xi = beat->getTOnset();
        yi_mV = static_cast<float>(sample.voltage_mV);
        x = x0 - xi + xOffset;
        y_mV = y0_mV - yi_mV + cOffset_mV;
        vPoints.push_back(x);
        vPoints.push_back(y_mV);
        updateMinMax(min, max, index++, x);
        updateMinMax(min, max, index++, y_mV);

        //T-Peak
        reader.seek(beat->getTPeak());
        sample = reader.read();
        xi = beat->getTPeak();
        yi_mV = static_cast<float>(sample.voltage_mV);
        x = x0 - xi + xOffset;
        y_mV = y0_mV - yi_mV + cOffset_mV;
        vPoints.push_back(x);
        vPoints.push_back(y_mV);
        updateMinMax(min, max, index++, x);
        updateMinMax(min, max, index++, y_mV);

        //T-Offset
        xi = beat->getTOffset();
        x = x0 - xi + xOffset;
        vPoints.push_back(x);
        updateMinMax(min, max, index++, x);

        //Previous R-R (Time)
        prevRR = beatIndex == 0 ? 0 : (*beatSegments)[beatIndex]->getRPeak() - (*beatSegments)[beatIndex - 1]->getRPeak();

        //next R-R (Time)
        nextRR = beatIndex + 1 < beatSegmentsVectorSize ? (*beatSegments)[beatIndex + 1]->getRPeak() - (*beatSegments)[beatIndex]->getRPeak() : 0;
        x = prevRR / nextRR;

        vPoints.push_back(x);
        updateMinMax(min, max, index++, x);

        pointsVector.push_back(vPoints);
    }

    featureMin = min;
    featureMax = max;
}

void clt::Clustering::updateMinMax(std::vector<float> *min, std::vector<float> *max, int index, float newValue)
{
    if (newValue > (*max)[index])
    {
        if (!(newValue == INFINITY) && !(newValue == NAN))
        {
            (*max)[index] = newValue;
        }
    }

    if (newValue < (*min)[index])
    {
        if (!(newValue == INFINITY) && !(newValue == NAN))
        {
            (*min)[index] = newValue;
        }
    }
}

void clt::Clustering::normalizeVector(std::vector<std::vector<float> > *vector)
{
    // Normalize all features found in pointsVector to values between [0,1]

    double currFeature;
    float normalizedFeatureValue;

    for(size_t beatIndex = 0; beatIndex < vector->size(); ++beatIndex)
    {
        for (int featureIndex = 0; featureIndex < numFeatures; ++featureIndex)
        {
            currFeature = (*vector)[beatIndex][featureIndex];

            normalizedFeatureValue = currFeature == INFINITY ?
                                                              0.5 :
                                                              ( ((*featureMax)[featureIndex] - (*featureMin)[featureIndex]) == 0 ?
                                                                    0 :
                                                                    (currFeature - (*featureMin)[featureIndex]) / ((*featureMax)[featureIndex] - (*featureMin)[featureIndex]));

            if ((normalizedFeatureValue > 1.0) || (normalizedFeatureValue < 0.0))
            {
                // bad feature value
                Q_ASSERT(false);
            }
            else
            {
                (*vector)[beatIndex][featureIndex] = normalizedFeatureValue;
            }
        }
    }
}

std::vector<std::vector<size_t> > *clt::Clustering::getKClusters(int k)
{

    Q_UNUSED(k)


    std::vector<std::vector<size_t>> *clusters;

    int dimensionality = (int)pointsVector[0].size();
    int branchingFactor = 100;

    BIRCH birch(dimensionality, branchingFactor, 0.005f);



    size_t totalPoints = (int)pointsVector.size();


    for(size_t i = 0; i < totalPoints; ++i)
    {
        birch.add(pointsVector[(int)i], i);
    }


    //BIRCH give us automatically k clusters (unknown the final value)
    clusters = birch.getBirchClusters();

    //uncomment this to use partition method to get k clusters
    //birch.partition(k, clusters);

    return clusters;
}
