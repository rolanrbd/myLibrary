#ifndef CLUSTERING_H
#define CLUSTERING_H
#include <QObject>
#include <QString>
#include <QVector>


using namespace std;

#define FLOAT_MAX_VALUE std::numeric_limits<float>::max();
#define INT_MAX_VALUE std::numeric_limits<int>::max();

namespace clt
{

typedef unsigned long long int longInt;

enum DistanceType
{
    UNSET_DISTANCE_TYPE         = 0,
    EUCLIDEAN_DISTANCE          = 1,
    SQUARED_EUCLIDEAN_DISTANCE  = 2
};

template <class T>
class Utils;
class IntHeapSelect;
class Linkage;
class WardLinkage;
class UPGMCLinkage;
class WPGMCLinkage;
class FastPair;
class HierarchicalClustering;
class Node;
class Leaf;
class TestBeatClassifier2;

template <class T>
class Utils
{
private:
     Utils();
public:
     /**
      * Swap two positions.
      */
     static void swap(std::vector<T> &arr, int i, int j);
    /**
     * To restore the max-heap condition when a node's priority is increased. We move up the heap, exchaning the node at position k with its parent
     * (at postion k/2) if necessary, continuing as long as a[k/2] &lt; a[k] or until we reach the top of the heap.
     */
    static void siftUp(std::vector<T> &arr, int k);
    /**
     * To restore the max-heap condition when a node's priority is decreased. We move down the heap, exchanging the node at position k with the larger
     * of that node's two children if necessary and stopping when the node at k is not smaller than either child or the bottom is reached. Note that
     * if n is even and k is n/2, then the node at k has only one child -- this case must be treated properly.
     */
    static void siftDown(std::vector<T> &arr, int k, int n);
    /**
     * The Euclidean distance.
     */
    static double euclideanDistance(std::vector<T> *x, std::vector<T>  *y);
    static double squaredEuclideanDistance(std::vector<T> *x, std::vector<T>  *y);

    static double centroid(std::vector<T> *x);
    static double radius(std::vector<T> *x);
    static double diameter(std::vector<T> *x);
};

class IntHeapSelect
{
private:
    // The heap size.
    int k;
    // The number of objects that have been added into heap.
    int n;
    // True if the heap is fully sorted.
    bool sorted;
    // The heap array.
    std::vector<int> heap;
public:
    /**
     * Constructor.
     * @param k the heap size.
     */
    IntHeapSelect(int k);
    /**
     * Constructor.
     * @param heap the array to store smallest values to track.
     */
    IntHeapSelect(std::vector<int> &hp);

    ~IntHeapSelect();
    /**
    * Assimilate a new value from the stream.
    */
    void add(int datum);
    /**
     * Returns the k-<i>th</i> smallest value seen so far.
     */
    int peek();
    /**
     * Returns the i-<i>th</i> smallest value seen so far. i = 0 returns the smallest value seen, i = 1 the second largest, ..., i = k-1 the last position
     * tracked. Also, i must be less than the number of previous assimilated.
     */
    int get(int i);
    /**
     * Sort the smallest values.
     */
    void sort();
    /**
     * Place the array in max-heap order. Note that the array is not fully sorted.
     */
    static void heapify(std::vector<int> &arr);
    /**
     * Sorts the specified array into descending order. It is based on Shell
     * sort, which is very efficient because the array is almost sorted by
     * heapifying.
     */
    static void sort(std::vector<int> a, int n);
};

class Linkage
{
protected:
    //The data size.
    clt::longInt linkageSize;
    /**
     * Linearized proximity matrix to store the pair-wise distance measure as dissimilarity between clusters. To save space, we only need the
     * lower half of matrix. And we use float instead of double to save  more space, which also help speed performance. During the
     * clustering, this matrix will be updated to reflect the dissimilarity of merged clusters.
     */
    std::vector<double> proximity;
    /** Initialize the linkage with the lower triangular proximity matrix. */
    void init(std::vector<std::vector<double> > *_proximity, size_t size);
    clt::longInt index(clt::longInt i, clt::longInt j);
public:
    virtual ~Linkage();
    /** Returns the proximity matrix size. */
    clt::longInt size();
    /**
     * Returns the distance/dissimilarity between two clusters/objects, which
     * are indexed by integers.
     */
    float d(clt::longInt i, clt::longInt j);
    /**
     * Merge two clusters into one and update the proximity matrix.
     * @param i cluster id.
     * @param j cluster id.
     */
    virtual void merge(clt::longInt i, clt::longInt j) = 0;
};

class WardLinkage: public Linkage
{
private:
    /**
    * The number of samples in each cluster.
    */
    std::vector<int> *n;
    std::vector<std::vector<double>> *centers;
public:
    ~WardLinkage();
    WardLinkage(std::vector<std::vector<double> > *_proximity, DistanceType distType, std::vector<std::vector<double>> *_centers);

    void merge(longInt i, longInt j) ;
};

class UPGMCLinkage: public Linkage
{
private:
    /**
     * The number of samples in each cluster.
     */
     std::vector<int> n;
public:
     UPGMCLinkage(std::vector<std::vector<double> > *_proximity);
     void merge(longInt i, longInt j);
};

class WPGMCLinkage: public Linkage
{
public:
     WPGMCLinkage(std::vector<std::vector<double> > *_proximity);
     void merge(longInt i, longInt j);
};

class FastPair
{
private:
    std::vector<clt::longInt> *points;     // points currently in set
    std::vector<clt::longInt> index;      // indices into points
    clt::longInt npoints;             // how much of array is actually used?
    std::vector<clt::longInt> neighbor;
    std::vector<float> distance;
    Linkage *linkage;
public:
    FastPair(std::vector<clt::longInt> *_points, Linkage *_linkage);
    /**
     * Find nearest neighbor of a given point.
     */
    void findNeighbor(longInt p);
    /**
     *  Add a point and find its nearest neighbor.
     */
    void add(int p);
    /**
     * Remove a point and update neighbors of points for which it had been nearest
     */
    void remove(int p);
    /**
     * Find closest pair by scanning list of nearest neighbors
     */
    double getNearestPair(std::vector<long long int> &pair);
    /**
     * All distances to point have changed, check if our structures are ok Note that although we completely recompute the neighbors of p,
     * we don't explicitly call findNeighbor, since that would double the number of distance computations made by this routine.
     * Also, like deletion, we don't change any other point's neighbor to p.
     */
    void updatePoint(int p);
    /**
     * Single distance has changed, check if our structures are ok.
     */
    void updateDistance(int p, int q);
};

class HierarchicalClustering
{
private:
    /**
     * An n-1 by 2 matrix of which row i describes the merging of clusters at  step i of the clustering. If an element j in the row is less than n, then
     * observation j was merged at this stage. If j &ge; n then the merge was with the cluster formed at the (earlier) stage j-n of the algorithm.
     */
     std::vector<std::vector<long long int>> merge;
    /**
     * A set of n-1 non-decreasing real values, which are the clustering height, i.e., the value of the criterion associated with the clustering method
     * for the particular agglomeration.
     */
     std::vector<double> height;
public:
    /**
     * Constructor.
     * Learn the Agglomerative Hierarchical Clustering with given linkage method, which includes proximity matrix.
     * @param linkage a linkage method to merge clusters. The linkage object includes the proximity matrix of data.
     */
    HierarchicalClustering(Linkage *linkage);
    ~HierarchicalClustering();
    /**
     * Returns an n-1 by 2 matrix of which row i describes the merging of clusters at  step i of the clustering. If an element j in the row is less than n, then
     * observation j was merged at this stage. If j &ge; n then the merge was with the cluster formed at the (earlier) stage j-n of the algorithm.
     */
    std::vector<std::vector<long long int> > * getTree();
    /**
     * Returns a set of n-1 non-decreasing real values, which are the clustering height, i.e., the value of the criterion associated with the clustering method
     * for the particular agglomeration.
     */
    std::vector<double> * getHeight();
    /**
    * Cuts a tree into several groups by specifying the desired number.
    * @param k the number of clusters.
    * @return the cluster label of each sample.
    */
    std::vector<int> partition(int k);
    /**
     * Cuts a tree into several groups by specifying the cut height.
     * @param h the cut height.
     * @return the cluster label of each sample.
     */
    std::vector<int> partition(double h);
    /**
    * BFS the merge tree and identify cluster membership.
    * @param membership the cluster membership array.
    * @param cluster the last merge point of cluster.
    * @param id the cluster ID.
    */
    void bfs(std::vector<int> &membership, int cluster, int id);
};

template <class T>
class Clustering: public QObject
{
protected:
    /**
    * Cluster label for outliers or noises.
    */
    quint64 OUTLIER;

    /**
    * Cluster a new instance.
    * @param x a new instance.
    * @return the cluster label.
    */
public:
    virtual ~Clustering(){}
    virtual int predict(T x) = 0;
};

class Node
{
protected:
    // Refenrence to the Branching factor. Maximum number of children nodes.
    int* refTo_B;
    // Reference to the dimensionality of data.
    int* refTo_d;
    // Refenrece to the Maximum radius of a sub-cluster.
    float* refTo_T;
    // The number of observations
    int n;
    // The sum of the observations
    std::vector<float> sum;
    // The number of children.
    int numChildren;
    // Children nodes.
    std::vector<Node*> children;
    // Parent node.
    Node * parent;
    void releaseMemory();
public:
    Node(int* _B, int* _d, float *_T);
    virtual ~Node();

    int getNumChildren();
    int getNumObservations();
    void setNumObservations(int _n);
    std::vector<float> * getSumOfObservations();
    std::vector<Node*> * getChildren();

    //Calculates the distance between x and CF center
    double distance(std::vector<float> x);
    //Calculates the distance between CF centers
    double distance(Node *node);
    //Returns the leaf node closest to the given data
    Leaf * search(std::vector<float> x);
    //Adds data to this node.
    void update(std::vector<float> x);
    //Adds data to the node
    Node *add(std::vector<float> x, size_t id);
    //Add a node as children. Split this node if the number of children reach the Branch Factor
    Node *add(Node *node);
    //Split the node and return a new node to add into the parent
    Node * split(Node *node);
    bool leaf();
    bool parentOfLeaves();
};

class Leaf: public Node
{
private:
    // To identify the node
    size_t id;
    //The cluster label of the leaf node.
    int clusterLabel;
public:
    Leaf(std::vector<float> x, int* _B, int* _d, float* _T, size_t _id);
    ~Leaf();

    void add(std::vector<float> x);
    void setClusterLabel(int clLbl);

    int getClusterLabel();
    size_t getID();
};

class BIRCH: public Clustering<std::vector<float>>
{
    Q_OBJECT

    friend class Node;
private:
    // Branching factor. Maximum number of children nodes.
    int B;
    // Maximum radius of a sub-cluster.
    float T;
    // The dimensionality of data.
    int d;
    // The root of CF tree.
    Node *root;
    // Leaves of CF tree as representatives of all data points.
    std::vector<std::vector<float>> centroids;
    void updateRoot();
    void outlierAnalysis(std::vector<std::vector<size_t> > *&clusters, std::vector<Leaf*> *outliers, std::vector<std::vector<double>> *OutlierCenters);
public:
    /**
    * Constructor.
    * @param d the dimensionality of data.
    * @param B the branching factor. Maximum number of children nodes.
    * @param T the maximum radius of a sub-cluster.
    */
    BIRCH(int _d, int _B, float _T);

    ~BIRCH();

    //Add a data point into CF tree.
    void add(std::vector<float> x, size_t id);
    //Returns the branching factor, which is the maximum number of children nodes.
    int getBrachingFactor();
    //Returns the dimensionality of data.
    float getMaxRadius();
    //Returns the dimensionality of data.
    int dimension();
    /**
     * Clustering leaves of CF tree into k clusters.
     * @param k the number of clusters.
     * @return the number of non-outlier leaves.
     */
    int partition(int k, std::vector<std::vector<size_t> > *&clusters);
    /**
    * Clustering leaves of CF tree into k clusters.
    * @param k the number of clusters.
    * @param minPts a CF leaf will be treated as outlier if the number of its points is less than minPts.
    * @return the number of non-outlier leaves.
    */
    int partition(int k, int minPts, std::vector<std::vector<size_t> > *&clusters);
    /**
     * Cluster a new instance to the nearest CF leaf. After building the CF tree, the user should call {@link #partition(int)} method first
     * to clustering leaves. Then they call this method to clustering new data.
     * @param x a new instance.
     * @return the cluster label, which is the label of nearest CF leaf.
     * Note that it may be {@link #OUTLIER}.
     */
    int predict(std::vector<float> x);

    std::vector<std::vector<size_t>> * getBirchClusters();

    /**
     * Returns the representatives of clusters.
     * @return the representatives of clusters
     */
    std::vector<std::vector<float> > *getCentroids();

signals:
    void notifyIncrementProgressBarValue();
    void notifyProgressText(QString text);
};

class Clustering
{
    //Members
private:
    std::vector<std::vector<float>> pointsVector;
    std::vector<float> *featureMax;
    std::vector<float> *featureMin;
    // Branching factor. Maximum number of children nodes.
    int B;
    // Maximum radius of a sub-cluster.
    float T;
    // The dimensionality of data. This value will be (d) in the birch algorithm
    int numFeatures;
    // Amount of clusters required
    int clusterSize;

    std::vector<float> signalValuesMV;

    std::vector<std::vector<size_t>> *kClusters;

    //Methods
public:
    ~Clustering();

protected:
    Clustering();

    virtual void resetMemberVariables();

    virtual void _classifyDraws();
    virtual void _calculateModelDraw();
private:
    void readItems();
    void updateMinMax(std::vector<float> *min, std::vector<float> *max, int index, float newValue);
    void normalizeVector(std::vector<std::vector<float>> *vector);

    std::vector<std::vector<size_t> > *getKClusters(int k);

};

} // end of namespace hlt

#endif // CLUSTERING_H
