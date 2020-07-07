#ifndef COMBINATIONSTREE_H
#define COMBINATIONSTREE_H

#include <QVector>
#include <QList>
#include <QPair>
#include <QDebug>

template<class T>
class GTree
{
private:
    quint16 level;
    T value;
    QString key;
    QList<GTree<T>*> subTrees;
    GTree<T> *parent;

    void getAllBranches(QList<QList<QPair<T, QString>>> &path)
    {
        if(leaf())
        {
            QList<QPair<T, QString>> &lastPath = path.last();
            QPair<T, QString> p(value, key);
            lastPath.push_back(p);
            return;
        }
        else
        {
            QList<QPair<T, QString>> &lastPath = path.last();
            QPair<T, QString> p(value, key);
            lastPath.push_back(p);

            QList<QPair<T, QString>> tempPath =  path.last();

            for(int i = 0; i < subTrees.size(); i++)
            {
                T x = subTrees[i]->getValue();
                if(subTrees[i]->leaf())
                {
                    QList<QPair<T, QString>> &lastPath = path.last();
                    QPair<T, QString> p(subTrees[i]->getValue(), subTrees[i]->getKey());
                    lastPath.push_back(p);

                    if( i+1 <  subTrees.size())
                    {
                        QList<QPair<T, QString>> newPath;
                        std::copy(tempPath.begin(), tempPath.end(), std::back_inserter(newPath));
                        path.push_back(newPath);
                    }
                }
                else
                    subTrees[i]->getAllBranches(path);
            }
        }
    }
    void setParent(GTree<T> *root){parent = root;}

    GTree(QPair<T, QString> x, quint16 _level)
    {
        level = _level;
        value = x.first;
        key = x.second;
        parent = nullptr;
    }

public:
    GTree()
    {
        parent = nullptr;
        key = "";
    }

    GTree(T x, quint16 _level = 0)
    {
        level = _level;
        value = x;
        parent = nullptr;
    }

    quint16 getLevel(){return level;}

    T getValue(){return value;}

    QString getKey(){return key;}

    GTree<T> * getParent(){return parent;}

    bool leaf(){return subTrees.empty();}

    bool addItem(quint16 _level, QPair<T, QString> root, T x)
    {
        if(level == _level - 1)
        {
            if(key == root.second)
            {
                QString _key = root.second == "" ? QString::number(x) : root.second + "-" + QString::number(x);
                QPair<T, QString> newItem(x,_key );
                GTree<T>* newNode = new GTree<T>(newItem, _level);
                newNode->setParent(this);
                subTrees.push_back( newNode );
                return true;
            }
            return false;
        }
        else
        {
            for(int i = 0; i < subTrees.size(); i++)
                if(subTrees[i]->addItem(_level, root, x))
                    return true;
            return false;
        }

        return false;
    }

    void addChild(GTree<T> *leaf, quint16 _level, QPair<T, QString> root, T x)
    {

        QString _key = root.second == "" ? QString::number(x) : root.second + "-" + QString::number(x);
        QPair<T, QString> newItem(x,_key );
        GTree<T>* newNode = new GTree<T>(newItem, _level);
        newNode->setParent(leaf);
        leaf->subTrees.push_back( newNode );
    }
    QList<QList<QPair<T, QString> > > getAllBranches()
    {
        QList<QList<QPair<T, QString>>> branches;

        if(subTrees.empty())
            return branches;

        QList<QList<QPair<T, QString>>> allPaths;

        for(int i = 0; i < subTrees.size(); i++)
        {
            QList<QPair<T, QString>> l;
            allPaths.push_back(l);
            subTrees[i]->getAllBranches(allPaths);
        }

        return allPaths;
    }

    QList<QList<T>>     getBranches(T root, quint16 startLevel = 0)
    {
        QList<QList<T>> branches;

        if(level == startLevel - 1)
        {
            if(subTrees.empty())
                return branches;

            for(int i = 0; i < subTrees.size(); i++)
            {
                if(subTrees.at(i)->getValue() == root)
                {
                    QList<QList<T>> subBranches = subTrees.at(i)->getAllBranches();

                    T x = subTrees.at(i)->getValue();

                    for(int j = 0; j < subBranches.size(); j++)
                        for(int k = 0; k < subBranches[j].size(); k++)
                        {
                            branches[i].push_back(x);
                            branches[i].push_back(subBranches.at(j).at(k));
                        }
                }
            }

        }
        return branches;
    }

    QList<GTree<T>*>&   getSubTrees()
    {
        return &subTrees;
    }

    void getAllLeaves(QList<QString> &path)
    {
        if(leaf())
        {
            path.push_back(key);
            return;
        }
        else
        {
            for(int i = 0; i < subTrees.size(); i++)
                subTrees[i]->getAllLeaves(path);
        }
    }

    void getLeavesRefAndPath(QList<GTree<T>*> &leavesRef, QList<QString> &path)
    {        if(leaf())
        {
            path.push_back(key);
            leavesRef.push_back(this);
            return;
        }
        else
        {
            for(int i = 0; i < subTrees.size(); i++)
                subTrees[i]->getLeavesRefAndPath(leavesRef, path);
        }
    }
    quint16 childSize(){ return subTrees.size();}


    Combinations * getAllCombinations(quint16 size, bool ascedOrder = true)
    {
        Q_UNUSED(ascedOrder)
        Combinations *cmb = new Combinations();

        QList<QString> subResults;
        getAllLeaves(subResults);

        for(quint16 i = 0; i < subResults.size(); ++i)
        {
           LotteryPlay *l = new LotteryPlay();

            QString s = subResults[i];
            QStringList numbers = s.split("-");

            if(numbers.size() != size)
                continue;

            QVector<quint16> comb;
            for(quint16 j = 0; j < size; ++j)
            {
                quint16 n = numbers[j].toInt();
                l->push_back(n);
            }

            if(l && l->size() != size)
            {
                l->clear();
                delete l;
                l = nullptr;
                continue;
            }
            else if(!l)
                continue;

            cmb->push_back(l);
        }

        subResults.clear();

        return cmb;
    }
//    GTree<T>*  getSubTree(int);
//    QList<T>  preOrderTree();
//    GTree<T>*  postOrderTree();
//    GTree<T>*  betweenOrderTree();
//    GTree<T>*  broadwaysTree();
};

#endif // COMBINATIONSTREE_H
