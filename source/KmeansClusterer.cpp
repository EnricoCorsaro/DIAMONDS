#include "KmeansClusterer.h"



// KmeansClusterer::KmeansClusterer()
//
// PURPOSE: 
//      Constructor.
//
// INPUT: 
//      metric: this class is used to compute the squared distance between two points
//      minNclusters: the minimum number of clusters that is to be fitted
//      maxNclusters: the maximum number of clusters that is to be fitted
//      Ntrials: as k-means is sensitive to the initial choice of the cluster centers, 
//               repeat k-means with different trials of initial centers, and pick the best one.
//      relTolerance: fitting the clusters converged when S < relTolerance, where S is the relative 
//                    change in total sum of distances of all points to their cluster center. 
// 

KmeansClusterer::KmeansClusterer(Metric &metric, unsigned int minNclusters, unsigned int maxNclusters, unsigned int Ntrials, double relTolerance)
: Clusterer(metric), 
  minNclusters(minNclusters), 
  maxNclusters(maxNclusters), 
  Ntrials(Ntrials), 
  relTolerance(relTolerance),
  engine(time(0))
{
    assert(minNclusters <= maxNclusters);
}









// KmeansClusterer::~KmeansClusterer()
//
// PURPOSE: 
//      Destructor.
//

KmeansClusterer::~KmeansClusterer()
{

}








// KmeansClusterer::chooseInitialClusterCenters()
//
// PURPOSE: 
//      Choose semi-randomly the initial centers of the clusters among the sample
//      of points given
//
// INPUT:
//      sample(Ndimensions, Npoints): sample of N-dimensional points
//      center(Ndimensions, Nclusters): set of N-dimensional coordinates of the cluster centers 
//      Nclusters: number of clusters considered
// 
// OUTPUT: 
//      void
//

void KmeansClusterer::chooseInitialClusterCenters(RefArrayXXd sample, RefArrayXXd centers, unsigned int Nclusters)
{
    unsigned int Npoints = sample.cols();


    // Set up a some random generators
    
    uniform_real_distribution<> uniform01(0.0, 1.0);
    uniform_int_distribution<> uniform(0, Npoints-1);


    // Picking the initial centers randomly is prone to leading to a local rather than 
    // the global minumum. Choosing k centers as far away from each other as possible 
    // is prone to outliers. We therefore adopt the method of Arthur & Vassilvitskii (2007)
    // which instead of always choosing a new center farthest from those picked so far, 
    // chooses each center at random with a probability proportional to its (squared) 
    // distance from the centers chosen already.


    // Choose the first center randomly

    int randomPointIndex = uniform(engine);
    int k;
    double distanceToClosestCenter[Npoints];
    double sumOfDistancesToClosestCenters;
    double distance;
    double uniform01Number;
    double cumulativeDistance;
    centers.col(0) = sample.col(randomPointIndex);

    
    // Select the other initial centers probabilistically 

    for (int n = 1; n < Nclusters; ++n)
    {
        // For each of the points in the sample, determine the distance to
        // its closest center 
    
        sumOfDistancesToClosestCenters = 0.0;
    
        for (int k = 0; k < Npoints; ++k)
        {
            for (int j = 0; j < n; ++j)
            {
                distance = metric.distance(sample.col(k), centers.col(j));
                if (j == 0)
                {
                    // This is the very first center we're checking, so let's adopt
                    // this initially as the closest center. 
                
                    distanceToClosestCenter[k] = distance;
                }
                else
                {
                    // Compare with the previous distance to the closest center.
                    // If it is even smaller, we've found an even closer center.
                    // In that case, keep the distance value.
                
                    if (distance < distanceToClosestCenter[k])
                    {
                        distanceToClosestCenter[k] = distance;
                    }
                }
            }
        
            sumOfDistancesToClosestCenters += distanceToClosestCenter[k];
        } 
    

        // Normalize the distances
    
        for (int k = 0; k < Npoints; ++k)
        { 
            distanceToClosestCenter[k] /= sumOfDistancesToClosestCenters;
        }
    

        // Generate a uniform random number between 0 and 1
    
        uniform01Number = uniform01(engine);
    

        // Select the point that makes the cumulative distance greater than the random
        // number. Those points with a larger distance to their closest center, will have 
        // a greater chance to be chosen as the next cluster center point, than the others.
        
        cumulativeDistance = distanceToClosestCenter[0];
        k = 0;
        
        while (cumulativeDistance < uniform01Number)
        {
            k++;
            cumulativeDistance += distanceToClosestCenter[k];
        }

        centers.col(n) = sample.col(k);

    } // end loop of selecting initial cluster centers
}












// KmeansClusterer::updateClusterCentersUntilConverged()
//
// PURPOSE: 
//      Evolve the cluster centers according to the k-means algorithm. Given a set of cluster
//      centers, gather for each center the sample points that are closest to it, and use 
//      the barycenter of these gathered points as the updated cluster center. And so on.
//
// INPUT:
//      sample(Ndimensions, Npoints): sample of N-dimensional points
//      center(Ndimensions, Nclusters): set of N-dimensional coordinates of the cluster centers 
//      clusterSizes(Nclusters): for each cluster, the number of points it contains
//      clusterIndices(Npoints): for each point the index of the cluster it belongs to. This index
//                               runs from 0 to Nclusters-1.
//      sumOfDistancesToClosestCenter: the sum over all points of their distance to the closest 
//                                     cluster center (i.e. the center of the cluster to which they
//                                     belong).
//      relTolerance: fitting the clusters converged when S < relTolerance, where S is the relative 
//                    change in total sum of distances of all points to their cluster center.  
// 
// OUTPUT:
//      True if the convergence of k-means was successfully, false otherwise. A convergence is
//      considered successfully if the relTolerance-criterion is satisfied (see above), _and_ if 
//      all clusters contain at least 2 points.


bool KmeansClusterer::updateClusterCentersUntilConverged(RefArrayXXd sample, RefArrayXXd centers, 
                                                         RefArrayXd clusterSizes, vector<int> &clusterIndices,
                                                         double &sumOfDistancesToClosestCenter, double relTolerance)
{
    unsigned int Npoints = sample.cols();
    unsigned int Ndimensions = sample.rows();
    unsigned int Nclusters = centers.cols();
    ArrayXXd updatedCenters = ArrayXXd::Zero(Ndimensions, Nclusters);   // coordinates of each of the new cluster centers


    // Perform the k-means clustering iteration, each time improving the cluster centers,
    // and redetermining which points belongs to which cluster
    
    bool stopIterations = false;
    bool convergenceReached;
    unsigned int indexOfClosestCenter;
    double oldSumOfDistances = 0.0;
    double newSumOfDistances = 0.0;
    double distanceToClosestCenter;
    double distance; 

    while (!stopIterations)
    {
        // Find for each point the closest cluster center.
        // At the same time recompute/update the new cluster centers, which is simply
        // the barycenter of all points belonging to the cluster. 
    
        clusterSizes.setZero();
        updatedCenters.setZero();
    
        for (int n = 0; n < Npoints; ++n)
        {
            distanceToClosestCenter = DBL_MAX;
        
            for (int i = 0; i < Nclusters; ++i)
            {
                distance = metric.distance(sample.col(n), centers.col(i));
                
                if (distance < distanceToClosestCenter)
                {
                    indexOfClosestCenter = i;
                    distanceToClosestCenter = distance;
                }
            }
        
            newSumOfDistances += distanceToClosestCenter;
            updatedCenters.col(indexOfClosestCenter) += sample.col(n);
            clusterSizes(indexOfClosestCenter) += 1; 
            clusterIndices[n] = indexOfClosestCenter;        
        }
    

        // Assert that all clusters contain at least 2 points. If not we probably started
        // with an unfortunate set of initial cluster centers. Flag this by immediately 
        // returning false.
       
        if (!(clusterSizes > 1).all())
        {
            convergenceReached = false;
            return convergenceReached;
        }
       

        // Finish computing the new updated centers. Given the check above, we are sure
        // that none of the clusters is empty. 
        
        updatedCenters.rowwise() /= clusterSizes.transpose();
        centers = updatedCenters;
    

        // A new set of clusters has been determined.
        // Decide whether the algorithm has converged. Convergence occurs when
        // the sum of all distances of all points to their cluster center does
        // not change significantly anymore. 
        // Note: in order for this criterion to work properly, the coordinate
        //       space should be normalized, so that one particular coordinate
        //       cannot numerically dominate all other coordinates.
    
        if (oldSumOfDistances == 0.0)
        {
            // This is the first center-updating iteration, so there is nothing to compare yet.
            // Simply set the variables.
    
            oldSumOfDistances = newSumOfDistances;
            newSumOfDistances = 0.0;
        }
        else
        {
            // If the relative change in sumOfDistances between old and new was smaller than
            // the threshold set by the user, stop the iteration loop.
        
            if (fabs(newSumOfDistances - oldSumOfDistances) / oldSumOfDistances < relTolerance)
            {
                sumOfDistancesToClosestCenter = newSumOfDistances;      // will be returned to user
                stopIterations = true;
            }
            else
            {
                oldSumOfDistances = newSumOfDistances;
                newSumOfDistances = 0.0;
            }   
        }
    }  // end k-means center-updating loop 
    
    
    // Convergence was properly reached, so return
    
    convergenceReached = true;
    
    return convergenceReached;  
}











// KmeansClusterer::evaluateBICvalue()
//
// PURPOSE:
//      Given a cluster configuration, evaluate its Bayesian Information Criterion value
//      assuming that the clusters are a gaussian mixture of spherical Gaussians.
//
// INPUT:
//      sample(Ndimensions, Npoints): sample of N-dimensional points
//      center(Ndimensions, Nclusters): set of N-dimensional coordinates of the cluster centers 
//      clusterSizes(Nclusters): for each cluster, the number of points it contains
//      clusterIndices(Npoints): for each point the index of the cluster it belongs to. This index
//                               runs from 0 to Nclusters-1.
// 
// OUTPUT:
//      BIC value: - 2 ln(L) + k ln(N) up to constant terms
//

double KmeansClusterer::evaluateBICvalue(RefArrayXXd sample, RefArrayXXd centers, 
                                         RefArrayXd clusterSizes, vector<int> &clusterIndices)
{
    unsigned int Npoints = sample.cols();
    unsigned int Ndimensions = sample.rows();
    unsigned int Nclusters = centers.cols();
        

    // Compute the intra-cluster variance for each cluster, assuming that each cluster 
    // is spherical, making it a one-dimensional problem.
    
    ArrayXd intraClusterVariances(Nclusters);
    intraClusterVariances.setZero();
    
    for (int n = 0; n < Npoints; ++n)
    {
        intraClusterVariances(clusterIndices[n]) += metric.distance(sample.col(n), centers.col(clusterIndices[n]));
    }
    
    intraClusterVariances /= (clusterSizes-1); 
    

    // Initialize the cluster priors, i.e. the prior probability that a point belongs 
    // to a particular cluster, which we set to the relative cluster size.
    
    ArrayXd clusterPriors = clusterSizes / double(Npoints);
   

    // Compute the log-likelihood (up to constant terms) which we assume is a mixture
    // of spherical gaussians. The log-likelihood of the entire sample is the sum of the
    // log-likelihood of each of the points.
    
    double logLikelihood = 0.0;
    int cluster;

    for (int n = 0; n < Npoints; ++n)
    {
        cluster = clusterIndices[n]; 
        logLikelihood +=   log(clusterPriors(cluster))
                         - Ndimensions / 2.0 * log(intraClusterVariances(cluster))
                         - metric.distance(sample.col(n), centers.col(cluster)) / 2.0 / intraClusterVariances(cluster);
    } 
    
    
    // The number of free parameters is the sum of 
    // 1) Ndimensions*Nclusters center coordinates
    // 2) Nclusters intra-cluster variances
    // 3) Nclusters-1 cluster priors (-1 because the priors should add up to one).
            
    unsigned int NfreeParameters = Ndimensions * Nclusters + Nclusters + Nclusters-1;
    
    
    // Compute the BIC value of the current cluster configuration, and compare it to
    // the best BIC value we currently have. The lower the BIC value, the better the
    // model. Store the better cluster configuration.
    
    double BICvalue = -2.0 * logLikelihood + NfreeParameters * log(Npoints);

    return BICvalue;
}












// KmeansClusterer::cluster()
//
// PURPOSE: 
//      Given a sample of N-dimensional points, use the k-means clustering algorithm + BIC
//      to determine the optimal number of clusters that can be fitted to the data.
//
// INPUT:
//      printOnTheScreen: a boolean value specifying whether the BIC values and corresponding number of clusters 
//                        are to be printed on the screen while the process is running.
//      sample(Ndimensions, Npoints): sample of N-dimensional points
//      optimalClusterIndices(Npoints): for each point the index of the cluster it belongs to. This index
//                                      runs from 0 to Nclusters-1.
// 
// OUTPUT:
//      The optimal number of clusters
//

int KmeansClusterer::cluster(const bool printOnTheScreen, RefArrayXXd sample, vector<int> &optimalClusterIndices, const bool sortSample)
{
    bool convergedSuccessfully;
    unsigned int Npoints = sample.cols();
    unsigned int Ndimensions = sample.rows();
    unsigned int optimalNclusters;    
    double bestBICvalue = DBL_MAX;
    double BICvalue; 
    double sumOfDistancesToClosestCenter;
    double bestSumOfDistancesToClosestCenter = DBL_MAX;
    vector<int> clusterIndices(Npoints);            // for each point the index of the cluster to ...
    vector<int> bestClusterIndices(Npoints);        // ... which it belongs
    ArrayXd clusterSizes;
    ArrayXd bestClusterSizes;
    ArrayXXd centers;
    ArrayXXd bestCenters;
    

    // As we don't know a prior the optimal number of clusters, loop over a
    // user-specified range of clusters, and determine which number gives the
    // optimal clustering
    
    if (printOnTheScreen)
    {
        cerr << "=========================================" << endl; 
        cerr << "Information on X-means clustering" << endl;
        cerr << "=========================================" << endl; 
    }

    for (unsigned int Nclusters = minNclusters; Nclusters <= maxNclusters; ++Nclusters)
    {
        centers = ArrayXXd::Zero(Ndimensions, Nclusters);          // coordinates of each of the old cluster centers
        bestCenters = ArrayXXd::Zero(Ndimensions, Nclusters);      // coordinates of the best centers (over all trials)
        clusterSizes = ArrayXd::Zero(Nclusters);                    // # of points belonging to each cluster...
        bestClusterSizes = ArrayXd::Zero(Nclusters);                // ... 'double', to avoid casting problems.               
       

        // The k-means algorithm is sensitive to the choice of the initial centers. 
        // We therefore run the algorithm 'Ntrial' times, and take the best clustering.
        
        bestSumOfDistancesToClosestCenter = DBL_MAX;
        
                    
        for (int m = 0; m < Ntrials; ++m)
        {
            chooseInitialClusterCenters(sample, centers, Nclusters);
            convergedSuccessfully = updateClusterCentersUntilConverged(sample, centers, clusterSizes, clusterIndices, 
                                                                       sumOfDistancesToClosestCenter, relTolerance);
   

            // If the convergence was not successfull (e.g. because some clusters contain 0 or 1 points),
            // we likely had an unfortunate set of initial cluster centers. In this case, simply continue
            // with the next 'trial'.
            
            if (!convergedSuccessfully) continue;
   

            // If we did obtain a successful convergence, compare it with the previous clusterings 
            // (all of them with the same number of clusters), and keep the best one.
            
            if (sumOfDistancesToClosestCenter < bestSumOfDistancesToClosestCenter)
            {
                bestSumOfDistancesToClosestCenter = sumOfDistancesToClosestCenter;
                bestCenters = centers;
                bestClusterIndices = clusterIndices;
                bestClusterSizes = clusterSizes;  
            }               
        } // end loop over Ntrials to determine the best clustering trying different initial centers
       

        // Evaluate the current number of clusters, using the BIC value. Note that this is only necessary 
        // if the user selected more than one particular number of clusters.

        if (maxNclusters - minNclusters > 1)
        {
            BICvalue = evaluateBICvalue(sample, bestCenters, bestClusterSizes, bestClusterIndices);
            
            if (printOnTheScreen)
                cerr << "Nclusters: " << Nclusters << " " << "BIC: " << BICvalue << endl;
                        
            if (BICvalue < bestBICvalue)
            {
                // In what follows 'best' refers to the best cluster configuration given a specific
                // number of clusters. 'optimal' refers to the optimal configuration over all possible
                // values for the number of clusters.
                
                bestBICvalue = BICvalue;
                optimalNclusters = Nclusters;
                optimalClusterIndices = bestClusterIndices;
            }
        }
        else
        {
            // User allowed only 1 particular number of clusters. Computation of BIC
            // to compare is no longer required.
                         
            optimalNclusters = Nclusters;
            optimalClusterIndices = bestClusterIndices;
        }
        
          
    } // end loop over Nclusters
    
    // That's it!

    if (printOnTheScreen)
    {
        cerr << "Optimal Nclusters: " << optimalNclusters << endl;
        cerr << endl;
    }

    return optimalNclusters;
}

