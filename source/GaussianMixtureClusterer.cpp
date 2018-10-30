#include "GaussianMixtureClusterer.h"



// GaussianMixtureClusterer::GaussianMixtureClusterer()
//
// PURPOSE: 
//      Constructor.
//
// INPUT: 
//      metric: this class is used to compute the squared distance between two points
//      featureProjector: this class is used to compute a dimensionality reduction of the input sample of points
//      featureProjectionActivated: a boolean to activate or deactivate the dimensionality reduction through the featureProjector class
//      minNclusters: the minimum number of clusters that is to be fitted
//      maxNclusters: the maximum number of clusters that is to be fitted
//      Ntrials: as k-means is sensitive to the initial choice of the cluster centers, 
//               repeat k-means with different trials of initial centers, and pick the best one.
//      relTolerance: fitting the clusters converged when dp < relTolerance, where dp is the relative 
//                    change in the total likelihood of the Gaussian Mixture Model as a function of the iterations. 
// 

GaussianMixtureClusterer::GaussianMixtureClusterer(Metric &metric, Projector &featureProjector, bool featureProjectionActivated, 
unsigned int minNclusters, unsigned int maxNclusters, unsigned int Ntrials, 
double relTolerance)
: Clusterer(metric, featureProjector, featureProjectionActivated),
  minNclusters(minNclusters), 
  maxNclusters(maxNclusters), 
  Ntrials(Ntrials), 
  relTolerance(relTolerance)
{
    // Set the seed of the random generator using the clock

    clock_t clockticks = clock();
    engine.seed(clockticks);

    // Do some sanity check(s)

    assert(minNclusters <= maxNclusters);
}









// GaussianMixtureClusterer::~GaussianMixtureClusterer()
//
// PURPOSE: 
//      Destructor.
//

GaussianMixtureClusterer::~GaussianMixtureClusterer()
{

}








// GaussianMixtureClusterer::chooseInitialClusterCenters()
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

void GaussianMixtureClusterer::chooseInitialClusterCenters(RefArrayXXd sample)
{
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










// GaussianMixtureClusterer::chooseInitialClusterCovarianceMatrices()
//
// PURPOSE: 
//      Choose the initial covariance matrices for each cluster by using the set of cluster
//      centers found previously.
//
// INPUT:
//      sample(Ndimensions, Npoints): sample of N-dimensional points
//      center(Ndimensions, Nclusters): set of N-dimensional coordinates of the cluster centers 
//      covarianceMatrices(Ndimensions,Nclusters*Ndimensions): set of Nclusters N-dimensional covariance matrices for the clusters
//      Nclusters: number of clusters considered
// 
// OUTPUT: 
//      void
//
void GaussianMixtureClusterer::chooseInitialClusterCovarianceMatrices(RefArrayXXd sample)
{
    double biasFactor = 1.0/(Npoints-1.0);

    for (int i = 0; i < Nclusters; i++)
    {
        differenceFromCenters.block(0,i*Npoints,Ndimensions,Npoints) = sample.colwise() - centers.col(i);

        covarianceMatrices.block(0, i*Ndimensions,Ndimensions,Ndimensions) = 
                differenceFromCenters.block(0, i*Npoints,Ndimensions,Npoints).matrix() * 
                differenceFromCenters.block(0, i*Npoints,Ndimensions,Npoints).matrix().transpose() * biasFactor;
        
        determinantOfCovarianceMatrices(i) = covarianceMatrices.block(0, i*Ndimensions,Ndimensions,Ndimensions).matrix().determinant();
       
        inverseOfCovarianceMatrices.block(0, i*Ndimensions,Ndimensions,Ndimensions).matrix() = 
                covarianceMatrices.block(0, i*Ndimensions,Ndimensions,Ndimensions).matrix().inverse();
    }
}








// GaussianMixtureClusterer::computeGaussianMixtureModel()
//
// PURPOSE: 
//      Choose the initial covariance matrices for each cluster by using the set of cluster
//      centers found previously. Assume amplitudes of individual Gaussians are defined from outside.
//
// INPUT:
//      sample(Ndimensions, Npoints): sample of N-dimensional points
//      center(Ndimensions, Nclusters): set of N-dimensional coordinates of the cluster centers 
//      covarianceMatrices(Ndimensions,Nclusters*Ndimensions): set of Nclusters N-dimensional covariance matrices for the clusters
//      Nclusters: number of clusters considered
// 
// OUTPUT: 
//      void
//
void GaussianMixtureClusterer::computeGaussianMixtureModel(RefArrayXXd sample)
{
    // Compute the multivariate Gaussian for each cluster selected

    modelProbability.setZero();
    ArrayXXd argument(Npoints, Npoints);

    for (int i = 0; i < Nclusters; i++)
    {
        argument = differenceFromCenters.block(0, i*Npoints,Ndimensions, Npoints).matrix().transpose() * 
                            inverseOfCovarianceMatrices.block(0, i*Ndimensions,Ndimensions, Ndimensions).matrix() * 
                            differenceFromCenters.block(0, i*Npoints,Ndimensions, Npoints).matrix();
        multivariateGaussians.col(i) = (-0.5*Ndimensions*log(2.0*Functions::PI) -0.5*log(fabs(determinantOfCovarianceMatrices(i))) - 
                            0.5*argument).matrix().diagonal().array().exp();
    }

    // Compute total Gaussian Mixture Model

    modelProbability = multivariateGaussians.matrix()*amplitudes.matrix(); 


    // Compute assignment probabilities for each data point and the cluster responsibility

    assignmentProbabilities = (multivariateGaussians.rowwise()*amplitudes.transpose()).colwise()/modelProbability;
    responsibilities = assignmentProbabilities.colwise().sum();

}









// GaussianMixtureClusterer::updateClustersUntilConverged()
//
// PURPOSE: 
//      Evolve the cluster centers and covariance matrices according to the Gaussian Mixture Model. 
//      Given a set of cluster centers and covariance matrices, compute the corresponding GMM.
//      Then perform a weighted sum of the differences (between each point and the cluster centers) for each cluster, 
//      so that the cluster center and covariance matrix can be updated. 
//      The weights are in this case the assignment probabilities for each point to belong
//      to a given cluster, according to the definition of the GMM.
//
// INPUT:
//      sample(Ndimensions, Npoints):   sample of N-dimensional points
//      center(Ndimensions, Nclusters): set of N-dimensional coordinates of the cluster centers 
//      clusterSizes(Nclusters):        for each cluster, the number of points it contains
//      clusterIndices(Npoints):        for each point the index of the cluster it belongs to. This index
//                                      runs from 0 to Nclusters-1.
//      relTolerance:                   fitting the clusters converged when S < relTolerance, where S 
//                                      is the relative increment in total logarithm of GMM probability
// 
// OUTPUT:
//      True if the convergence of k-means was successful, false otherwise. A convergence is
//      considered successful if the relTolerance-criterion is satisfied (see above), _and_ if 
//      all clusters contain at least 2 points.


bool GaussianMixtureClusterer::updateClustersUntilConverged(RefArrayXXd sample)
{
    // Perform the GMM EM clustering iteration, each time improving the cluster centers and covariance matrices.
    // Once convergence is reached, determe which points belongs to which cluster
    
    double relativeIncrement;
    bool stopIterations = false;
    bool convergenceReached;
    unsigned int loop = 0;


    while (!stopIterations)
    {
        // Recompute/update the new cluster centers, which is simply
        // the barycenter of all points belonging to the cluster according to
        // the probabilities computed from each cluster (weights). 
    
        amplitudes = responsibilities/Npoints; 
        centers = sample.matrix() * assignmentProbabilities.matrix();
        centers = centers.rowwise()/responsibilities.transpose();


        // Update the covariance matrices for each cluster

        for (int i = 0; i < Nclusters; i++)
        {
            differenceFromCenters.block(0,i*Npoints,Ndimensions, Npoints) = sample.colwise() - centers.col(i);

            covarianceMatrices.block(0, i*Ndimensions,Ndimensions, Ndimensions) = 
                    differenceFromCenters.block(0, i*Npoints,Ndimensions, Npoints).matrix() * 
                    assignmentProbabilities.col(i).matrix().asDiagonal() * 
                    differenceFromCenters.block(0, i*Npoints,Ndimensions, Npoints).matrix().transpose();

            covarianceMatrices.block(0, i*Ndimensions,Ndimensions, Ndimensions) /= responsibilities(i);

            determinantOfCovarianceMatrices(i) = covarianceMatrices.block(0, i*Ndimensions,Ndimensions,Ndimensions).matrix().determinant();
           
            
            // If at least one of the covariance matrices determinant is zero, then the algorithm must be stopped and 
            // a new trial has to be done
            
            if (determinantOfCovarianceMatrices(i) == 0)
            {
                convergenceReached = false;
                return convergenceReached;
            }
            
            inverseOfCovarianceMatrices.block(0, i*Ndimensions,Ndimensions, Ndimensions).matrix() = 
            covarianceMatrices.block(0, i*Ndimensions,Ndimensions, Ndimensions).matrix().inverse();
        }
        
            
        // If this step has been reached, then covariance matrices are invertible and the GMM can be computed
            
        computeGaussianMixtureModel(sample);
       

        // A new set of clusters has been determined.
        // Decide whether the algorithm has converged. Convergence occurs when
        // the increment in the total model probability becomes lower than an input threshold.
        
        relativeIncrement = fabs(modelProbability.log().sum() - totalLogOfModelProbability)/fabs(minTotalLogOfModelProbability); 
        totalLogOfModelProbability = modelProbability.log().sum();

        if (relativeIncrement <= relTolerance)
        {
            stopIterations = true;
        }
       
        if (loop >= 30)                 // !!!!!!!!!!!!!
        {
            stopIterations = true;
        }

        loop++;

    }  // end GMM updating loop 
    
    
    // Convergence was properly reached, so return
    
    convergenceReached = true;
    
    return convergenceReached;  
}









// GaussianMixtureClusterer::obtainClusterMembership()
//
// PURPOSE: 
//      Given a sample of N-dimensional points, use the optimal solution from the Gaussian Mixture Model clustering algorithm 
// based on Expectation Maximization to determine the cluster membership of each sample point, as well as
// the number of points contained in each cluster.
//
// INPUT:
//      sample(Ndimensions, Npoints): sample of N-dimensional points
//      optimalClusterIndices(Npoints): for each point the index of the cluster it belongs to. This index
//                                      runs from 0 to Nclusters-1.
//      optimalClusterSizes(Nclusters): for each of the clusters, this vector contains the number of points
//      optimalAssignmentProbabilities(Npoints, Nclusters): assignment probabilities for each sample point according to each cluster
// 
// OUTPUT:
//      void
//

void GaussianMixtureClusterer::obtainClusterMembership(vector<int> &optimalClusterIndices, vector<int> &optimalClusterSizes,
                                                        RefArrayXXd optimalAssignmentProbabilities)
{
    Index maxIndexOfCluster;        // Store subscript containing index of cluster
    double maxAssignmentProbability;

    for (int i = 0; i < Npoints; ++i)
    {
        maxAssignmentProbability = optimalAssignmentProbabilities.row(i).maxCoeff(&maxIndexOfCluster);
        optimalClusterIndices[i] = maxIndexOfCluster;
        optimalClusterSizes[maxIndexOfCluster]++;
    }
}













// GaussianMixtureClusterer::cluster()
//
// PURPOSE: 
//      Given a sample of N-dimensional points, use the Gaussian Mixture Model clustering algorithm based on
//      Expectation Maximization to determine the optimal number of clusters that can be fitted to the data.
//
// INPUT:
//      sample(Ndimensions, Npoints): sample of N-dimensional points
//      optimalClusterIndices(Npoints): for each point the index of the cluster it belongs to. This index
//                                      runs from 0 to Nclusters-1.
//      optimalClusterSizes(Nclusters): for each of the clusters, this vector contains the number of points
// 
// OUTPUT:
//      The optimal number of clusters
//

int GaussianMixtureClusterer::cluster(RefArrayXXd sample, vector<int> &optimalClusterIndices, vector<int> &optimalClusterSizes)
{
    Npoints = sample.cols();
    ArrayXXd optimizedSample;

    // If activated, apply a dimensionality reduction of the data sample according to the chosen feature projector (e.g. PCA)

    if (featureProjectionActivated)
    {
        // If a dimensionality reduction is applied, from now on the clustering analysis 
        // will be performed on the sample of reduced dimensionality. 
        // The new sample will consist of a number of points equal to the original sample but with a number of
        // features (i.e. dimensions) <= than that of the original sample.
        
        optimizedSample = featureProjector.projection(sample);
    }
    else
    {
        optimizedSample = sample;
    }
    
    Ndimensions = optimizedSample.rows();
    modelProbability.resize(Npoints);

    bool convergedSuccessfully;
    unsigned int optimalNclusters;    


    // Best total model probability for a given number of clusters
    
    double bestTotalLogOfModelProbability;          
    
    
    // Best total model probability wrt number of clusters 
    
    double bestBICvalue = -1.0*numeric_limits<double>::max();
    double BICvalue;

    ArrayXXd bestCenters;
    ArrayXXd bestCovarianceMatrices;
    ArrayXXd bestAssignmentProbabilities;
    ArrayXXd optimalCenters;
    ArrayXXd optimalCovarianceMatrices;
    ArrayXXd optimalAssignmentProbabilities;


    // As we don't know a priori the optimal number of clusters, loop over a
    // user-specified range of clusters, and determine which number gives the
    // optimal clustering

    for (unsigned int c = minNclusters; c <= maxNclusters; ++c)
    {
        Nclusters = c;


        // Give an initial random set of amplitudes for the selected number of clusters
       
        uniform_real_distribution<> uniform01(0.0, 1.0);
        amplitudes.resize(Nclusters);
        
        for (unsigned int i = 0; i < Nclusters; i++)
        {
            amplitudes(i) = uniform01(engine);
        }

        amplitudes /= amplitudes.sum();

        
        // Resize relevant arrays according to the number of clusters

        differenceFromCenters.resize(Ndimensions, Npoints*Nclusters);
        differenceFromCenters.setZero();
        
        covarianceMatrices.resize(Ndimensions, Ndimensions*Nclusters);
        covarianceMatrices.setZero();

        centers.resize(Ndimensions, Nclusters);          // coordinates of each of the old cluster centers
        centers.setZero();
        
        inverseOfCovarianceMatrices.resize(Ndimensions, Ndimensions*Nclusters);
        inverseOfCovarianceMatrices.setZero();
        
        determinantOfCovarianceMatrices.resize(Nclusters);
        determinantOfCovarianceMatrices.setZero();

        multivariateGaussians.resize(Npoints, Nclusters);
        multivariateGaussians.setZero();

        assignmentProbabilities.resize(Npoints, Nclusters);
        assignmentProbabilities.setZero();

        responsibilities.resize(Nclusters);
        responsibilities.setZero();

        bestCenters.resize(Ndimensions, Nclusters);
        bestCovarianceMatrices.resize(Ndimensions, Ndimensions*Nclusters);
        bestAssignmentProbabilities.resize(Npoints, Nclusters);
       

        // The GMM EM algorithm is sensitive to the choice of the initial centers. 
        // We therefore run the algorithm 'Ntrial' times, and take the best clustering.
        
        bestTotalLogOfModelProbability = -1.0*numeric_limits<double>::max();
        
        for (int m = 0; m < Ntrials; ++m)
        {
            chooseInitialClusterCenters(optimizedSample);
            chooseInitialClusterCovarianceMatrices(optimizedSample);
            
            computeGaussianMixtureModel(optimizedSample);
            minTotalLogOfModelProbability = modelProbability.log().sum();
            totalLogOfModelProbability = minTotalLogOfModelProbability;

            convergedSuccessfully = updateClustersUntilConverged(optimizedSample);


            // If the convergence was not successfull (e.g. because some clusters contain 0 or 1 points),
            // we likely had an unfortunate set of initial cluster centers. In this case, simply continue
            // with the next 'trial'.
           

            if (!convergedSuccessfully) continue;
   

            // If we did obtain a successful convergence, compare it with the previous clusterings 
            // (all of them with the same number of clusters), and keep the best one.
            
            if (modelProbability.log().sum() > bestTotalLogOfModelProbability)
            {
                bestTotalLogOfModelProbability = modelProbability.log().sum();
                bestCenters = centers;
                bestCovarianceMatrices = covarianceMatrices;
                bestAssignmentProbabilities = assignmentProbabilities;
            }

        } // end loop over Ntrials to determine the best clustering trying different initial centers
       

        // Identify the number of clusters that is best at reproducing the sample. This can be done by directly inspecting
        // the total model probability since the GMM is a probability function by definition.
        // Note that this step is only necessary if the user selected more than one particular number of clusters.

        if (maxNclusters - minNclusters > 0)
        {
            // Evaluate the BIC value by taking into account the number of free parameters of the model
            // Considering a multivariate Gaussian in Ndimensions one has Ndimensions + Ndimensions*(Ndimensions + 1)/2
            // free parameters. This has to be multiplied by the number of mixture components in the model and add the
            // number of weights for each component as free parameters as well.

            BICvalue = bestTotalLogOfModelProbability - 0.5*log(Npoints)*(Nclusters-1 + 
                        Nclusters*(Ndimensions + Ndimensions*(Ndimensions+1)*0.5));
            
            if (BICvalue > bestBICvalue)
            {
                // We found a cluster combination that is better than anything found before. Save it.
                // In what follows 'best' refers to the best cluster configuration given a specific
                // number of clusters. 'optimal' refers to the optimal configuration over all possible
                // values for the number of clusters.
                
                bestBICvalue = BICvalue;
                
                optimalNclusters = Nclusters;
                optimalCenters.resize(Ndimensions, Nclusters);
                optimalCenters = bestCenters;
                optimalCovarianceMatrices.resize(Ndimensions, Ndimensions*Nclusters);
                optimalCovarianceMatrices = bestCovarianceMatrices;
                optimalAssignmentProbabilities.resize(Npoints, Nclusters);
                optimalAssignmentProbabilities = bestAssignmentProbabilities;
            }
        }
        else
        {
            // User allowed only 1 particular number of clusters. There is no need to check for the best
            // number of clusters.
                         
            optimalNclusters = Nclusters;
            optimalCenters.resize(Ndimensions, Nclusters);
            optimalCenters = bestCenters;
            optimalCovarianceMatrices.resize(Ndimensions, Ndimensions*Nclusters);
            optimalCovarianceMatrices = bestCovarianceMatrices;
            optimalAssignmentProbabilities.resize(Npoints, Nclusters);
            optimalAssignmentProbabilities = bestAssignmentProbabilities;
        }
    } // end loop over Nclusters
   
    // Finally, obtain cluster membership for each sample point with the best solution identified for this sample 
    // and range for number of clusters. In the case a feature projection was used, the cluster indices will stick in
    // an identical way from the reduced sample to the original sample because the number of data points is unchanged
    // during the dimensionality reduction process.

    optimalClusterSizes.resize(optimalNclusters);
    obtainClusterMembership(optimalClusterIndices, optimalClusterSizes, optimalAssignmentProbabilities);

    // That's it!

    return optimalNclusters;
}

