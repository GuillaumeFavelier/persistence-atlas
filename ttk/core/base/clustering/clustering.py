def compute_bic(centroids, labels, n_clusters, X):
    """
    Computes the BIC metric for a given clusters
    Parameters:
    -----------------------------------------
    kmeans:  List of clustering object from scikit learn
    X     :  multidimension np array of data points
    Returns:
    -----------------------------------------
    BIC value
    """
    import numpy as np
    from scipy.spatial import distance
    # assign centers and labels
    centers = centroids
    # number of clusters
    m = n_clusters
    # size of the clusters
    n = np.bincount(labels)
    # size of data set
    N, d = X.shape

    # compute variance for all clusters beforehand
    cl_var = (1.0 / (N - m) / d) * sum([sum(distance.cdist(X[np.where(labels == i)], [centers[i]],
                                       'euclidean')**2) for i in range(m)])

    const_term = 0.5 * m * np.log(N) * (d + 1)

    BIC = (np.sum([n[i] * np.log(n[i]) -
                n[i] * np.log(N) -
                ((n[i] * d) / 2) * np.log(2 * np.pi * cl_var) -
                ((n[i] - 1) * d / 2) for i in range(m)]) - const_term)

    return BIC

def getValueKey(item):
    return item[1]

def getIdKey(item):
    return item[0]

def doIt(M, values, nmin, nmax, cgap):
    # for debug purpose:
    # print("[Clustering] Python: Input maps:")
    # print(M)
    import importlib
    import math
    import numpy as np
    # import matplotlib.pyplot as plt
    from sklearn.cluster import k_means

    print("[Clustering] Python: Input eigenvalues:")
    print(values)

    print("[Clustering] min: ", nmin)
    print("[Clustering] max: ", nmax)

    kmin = nmin
    kmax = nmax + 1

    k=kmin

    if(kmin!=nmax):
        print("[Clustering] nmin!=nmax")
        maxgap = -1000000
        maxid = 0

        # check if Numpy is installed
        loader = importlib.find_loader('numpy')
        found = loader is not None
        if found:
            print("[Clustering] Python: numpy module found.")
        else:
            print("[Clustering] Python error: numpy module not found.")
            return 0

        # check if scipy is installed
        loader = importlib.find_loader('scipy')
        found = loader is not None
        if found:
            print("[Clustering] Python: scipy module found.")
        else:
            print("[Clustering] Python error: scipy module not found.")
            return 0

        # check if Scikit-learn is installed
        loader = importlib.find_loader('sklearn')
        found = loader is not None
        if found:
            print("[Clustering] Python: sklearn module found.")
        else:
            print("[Clustering] Python error: sklearn module not found.")
            return 0

        from sklearn.cluster import k_means

        eigenvalues = sorted(values[0], reverse=True)
        #x=[]
        #for i in range(0, len(eigenvalues)):
        #    print("value: ", eigenvalues[i])
        #    print("Id: ", i)
        #    x.append(i+1)
        #plt.plot(x, eigenvalues)
        #plt.show(x)

        gapsPairs=[]
        gaps=[]
        x=[]

        if(len(eigenvalues) < kmax):
            kmax = len(eigenvalues)

        averageGap = 0.
        for i in range(2, kmax):
            gap = math.fabs(eigenvalues[i]-eigenvalues[i-1])
            gaps.append(gap)
            x.append(i)
            gapsPairs.append([i, gap])
            #print("Gap: ", gap)
            #print("Id: ", i)
            averageGap = averageGap + gap

        # export to CSV
        # import csv
        # csvfile = open('/tmp/clustering_gaps.csv', 'w')
        # csvwriter = csv.writer(csvfile, delimiter=';')
        # csvwriter.writerow(['eigengaps'])
        # for i in range(len(gaps)):
        #     csvwriter.writerow([gaps[i]])

        averageGap = averageGap/kmax


        sortedGapsPairs=sorted(gapsPairs, key=getValueKey, reverse=True)
        print("[Clustering] Python Sorted: ", sortedGapsPairs)

        maxgap=sortedGapsPairs[0][1]


        significantGapNb = 1
        for i in range(1, len(sortedGapsPairs)):
            currentGap=sortedGapsPairs[i][1]
            #print("percentage: ", currentGap/maxgap)
            if currentGap > 0.00001 and (currentGap/maxgap > 0.1):
                significantGap = i

        if(cgap>significantGap):
            significantGap = cgap
        significantGapsPairs = sortedGapsPairs[0:significantGap+1]
        #print("Sorted cut significant: ", significantGapsPairs)

        significantGapsPairs = sorted(significantGapsPairs, key=getIdKey)
        print("[Clustering] Python Significant Pairs: ", significantGapsPairs)

        # export to CSV
        # import csv
        # csvfile = open('/tmp/clustering_sortedGaps.csv', 'w')
        # csvwriter = csv.writer(csvfile, delimiter=';')
        # csvwriter.writerow(['k', 'gaps'])
        # for i in range(len(significantGapsPairs)):
        #     csvwriter.writerow([significantGapsPairs[i][0], significantGapsPairs[i][1]])

        #significantGapsPairs=sorted(significantGapsPairs, key=getIdKey)
        #print("Significant: ", significantGapsPairs)
        #print("Current gap: ", cgap)

        maxgap = significantGapsPairs[cgap-1][1]
        maxid = significantGapsPairs[cgap-1][0]

        print("[Clustering] Python maxgap: ", maxgap)
        print("[Clustering] Python maxid: ", maxid)
        print("[Clustering] Python Percentage average over max gap: ", averageGap/maxgap)

        k = maxid

        # if(k<0):
        #     k=nmax

        # if( averageGap/maxgap > 0.6 ):
        #     k = nmax
        #     print("[Clustering] Python not significant gap found in range: ", averageGap/maxgap)

        print("[Clustering] Python chosen k: ", k)

        # plt.plot(x, gaps)
        # plt.show(x)

    _, labels, _ = k_means(M, n_clusters=k, random_state=0)

    #for i in range(0, len(labels)):
    #    if(labels[i] < 1):
    #        labels[i]=1

    result = np.asarray(labels)

    # for debug purpose:
    #print("[Clustering] Python: Output assignation:")
    #print(result)

    return result
