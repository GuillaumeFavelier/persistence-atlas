def _graph_connected_component(graph, node_id):
    """Find the largest graph connected components that contains one
    given node
    Parameters
    ----------
    graph : array-like, shape: (n_samples, n_samples)
        adjacency matrix of the graph, non-zero weight means an edge
        between the nodes
    node_id : int
        The index of the query node of the graph
    Returns
    -------
    connected_components_matrix : array-like, shape: (n_samples,)
        An array of bool value indicating the indexes of the nodes
        belonging to the largest connected components of the given query
        node
    """
    import numpy as np
    from scipy import sparse

    n_node = graph.shape[0]
    if sparse.issparse(graph):
        # speed up row-wise access to boolean connection mask
        graph = graph.tocsr()
    connected_nodes = np.zeros(n_node, dtype=np.bool)
    nodes_to_explore = np.zeros(n_node, dtype=np.bool)
    nodes_to_explore[node_id] = True
    for _ in range(n_node):
        last_num_component = connected_nodes.sum()
        np.logical_or(connected_nodes, nodes_to_explore, out=connected_nodes)
        if last_num_component >= connected_nodes.sum():
            break
        indices = np.where(nodes_to_explore)[0]
        nodes_to_explore.fill(False)
        for i in indices:
            if sparse.issparse(graph):
                neighbors = graph[i].toarray().ravel()
            else:
                neighbors = graph[i]
            np.logical_or(nodes_to_explore, neighbors, out=nodes_to_explore)
    return connected_nodes

def _graph_is_connected(graph):
    """ Return whether the graph is connected (True) or Not (False)
    Parameters
    ----------
    graph : array-like or sparse matrix, shape: (n_samples, n_samples)
        adjacency matrix of the graph, non-zero weight means an edge
        between the nodes
    Returns
    -------
    is_connected : bool
        True means the graph is fully connected and False means not
    """
    from scipy import sparse
    from scipy.sparse.csgraph import connected_components

    if sparse.isspmatrix(graph):
        # sparse graph, find all the connected components
        n_connected_components, _ = connected_components(graph)
        return n_connected_components == 1
    else:
        # dense graph, find all connected components start from node 0
        return _graph_connected_component(graph, 0).sum() == graph.shape[0]

def _set_diag(laplacian, value, norm_laplacian):
    """Set the diagonal of the laplacian matrix and convert it to a
    sparse format well suited for eigenvalue decomposition
    Parameters
    ----------
    laplacian : array or sparse matrix
        The graph laplacian
    value : float
        The value of the diagonal
    norm_laplacian : bool
        Whether the value of the diagonal should be changed or not
    Returns
    -------
    laplacian : array or sparse matrix
        An array of matrix in a form that is well suited to fast
        eigenvalue decomposition, depending on the band width of the
        matrix.
    """
    import numpy as np
    from scipy import sparse

    n_nodes = laplacian.shape[0]
    # We need all entries in the diagonal to values
    if not sparse.isspmatrix(laplacian):
        if norm_laplacian:
            laplacian.flat[::n_nodes + 1] = value
    else:
        laplacian = laplacian.tocoo()
        if norm_laplacian:
            diag_idx = (laplacian.row == laplacian.col)
            laplacian.data[diag_idx] = value
        # If the matrix has a small number of diagonals (as in the
        # case of structured matrices coming from images), the
        # dia format might be best suited for matvec products:
        n_diags = np.unique(laplacian.row - laplacian.col).size
        if n_diags <= 7:
            # 3 or less outer diagonals on each side
            laplacian = laplacian.todia()
        else:
            # csr has the fastest matvec and is thus best suited to
            # arpack
            laplacian = laplacian.tocsr()
    return laplacian

def my_spectral_embedding(adjacency, n_components=8, eigen_solver=None,
                       random_state=None, eigen_tol=0.0,
                       norm_laplacian=False, drop_first=True):
    """Project the sample on the first eigenvectors of the graph Laplacian.
    The adjacency matrix is used to compute a normalized graph Laplacian
    whose spectrum (especially the eigenvectors associated to the
    smallest eigenvalues) has an interpretation in terms of minimal
    number of cuts necessary to split the graph into comparably sized
    components.
    This embedding can also 'work' even if the ``adjacency`` variable is
    not strictly the adjacency matrix of a graph but more generally
    an affinity or similarity matrix between samples (for instance the
    heat kernel of a euclidean distance matrix or a k-NN matrix).
    However care must taken to always make the affinity matrix symmetric
    so that the eigenvector decomposition works as expected.
    Note : Laplacian Eigenmaps is the actual algorithm implemented here.
    Read more in the :ref:`User Guide <spectral_embedding>`.
    Parameters
    ----------
    adjacency : array-like or sparse matrix, shape: (n_samples, n_samples)
        The adjacency matrix of the graph to embed.
    n_components : integer, optional, default 8
        The dimension of the projection subspace.
    eigen_solver : {None, 'arpack', 'lobpcg', or 'amg'}, default None
        The eigenvalue decomposition strategy to use. AMG requires pyamg
        to be installed. It can be faster on very large, sparse problems,
        but may also lead to instabilities.
    random_state : int, RandomState instance or None, optional, default: None
        A pseudo random number generator used for the initialization of the
        lobpcg eigenvectors decomposition.  If int, random_state is the seed
        used by the random number generator; If RandomState instance,
        random_state is the random number generator; If None, the random number
        generator is the RandomState instance used by `np.random`. Used when
        ``solver`` == 'amg'.
    eigen_tol : float, optional, default=0.0
        Stopping criterion for eigendecomposition of the Laplacian matrix
        when using arpack eigen_solver.
    norm_laplacian : bool, optional, default=True
        If True, then compute normalized Laplacian.
    drop_first : bool, optional, default=True
        Whether to drop the first eigenvector. For spectral embedding, this
        should be True as the first eigenvector should be constant vector for
        connected graph, but for spectral clustering, this should be kept as
        False to retain the first eigenvector.
    Returns
    -------
    embedding : array, shape=(n_samples, n_components)
        The reduced samples.
    Notes
    -----
    Spectral Embedding (Laplacian Eigenmaps) is most useful when the graph
    has one connected component. If there graph has many components, the first
    few eigenvectors will simply uncover the connected components of the graph.
    References
    ----------
    * https://en.wikipedia.org/wiki/LOBPCG
    * Toward the Optimal Preconditioned Eigensolver: Locally Optimal
      Block Preconditioned Conjugate Gradient Method
      Andrew V. Knyazev
      http://dx.doi.org/10.1137%2FS1064827500366124
    """
    import warnings

    import numpy as np
    from scipy import sparse
    from scipy.linalg import eigh
    from scipy.sparse.linalg import eigsh, lobpcg

    from sklearn.base import BaseEstimator
    from sklearn.externals import six
    from sklearn.utils import check_random_state, check_array, check_symmetric
    from sklearn.utils.extmath import _deterministic_vector_sign_flip
    from sklearn.metrics.pairwise import rbf_kernel
    from sklearn.neighbors import kneighbors_graph

    adjacency = check_symmetric(adjacency)
    try:
        from pyamg import smoothed_aggregation_solver
    except ImportError:
        if eigen_solver == "amg":
            raise ValueError("The eigen_solver was set to 'amg', but pyamg is "
                             "not available.")
    if eigen_solver is None:
        eigen_solver = 'arpack'
    elif eigen_solver not in ('arpack', 'lobpcg', 'amg'):
        raise ValueError("Unknown value for eigen_solver: '%s'."
                         "Should be 'amg', 'arpack', or 'lobpcg'"
                         % eigen_solver)
    random_state = check_random_state(random_state)
    n_nodes = adjacency.shape[0]
    # Whether to drop the first eigenvector
    if drop_first:
        n_components = n_components + 1
    if not _graph_is_connected(adjacency):
        warnings.warn("Graph is not fully connected, spectral embedding"
                      " may not work as expected.")
    laplacian, dd = sparse.csgraph.laplacian(adjacency, normed=norm_laplacian,
                                             return_diag=True)
    if (eigen_solver == 'arpack' or eigen_solver != 'lobpcg' and
       (not sparse.isspmatrix(laplacian) or n_nodes < 5 * n_components)):
        # lobpcg used with eigen_solver='amg' has bugs for low number of nodes
        # for details see the source code in scipy:
        # https://github.com/scipy/scipy/blob/v0.11.0/scipy/sparse/linalg/eigen
        # /lobpcg/lobpcg.py#L237
        # or matlab:
        # http://www.mathworks.com/matlabcentral/fileexchange/48-lobpcg-m
        laplacian = _set_diag(laplacian, 1, norm_laplacian)

        # Here we'll use shift-invert mode for fast eigenvalues
        # (see http://docs.scipy.org/doc/scipy/reference/tutorial/arpack.html
        #  for a short explanation of what this means)
        # Because the normalized Laplacian has eigenvalues between 0 and 2,
        # I - L has eigenvalues between -1 and 1.  ARPACK is most efficient
        # when finding eigenvalues of largest magnitude (keyword which='LM')
        # and when these eigenvalues are very large compared to the rest.
        # For very large, very sparse graphs, I - L can have many, many
        # eigenvalues very near 1.0.  This leads to slow convergence.  So
        # instead, we'll use ARPACK's shift-invert mode, asking for the
        # eigenvalues near 1.0.  This effectively spreads-out the spectrum
        # near 1.0 and leads to much faster convergence: potentially an
        # orders-of-magnitude speedup over simply using keyword which='LA'
        # in standard mode.
        try:
            # We are computing the opposite of the laplacian inplace so as
            # to spare a memory allocation of a possibly very large array
            laplacian *= -1
            v0 = random_state.uniform(-1, 1, laplacian.shape[0])
            lambdas, diffusion_map = eigsh(laplacian, k=n_components,
                                           sigma=1.0, which='LM',
                                           tol=eigen_tol, v0=v0)
            embedding = diffusion_map.T[n_components::-1] * dd
        except RuntimeError:
            # When submatrices are exactly singular, an LU decomposition
            # in arpack fails. We fallback to lobpcg
            eigen_solver = "lobpcg"
            # Revert the laplacian to its opposite to have lobpcg work
            laplacian *= -1
    if eigen_solver == 'amg':
        # Use AMG to get a preconditioner and speed up the eigenvalue
        # problem.
        if not sparse.issparse(laplacian):
            warnings.warn("AMG works better for sparse matrices")
        # lobpcg needs double precision floats
        laplacian = check_array(laplacian, dtype=np.float64,
                                accept_sparse=True)
        laplacian = _set_diag(laplacian, 1, norm_laplacian)
        ml = smoothed_aggregation_solver(check_array(laplacian, 'csr'))
        M = ml.aspreconditioner()
        X = random_state.rand(laplacian.shape[0], n_components + 1)
        X[:, 0] = dd.ravel()
        lambdas, diffusion_map = lobpcg(laplacian, X, M=M, tol=1.e-12,
                                        largest=False)
        embedding = diffusion_map.T * dd
        if embedding.shape[0] == 1:
            raise ValueError

    elif eigen_solver == "lobpcg":
        # lobpcg needs double precision floats
        laplacian = check_array(laplacian, dtype=np.float64,
                                accept_sparse=True)
        if n_nodes < 5 * n_components + 1:
            # see note above under arpack why lobpcg has problems with small
            # number of nodes
            # lobpcg will fallback to eigh, so we short circuit it
            if sparse.isspmatrix(laplacian):
                laplacian = laplacian.toarray()
            lambdas, diffusion_map = eigh(laplacian)
            embedding = diffusion_map.T[:n_components] * dd
        else:
            laplacian = _set_diag(laplacian, 1, norm_laplacian)
            # We increase the number of eigenvectors requested, as lobpcg
            # doesn't behave well in low dimension
            X = random_state.rand(laplacian.shape[0], n_components + 1)
            X[:, 0] = dd.ravel()
            lambdas, diffusion_map = lobpcg(laplacian, X, tol=1e-15,
                                            largest=False, maxiter=2000)
            embedding = diffusion_map.T[:n_components] * dd
            if embedding.shape[0] == 1:
                raise ValueError
    embedding = _deterministic_vector_sign_flip(embedding)
    if drop_first:
        vectors = embedding[1:n_components].T
    else:
        vectors = embedding[:n_components].T

    return (lambdas, vectors)


def computeStress(k, D, M, C):
    import numpy as np
    from scipy.spatial import distance

    # normalize D
    # d = distance.squareform(D)
    d = np.array(D, copy=True)
    dmin = np.min(d)
    d = np.add(d, -dmin)
    dmax = np.max(d)
    d = np.divide(d, dmax)

    # normalize M
    f = distance.pdist(M)
    f = distance.squareform(f)
    fmin = np.min(f)
    f = np.add(f, -fmin)
    fmax = np.max(f)
    f = np.divide(f, fmax)

    # prepare final normalization
    n = np.sum(np.power(d, 2))

    # follow connectivity rules
    n = d.shape[0]
    for i in range(n):
        for j in range(n):
            if C[i, j] == 0:
                d[i, j] = 1.0
                f[i, j] = 1.0

    e = np.sum(np.power(f - d, 2)) / n
    s = np.sqrt(e)

    return (k, s)


def parallelStress(r, D, M, C, q):
    result = list()
    for k in r:
        m = M[:, 0:k]
        ret = computeStress(k, D, m, C)
        result.append(ret)

    q.put(result)


def parallelEmbedding(kmin, kmax, B, D, C, njobs):
    import numpy as np
    # import matplotlib.pyplot as plt
    from multiprocessing import Process, Queue
    l, M = my_spectral_embedding(B, n_components=kmax, random_state=0)

    # create processes container and shared structure
    P = list()
    q = Queue()
    R = np.array_split(range(kmin, kmax + 1), njobs)
    for r in R:
        p = Process(target=parallelStress, args=(r, D, M, C, q))
        p.start()
        P.append(p)

    result = []
    for p in P:
        while p.is_alive():
            p.join(timeout=1)
            while not q.empty():
                ret = q.get(block=False)
                result += ret

    dimensionality, stress = zip(*result)

    x = stress.index(min(stress))
    x = dimensionality[x]

    dimensionality, stress = zip(*sorted(zip(dimensionality, stress)))

    # export to CSV
    # import csv
    # csvfile = open('/tmp/spectralEmbedding.csv', 'w')
    # csvwriter = csv.writer(csvfile, delimiter=';')
    # csvwriter.writerow(['dimensionality', 'stress'])
    # for i in range(1, len(dimensionality)):
    #     csvwriter.writerow([dimensionality[i], stress[i]])

    # plt.plot(dimensionality, stress)
    # plt.show()

    return (l, M[:, 0:x], M)


def sequentialEmbedding(kmin, kmax, B, D, C):
    # import matplotlib.pyplot as plt
    l, M = my_spectral_embedding(B, n_components=kmax, random_state=0)

    x = 1
    minStress = -1
    stress = list()
    dimensionality = list()
    for k in range(kmin, kmax + 1):
        m = M[:, 0:k]
        ret = computeStress(k, D, m, C)
        # get stress minimizer
        if minStress == -1:
            minStress = ret[1]
            x = ret[0]
        elif ret[1] < minStress:
            minStress = ret[1]
            x = ret[0]
        # append for plot
        dimensionality.append(ret[0])
        stress.append(ret[1])

    # export to CSV
    # import csv
    # csvfile = open('/tmp/spectralEmbedding.csv', 'w')
    # csvwriter = csv.writer(csvfile, delimiter=';')
    # csvwriter.writerow(['dimensionality', 'stress'])
    # for i in range(1, len(dimensionality)):
    #     csvwriter.writerow([dimensionality[i], stress[i]])

    # plt.plot(dimensionality, stress)
    # plt.show()

    # compute spectral clustering
    return (l, M[:, 0:x], M)

def getValueKey(item):
    return item[1]

def getIdKey(item):
    return item[0]

def sequentialGapEmbedding(kmin, kmax, B, D, C):
    import importlib
    import math

    l, M = my_spectral_embedding(B, n_components=kmax, random_state=0)

    dimension = kmax
    gapsPairs=[]

    eigenvalues = sorted(l, reverse=True)

    averageGap = 0.
    for i in range(2, len(eigenvalues)):
        gap = math.fabs(eigenvalues[i]-eigenvalues[i-1])
        gapsPairs.append([i, gap])
        averageGap = averageGap + gap


    averageGap = averageGap/len(eigenvalues)
    sortedGapsPairs=sorted(gapsPairs, key=getValueKey, reverse=True)

    maxgap=sortedGapsPairs[0][1]

    significantGapNb = 1
    for i in range(1, len(sortedGapsPairs)):
        currentGap=sortedGapsPairs[i][1]

        if maxgap > 0.00001 and (currentGap/maxgap > 0.1):
            significantGap = i

    significantGapsPairs = sortedGapsPairs[0:significantGap+1]

    significantGapsPairs = sorted(significantGapsPairs, key=getIdKey)

    dimension=significantGapsPairs[0][0] - 1

    print("[SpectralEmbedding] Python: Found dimension from gap: ", dimension)
    return (l, M[:, 0:dimension], M)

def doIt(D, mincomponents, maxcomponents, nneighbors, sigma, njobs):
    # for debug purpose:
    print("[SpectralEmbedding] Python: Input distance matrix:")
    print(D)

    import importlib

    # check if Numpy is installed
    loader = importlib.find_loader('numpy')
    found = loader is not None
    if found:
        print("[SpectralEmbedding] Python: numpy module found.")
    else:
        print("[SpectralEmbedding] Python error: numpy module not found.")
        return 0

    # check if scipy is installed
    loader = importlib.find_loader('scipy')
    found = loader is not None
    if found:
        print("[SpectralEmbedding] Python: scipy module found.")
    else:
        print("[SpectralEmbedding] Python error: scipy module not found.")
        return 0

    # check if Scikit-learn is installed
    loader = importlib.find_loader('sklearn')
    found = loader is not None
    if found:
        print("[SpectralEmbedding] Python: sklearn module found.")
    else:
        print("[SpectralEmbedding] Python error: sklearn module not found.")
        return 0

    from sklearn.neighbors import kneighbors_graph

    forcedDimension = sigma
    print("[SpectralEmbedding] Python: maxcomponents: ", maxcomponents)

    connectivity = kneighbors_graph(D, nneighbors, include_self=True)
    B = 0.5 * (connectivity + connectivity.T)
    print("Weight matrix:")
    print(B.toarray())

    maxcomponents = max(maxcomponents, 3)

    if forcedDimension > 0:
        values, M = my_spectral_embedding(B, n_components=maxcomponents, random_state=0)
        vectors = M[:, 0:forcedDimension]
    else:

        values, vectors, M = sequentialGapEmbedding(mincomponents, maxcomponents, B, D, connectivity)

    coords = M[:, 0:3]

    # for debug purpose:
    print("[SpectralEmbedding] Python: Output eigenvectors:")
    print(vectors)
    print("[SpectralEmbedding] Python: Output eigenvalues:")
    print(values)
    print("[SpectralEmbedding] Python: Output coords:")
    print(coords)

    # padLength = len(D) - len(values)
    # Z = np.pad(values, (0, padLength), mode='constant')
    # Z = np.expand_dims(Z, axis=0)
    # result = np.concatenate((vectors, Z.T), axis=1)

    L = list()
    L.append(values)
    L.append(vectors)
    L.append(coords)

    return L
