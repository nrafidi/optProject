Code for running smoothing proximal gradient method for graph-structured input-output lasso, with fused and group penalties

Run example.m to see an example with synthetic data

Usage is Bhat = graph_siol_spg(X,Y,Psi,Theta,lambda,gamma1,gamma2,penFlag)

Inputs:
   X: n-by-p matrix of input features, one feature per column [required]
   Y: n-by-q matrix of output tasks, one task per column [required]
   Psi: p-by-p matrix encoding graph structure among inputs [optional; default is to use sample covariance matrix]
   Theta: q-by-q matrix encoding graph structure among outputs [optional; default is to use sample covariance matrix]
   lambda: weight for L1 norm penalty [optional; default is 1]
   gamma1: weight for input structured sparsity penalty [optional; default is 1]
   gamma2: weight for output structured sparsity penalty [optional; default is 1]
   penFlag: indicates whether to use fused lasso penalty ('fused') or 
       group lasso penalty ('group') [optional; default is 'fused']

To use the default value for any argument, pass in [] or nothing, e.g.
   Bhat = graph_siol_spg(X,Y,[],[],lambda,gamma1,gamma2,penFlag)
   Bhat = graph_siol_spg(X,Y)