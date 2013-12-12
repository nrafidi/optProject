Smoothing Proximal Gradient Method for Multi-task (response) Graph-guided Fused Lasso

Please run install_mex.m to mex .c files (I have mexed the .c files under windows for both 32 and 64 bit machines)

Run test.m for testing multi-task graph-guided fused lasso

The files for main optimization procedure include

(1) SPG_multi.m  Smoothing Proximal Gradient (SPG) with a pre-computed Lipschitz constant, suitable for medium-scale problem

(2) SPG_multi_linesearch.m  Smoothing Proximal Gradient (SPG) with a line-search scheme on Lipschitz constant, suitable for large-scale problem

Please cite if you use the code in your research. 

@INPROCEEDINGS{xichen:UAI:11,
  author = {Xi Chen and Qihang Lin and Seyoung Kim and Jaime Carbonell and Eric
    Xing},
  title = {Smoothing Proximal Gradient Method for General Structured Sparse
    Learning},
  booktitle = {Proceedings of Uncertainty in Artificial Intelligence (UAI)},
  year = {2011}
}

If there is any comment or question, please contact: xichen@cs.cmu.edu