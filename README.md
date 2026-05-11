# Enhancing Protein-Ligand Binding Affinity Prediction via Hard Pair Retrieval and Semi-Supervised Graph Convolutional Networks
# Enhancing Protein-Ligand 
########################################################## 

The source code for
Enhancing Protein-Ligand Binding Affinity Prediction via Hard Pair Retrieval and Semi-Supervised Graph Convolutional Networks

##########################################################

**Requirements**
you can run the code under the following environment: Keras 2.14.0, Tensorflow 2.14.0, and Python 3.10.12.

##########################################################

**Data**


Download the data from the following link.

https://drive.google.com/open?id=1B72WnWMbywxK2M9RntquRWQ3cm6U9YoW

Download the folded data from the following link:

https://drive.google.com/open?id=15KotSJWknMOAnHM68RpOh_rqMISsMwsE
##########################################################

**Usage**

At first, Keras-gcn should be installed. The hyperparameters are set in the config.py file. You can select Davis, KIBA, PDBind, and BindingDB datasets in this file. Also, The other hyperparameters could be set to your desired values. Then, you can run the sgcn.py.
##########################################################

**Hint**

we have modified the ‘kegra’ package, and we have uploaded this package to our repository. This package has been downloaded from this link, and if you use our code, you should cite our paper and their paper. You should cite their paper as follows:
@inproceedings{kipf2017semi,
  title={Semi-Supervised Classification with Graph Convolutional Networks},
  author={Kipf, Thomas N. and Welling, Max},
  booktitle={International Conference on Learning Representations (ICLR)},
  year={2017}
}}

