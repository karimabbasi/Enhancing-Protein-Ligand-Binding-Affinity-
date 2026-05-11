import numpy as np



def enrichment_factor(protein_index, predictions, true_label, x_percent):
  #######################################################
  # protein_index: (n,) array where n shows the number of pairs. each element shows the protein index of pair.
  # predictions: (n,) array which is the output of the model for pairs
  # true_label: (n,) array which determines the activ and inactive pairs (ground truth knowledge)
  # x_percent: shows the percent that enrichment factor should be computed
  ########################################################

  number_of_unique_proteins = np.max(protein_index)+1
  num_all_samples = predictions.shape[0]
  ntb_x = np.zeros((number_of_unique_proteins, 1))
  for i in range(number_of_unique_proteins):
    idx = np.argwhere(protein_index == i)
    predicted_values = predictions[idx]
    predicted_values = predicted_values.flatten()
    sorted_predicted_idx = np.argsort(predicted_values)[::-1]
    
    num_x_percent = np.min([np.floor(x_percent* num_all_samples), len(sorted_predicted_idx)])
    ntb_x_percent = 0
    for j in range(int(num_x_percent)):
      if (true_label[idx[sorted_predicted_idx[j]]] == 1):
        ntb_x_percent = ntb_x_percent + 1
    tmp = np.sum(true_label[idx])
    ntb_x[i] = ntb_x_percent / np.sum(true_label[idx])

  return np.mean(ntb_x)

 

# Example to show how to run the enrichment_factor
#protein_index = np.array([0, 0, 1, 1, 0, 0, 0, 1, 1 , 1])
#true_label = np.array([1, 1, 1, 0, 0, 1, 0, 1, 0, 1])
#predictions = np.array([0.3, 0.4, 0.9, 0.4, 0.1, 0.6, 0.9, 0.3, 0.5, 0.7])
#x_percent = 0.1 
#print(enrichment_factor(protein_index, predictions, true_label, x_percent))