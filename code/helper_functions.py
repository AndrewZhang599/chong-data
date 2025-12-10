import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os


def quantile_list(data,nquantiles):
  """Function that finds the value of n quantiles in an array. \n
  data: 1-D array \n
  nquantiles: number of quantiles to return"""
  quantiles_list = []
  for i in range(nquantiles):
    q = i/nquantiles
    if len(data) == 0:
      data = [0,0]
    quantiles_list.append(np.quantile(data,q))
  return np.array(quantiles_list)

def quantile_data(data,nquantiles):
    """Function that returns a list of the data split into quantiles. \n
  data: 1-D array \n
  nquantiles: number of quantiles to split the data into"""
    data_list = []
    l = np.linspace(0,1,nquantiles+1)
    l = np.quantile(data.rxtime,l)
    for i in range(nquantiles):
        #print((data.rt>=l[i])*(data.rt<l[i+1]))
        data_list.append(data[(data.rxtime>=l[i])*(data.rxtime<l[i+1])])
    return data_list

def rate_by_quantile_no_subs(data,response,nquants):
  """Returns the rate of occurence for a certain response within each n quantile(s) \n
  data: dataframe of hddm-ready data \n
  response: The specifc response to provide the rate of occurence for \n
  nquants: number of quantiles to return rate for"""
  out = np.zeros(nquants)
  qd = quantile_data(data,nquants)
  for j in range(nquants):
    out[j] = np.mean(qd[j].response==response)
  return out


def rate_by_quantile_no_subs_mult_choices(data,responses,nquants):
  """Returns the rate of occurence for a certain response within each n quantile(s) \n
  data: dataframe of hddm-ready data \n
  response: a list of The specifc response to provide the rate of occurence for \n
  nquants: number of quantiles to return rate for"""
  out = np.zeros(nquants)
  qd = quantile_data(data,nquants)
  total = np.zeros((len(responses)))
  for j in range(nquants):
    for i in range(len(responses)):
      total[i] = np.sum(qd[j].response == responses[i])
      # total[i,1] = len(qd[j].response == responses[i])
    out[j] = np.sum(total)/len(qd[0])
    print(out[j])
  return out


def rate_by_quantile(data,response,nquants):
  """Returns the rate of occurence for a certain response within each n quantile(s) \n
  data: dataframe of hddm-ready data \n
  response: The specifc response to provide the rate of occurence for \n
  nquants: number of quantiles to return rate for"""

  subs = np.unique(data.subj)
  out = np.zeros((len(subs),nquants))
  for i in range(len(subs)):
      d = data[data.subj==subs[i]]
      qd = quantile_data(d,nquants)
      for j in range(nquants):
          out[i,j] = np.mean(qd[j].response==response)
          print(out[i,j])
  print(np.mean(out,axis = 0))
  return np.mean(out,axis=0)

def load_posterior_predictive(model, task = 0, coherence = 0, group = 0):
  """Loads the posterior predictive of a model with the specified paramters"""
  files = os.listdir("data/posterior_predictives/{}".format(model))
  isempty = True
  for file in files:
    if (model in file) and ("task_{}".format(task) in file) and ("coh_{}".format(coherence) in file) and ("group_{}".format(group) in file):
      data = pd.read_pickle("data/posterior_predictives/{}/{}".format(model,file))
      isempty = False
  if isempty:
    print("model with these arguments does not exist")
  else:
    return data


### Plotting functions ###

def err_plot(data,nquants, color = 'red'):
  """Plots the rate of high and low dimension errors in each quantile"""
  plt.plot(rate_by_quantile(data, 1, nquants) + rate_by_quantile(data,0, nquants) ,'-^', color = color)
  plt.plot(rate_by_quantile(data, 2, nquants) + rate_by_quantile(data,0, nquants),':v', color = color)
  plt.ylim((0, .8))

def err_plot_mean(dataset_list, nquants, color):
  """Goes thru each dataset in thte list. Means the rate. Then plots the rate of high and low dimension errors in each quantile"""
  high_dim_err = np.zeros((len(dataset_list), nquants))
  low_dim_err = np.zeros((len(dataset_list), nquants))
  for d in range(len(dataset_list)):
    high_dim_err[d,:] = rate_by_quantile(dataset_list[d], 1, nquants) + rate_by_quantile(dataset_list[d],0, nquants) 
    low_dim_err[d,:] = rate_by_quantile(dataset_list[d], 2, nquants) + rate_by_quantile(dataset_list[d],0, nquants) 
  plt.plot(np.mean(high_dim_err, axis = 0),'-^', color = color)
  plt.plot(np.mean(low_dim_err, axis = 0),':v', color = color)
  plt.ylim((0, .8))
  
def err_plot_new(data,nquants):
  plt.plot(rate_by_quantile(data,1,nquants) + rate_by_quantile(data,0,nquants),'-o')
def single_qpp(data_table, nquants, color, coh_level='highDim'):
  """Plots a single quantile probability plot. If coherence level is 'highDim', 
  then errors in the high dimension are used for the error rates. Otherwise, low dimension 
  errors are used."""
  if coh_level == 'highDim':
    for i in [-1,1]:

        data = data_table[data_table.highDimCoh == i]
        rate = np.mean(data.response==1)
        if i == -1:
          line_style = 'v'
        else:
          line_style = '^'
        plt.plot(np.repeat(rate,nquants),quantile_list(data.rxtime[data.response==1],nquants),':' + line_style, color = color)
        rate = sum(data.response==3)/len(data.rxtime)
        plt.plot(np.repeat(rate,nquants),quantile_list(data.rxtime[data.response==3],nquants),'-' + line_style, color = color)

  else:
    for i in [-1,1]:
        data = data_table[data_table.lowDimCoh == i]
        rate = sum(data.response==2)/len(data.rxtime)
        if i == -1:
          line_style = 'v'
        else:
          line_style = '^'
        plt.plot(np.repeat(rate,nquants),quantile_list(data.rxtime[data.response==2],nquants),':' + line_style, color = color)

        rate = sum(data.response==3)/len(data.rxtime)
        plt.plot(np.repeat(rate,nquants),quantile_list(data.rxtime[data.response==3],nquants),'-' + line_style, color = color)

def single_qpp_mean(dataset_list, nquants, color, coh_level = 'highDim'):
  qlist_array = np.zeros((len(dataset_list), 2,2, nquants))      #dataset, 2 #coherence = [-1,1],2 #incorrect, correct, nquants
  rate_array = np.zeros((len(dataset_list), 2,2))
  coherence_list = [-1,1]
  for d in range(len(dataset_list)):
    c = dataset_list[d]
    for i in range(2):
      if coh_level == 'highDim':
        data = c[c.highDimCoh == coherence_list[i]]
        rate_array[d,i,0] = sum(data.response==1)/(len(data.rxtime) + 1)
        qlist_array[d,i,0,:] = quantile_list(data.rxtime[data.response==1],nquants)
        rate_array[d,i,1] = sum(data.response==3)/(len(data.rxtime) + 1)
        qlist_array[d,i,1,:] = quantile_list(data.rxtime[data.response==3],nquants)
      else:
        data = c[c.lowDimCoh == coherence_list[i]]
        rate_array[d,i,0] = sum(data.response==2)/(len(data.rxtime) + 1)
        qlist_array[d,i,0,:] = quantile_list(data.rxtime[data.response==2],nquants)
        rate_array[d,i,1] = sum(data.response==3)/(len(data.rxtime) + 1)
        qlist_array[d,i,1,:] = quantile_list(data.rxtime[data.response==3],nquants)
  mean_rate = np.mean(rate_array, axis = 0)
  mean_qlist = np.mean(qlist_array, axis = 0)

  # if coh_level == 'highDim':
  plt.plot(np.repeat(mean_rate[0,0], nquants), mean_qlist[0,0,:], ':v', color = color)
  plt.plot(np.repeat(mean_rate[0, 1], nquants), mean_qlist[0,1,:], '-v', color = color)
  plt.plot(np.repeat(mean_rate[1,0], nquants), mean_qlist[1,0,:], ':^', color = color)
  plt.plot(np.repeat(mean_rate[1, 1], nquants), mean_qlist[1,1,:], '-^', color = color)



def single_coh_hist(data,nquants):
  """Plots a histogram of the data and idenifies n quantiles in the plot."""
  plt.hist(data.rxtime,bins=50,density=True,range=(0,8), histtype = 'step')

  q_list = quantile_list(data.rxtime,nquants)
  plt.vlines(q_list,0,1,linestyles='--',colors='b')

  lines = q_list
  lines = np.append(lines,5)

def single_coh_hist_line_opt(data, nquants, color = 'black', plot_line = False):
  """Plots a histogram of the data and idenifies n quantiles in the plot."""
  plt.hist(data.rxtime,bins=50,density=True,range=(0,8), histtype = 'step', color = color)

  if plot_line == True:
    q_list = quantile_list(data.rxtime,nquants)
    plt.vlines(q_list,0,1,linestyles='--',colors='b')
    lines = q_list
    lines = np.append(lines,5)

def single_coh_hist_line_opt_mean(dataset_list, nquants, color = 'black', plot_line = False):
  """Goes thru the dataset list. Plots a histogram of the data and idenifies n quantiles in the plot."""
  rxtime_list = []
  for d in dataset_list:
    rxtime_list.append(d.rxtime)
  mean_rxtime = np.mean(rxtime_list, axis = 0)

  plt.hist(mean_rxtime,bins=50,density=True,range=(0,8), histtype = 'step', color = color)

  if plot_line == True:
    q_list = quantile_list(data.rxtime,nquants)
    plt.vlines(q_list,0,1,linestyles='--',colors='b')
    lines = q_list
    lines = np.append(lines,5)
  # #locs = []
  # for i in range(nquants):
  #   #locs.append()
  #   plt.text(np.mean([lines[i],lines[i+1]])-0.1,0.95,'{}'.format(i+1))

    

def descriptive_plts(chong_data):
  """Function plotting plots that only require the data!"""

def create_list_by_condition(dataset_list, highdimcoh, lowdimcoh):
  sieved_list = []
  for d in dataset_list:
    sieved_list.append(d[(d['highDimCoh'] == highdimcoh ) & (d['lowDimCoh'] == lowdimcoh)])
  return sieved_list

def all_plots_with_real(ppc):
  """Just like all_plots, but loads in and plots the real data first"""
  #Load in the empirical dataset
  df = pd.read_csv('/users/avo2/data/avo2/LANminimal/experiments/fullWithHSSMCoding.csv')
  # data = df[['rt', 'subj_idx', 'response', 'highDimCoh', 'lowDimCoh', 'irrDimCoh']]

      #Remove behavioural NaNs
  df = df[df['rxtime'] > 0.2]
    #Remove extra participants
  df = df[df['subj'] != 13 ]
  df = df[df['subj'] != 16 ]
  #Edit Coherence from [1,2] to [-1,1]
  # data['highDimCoh'] = data['highDimCoh'].replace(1,-1)
  # data['highDimCoh'] = data['highDimCoh'].replace(2,1)
  # data['lowDimCoh'] = data['lowDimCoh'].replace(1,-1)
  # data['lowDimCoh'] = data['lowDimCoh'].replace(2,1)
  # data['irrDimCoh'] = data['irrDimCoh'].replace(1,-1)
  # data['irrDimCoh'] = data['irrDimCoh'].replace(2,1)
  # data = data.rename({'rt': 'rxtime',
  #                   'subj_idx': 'subj'}, axis = 'columns')

  chong_data = df

  
  plt.subplot(2,5,1)
  single_coh_hist_line_opt( chong_data[(chong_data['highDimCoh']==-1) & (chong_data['lowDimCoh']==-1)], 5, plot_line=False)
  single_coh_hist_line_opt_mean(create_list_by_condition(ppc, -1, -1), 5, color = 'orange')
  # single_coh_hist_line_opt(ppc[(ppc['highDimCoh']==-1) & (ppc['lowDimCoh']==-1)], 5, color = 'orange')
  plt.legend(['empirical', 'posterior predictive'])
  plt.title('Reaction Time Quantiles \n by Coherence: LL',fontsize=18)
  plt.xlabel('Reaction Time (seconds)',fontsize=14)
  plt.ylabel('Density',fontsize=14)

  plt.subplot(2,5,2)
  plt.title('Reaction Time Quantiles \n by Coherence: LH',fontsize=18)
  plt.xlabel('Reaction Time (seconds)',fontsize=14)
  plt.ylabel('Density',fontsize=14)
  single_coh_hist_line_opt( chong_data[(chong_data['highDimCoh']==-1) & (chong_data['lowDimCoh']==1)], 5, plot_line=False)
  single_coh_hist_line_opt_mean(create_list_by_condition(ppc, -1, 1), 5, color = 'orange')
  plt.legend(['empirical', 'posterior predictive'])

  plt.subplot(2,5,6)
  err_plot(chong_data[(chong_data['highDimCoh']==-1) & (chong_data['lowDimCoh']==-1)],5, color = 'black')
  err_plot_mean(create_list_by_condition(ppc, -1, -1) ,5, color = 'orange')
  plt.title('Error Rate in RT Quantiles \n by Coherence: LL',fontsize=18)
  plt.xlabel('Reaction Time Quantile',fontsize=14)
  plt.ylabel('Error Rate',fontsize=14)
  plt.xticks([0,1,2,3,4],[1,2,3,4,5])
  plt.legend(['high dim error (empirical)', 'low dim error (empirical)', 'high dim error (ppc)', 'low dim error (ppc)'])
  


  plt.subplot(2,5,4)
  plt.title('Reaction Time Quantiles \n by Coherence: HH',fontsize=18)
  plt.xlabel('Reaction Time (seconds)',fontsize=14)
  plt.ylabel('Density',fontsize=14)
  single_coh_hist_line_opt(chong_data[(chong_data['highDimCoh']==1) & (chong_data['lowDimCoh']==1)], 5, plot_line=False)
  single_coh_hist_line_opt_mean(create_list_by_condition(ppc, 1, 1), 5, color = 'orange')
  plt.legend(['empirical', 'posterior predictive'])
 
  plt.subplot(2,5,3)
  plt.title('Reaction Time Quantiles \n by Coherence: HL',fontsize=18)
  plt.xlabel('Reaction Time (seconds)',fontsize=14)
  plt.ylabel('Density',fontsize=14)
  single_coh_hist_line_opt(chong_data[(chong_data['highDimCoh']==1) & (chong_data['lowDimCoh']==-1)], 5, plot_line=False)
  single_coh_hist_line_opt_mean(create_list_by_condition(ppc, 1, -1), 5, color = 'orange')
  plt.legend(['empirical', 'posterior predictive'])


  plt.subplot(2,5,7)
  err_plot(chong_data[(chong_data['highDimCoh']==-1) & (chong_data['lowDimCoh']==1)],5, color = 'black')
  err_plot_mean(create_list_by_condition(ppc, -1, 1) ,5, color = 'orange')
  plt.title('Error Rate in RT Quantiles \n by Coherence: LH',fontsize=18)
  plt.xlabel('Reaction Time Quantile',fontsize=14)
  plt.ylabel('Error Rate',fontsize=14)
  plt.xticks([0,1,2,3,4],[1,2,3,4,5])
  plt.legend(['high dim error (empirical)', 'low dim error (empirical)', 'high dim error (ppc)', 'low dim error (ppc)'])


  plt.subplot(2,5,8)
  err_plot(chong_data[(chong_data['highDimCoh']==1) & (chong_data['lowDimCoh']==-1)],5, color = 'black')
  err_plot_mean(create_list_by_condition(ppc, 1, -1) ,5, color = 'orange')
  plt.title('Error Rate in RT Quantiles \n by Coherence: HL',fontsize=18)
  plt.xlabel('Reaction Time Quantile',fontsize=14)
  plt.ylabel('Error Rate',fontsize=14)
  plt.xticks([0,1,2,3,4],[1,2,3,4,5])
  plt.legend(['high dim error (empirical)', 'low dim error (empirical)', 'high dim error (ppc)', 'low dim error (ppc)'])

  plt.subplot(2,5,9)
  err_plot(chong_data[(chong_data['highDimCoh']==1) & (chong_data['lowDimCoh']==1)],5, color = 'black')
  err_plot_mean(create_list_by_condition(ppc, 1, 1) ,5, color = 'orange')
  plt.title('Error Rate in RT Quantiles \n by Coherence: HH',fontsize=18)
  plt.xlabel('Reaction Time Quantile',fontsize=14)
  plt.ylabel('Error Rate',fontsize=14)
  plt.xticks([0,1,2,3,4],[1,2,3,4,5])
  plt.legend(['high dim error (empirical)', 'low dim error (empirical)', 'high dim error (ppc)', 'low dim error (ppc)'])


  plt.subplot(2,5,5)
  plt.title('High Dimension Quantile\nProbability Plot', fontsize = 18)
  single_qpp(chong_data, 5, 'black')
  single_qpp_mean(ppc, 5, 'orange' )
  plt.legend(['low coherence,\n incorrect','low coherence,\n correct','high coherence,\n incorrect','high coherence,\n correct'])
  plt.xlabel('Response Probability',fontsize=14)
  plt.ylabel('Reaction Time (seconds)',fontsize=14)

  plt.subplot(2,5,10)
  plt.title('Low Dimension Quantile\nProbability Plot', fontsize = 18)
  plt.ylabel('Reaction Time (seconds)')
  single_ (chong_data,5,'black','lowDim')
  single_qpp_mean(ppc,5,'orange','lowDim')
  plt.legend(['low coherence,\n incorrect','low coherence,\n correct','high coherence,\n incorrect','high coherence,\n correct'])
  plt.xlabel('Response Probability',fontsize=14)
  plt.ylabel('Reaction Time (seconds)',fontsize=14)
  
def all_plots(chong_data):
  """Function that plots all plots relevant to chong data study"""
  plt.subplot(2,5,1)
  single_coh_hist(chong_data[(chong_data['highDimCoh']==-1) & (chong_data['lowDimCoh']==-1)], 5)
  plt.title('Reaction Time Quantiles \n by Coherence: LL',fontsize=18)
  plt.xlabel('Reaction Time (seconds)',fontsize=14)
  plt.ylabel('Density',fontsize=14)
  plt.subplot(2,5,2)
  plt.title('Reaction Time Quantiles \n by Coherence: LH',fontsize=18)
  plt.xlabel('Reaction Time (seconds)',fontsize=14)
  plt.ylabel('Density',fontsize=14)
  single_coh_hist(chong_data[(chong_data['highDimCoh']==-1) & (chong_data['lowDimCoh']==1)],5)
  plt.subplot(2,5,3)
  plt.title('Reaction Time Quantiles \n by Coherence: HL',fontsize=18)
  plt.xlabel('Reaction Time (seconds)',fontsize=14)
  plt.ylabel('Density',fontsize=14)
  single_coh_hist(chong_data[(chong_data['highDimCoh']==1) & (chong_data['lowDimCoh']==-1)],5)
  plt.subplot(2,5,4)
  plt.title('Reaction Time Quantiles \n by Coherence: HH',fontsize=18)
  plt.xlabel('Reaction Time (seconds)',fontsize=14)
  plt.ylabel('Density',fontsize=14)
  single_coh_hist(chong_data[(chong_data['highDimCoh']==1) & (chong_data['lowDimCoh']==1)],5)

  plt.subplot(2,5,6)
  err_plot(chong_data[(chong_data['highDimCoh']==-1) & (chong_data['lowDimCoh']==-1)],5)
  plt.title('Error Rate in RT Quantiles \n by Coherence: LL',fontsize=18)
  plt.xlabel('Reaction Time Quantile',fontsize=14)
  plt.ylabel('Error Rate',fontsize=14)
  plt.xticks([0,1,2,3,4],[1,2,3,4,5])
  plt.legend(['High Dimension Error','Low Dimension Error'])

  plt.subplot(2,5,7)
  err_plot(chong_data[(chong_data['highDimCoh']==-1) & (chong_data['lowDimCoh']==1)],5)
  plt.title('Error Rate in RT Quantiles \n by Coherence: LH',fontsize=18)
  plt.xlabel('Reaction Time Quantile',fontsize=14)
  plt.ylabel('Error Rate',fontsize=14)
  plt.xticks([0,1,2,3,4],[1,2,3,4,5])
  plt.legend(['High Dimension Error','Low Dimension Error'])


  plt.subplot(2,5,8)
  err_plot(chong_data[(chong_data['highDimCoh']==1) & (chong_data['lowDimCoh']==-1)],5)
  plt.title('Error Rate in RT Quantiles \n by Coherence: HL',fontsize=18)
  plt.xlabel('Reaction Time Quantile',fontsize=14)
  plt.ylabel('Error Rate',fontsize=14)
  plt.xticks([0,1,2,3,4],[1,2,3,4,5])
  plt.legend(['High Dimension Error','Low Dimension Error'])


  plt.subplot(2,5,9)
  err_plot(chong_data[(chong_data['highDimCoh']==1) & (chong_data['lowDimCoh']==1)],5)
  plt.title('Error Rate in RT Quantiles \n by Coherence: HH',fontsize=18)
  plt.xlabel('Reaction Time Quantile',fontsize=14)
  plt.ylabel('Error Rate',fontsize=14)
  plt.xticks([0,1,2,3,4],[1,2,3,4,5])
  plt.legend(['High Dimension Error','Low Dimension Error'])


  plt.subplot(2,5,5)
  plt.title('High Dimension Quantile\nProbability Plot', fontsize = 18)
  single_qpp(chong_data,5,'black')
  # single_qpp_mean(ppc,5,'orange')
  plt.legend(['Low coherence,\n incorrect','Low coherence,\n correct','High coherence,\n incorrect','High coherence,\n correct'])
  plt.xlabel('Response Probability',fontsize=14)
  plt.ylabel('Reaction Time (seconds)',fontsize=14)

  plt.subplot(2,5,10)
  plt.title('Low Dimension Quantile\nProbability Plot', fontsize = 18)
  plt.ylabel('Reaction Time (seconds)')
  single_qpp(chong_data,5, 'black', 'lowDim')
  # single_qpp_mean(ppc,5,'orange', 'lowDim')
  plt.legend(['Low coherence,\n incorrect','Low coherence,\n correct','High coherence,\n incorrect','High coherence,\n correct'])
  plt.xlabel('Response Probability',fontsize=14)
  plt.ylabel('Reaction Time (seconds)',fontsize=14)
    
    

