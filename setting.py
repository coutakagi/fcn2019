'''
Created on 2019/11/17

@author: Takagi
'''

class setting(object):
    "setting parameters and default values"
    
    #set default values
    
    #input/output parameters
    input_folder=""#input folder of connectivity data
    output_folder=""#output folder for simulation 
    matrix_size=177#input matrix size
    
    #threshld values for the connectivity matrix (Threshold=Average+n*SD)
    start_threshold=-1.0#start values for threshold
    end_threshold=2.0#start values for threshold
    threshold_step=0.1#interval step of threshold range
    
    # random activation signal input parameters
    batch_size=10#batch size for simulation
    repeat_num=100#repeat number for each batch
    activation_probability=0.05#probability for positive/negative activation
    
    deb=1#print out result if debug

    def __init__(self):
        '''
        Constructor
        '''
        
        
        
        
    
