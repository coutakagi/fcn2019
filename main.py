'''
Created on 2019/11/17

@author: Takagi
'''

if __name__ == '__main__':
    print("start")
    import os
    from setting import setting 
    set=setting() 
    set.input_folder=os.path.join(os.path.dirname(__file__),"inputdata")
    set.output_folder=os.path.join(os.path.dirname(__file__),"outputdata")
    
    from networkmodel import networksimulation
    networksimulation(set)
    
    print("End all process")
    
