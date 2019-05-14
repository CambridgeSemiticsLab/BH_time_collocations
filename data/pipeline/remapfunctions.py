import collections, os
from tf.fabric import Fabric

class RemapFunctions:
    
    '''
    -- Correcting Time Functions in BHSA -- 
    In this class, phrase functions are redrawn
    based on the manual corrections stored in the 
    dictionary named `newfunctions`
    
    The Class simply wraps the necessary changes and
    also provides features needed for a report before and 
    after changes are affected.
    
    OUTPUT:
        * function.tf
        * note.tf
    '''

    def __init__(self, bhsa_dir, export_dir, metadata):
    
        '''
        Stages the changes to be made.
        Provides attributes for reporting.
        '''

        # load TF api + functions
        TF = Fabric(locations=bhsa_dir, silent=True)
        self.api = TF.load('function', silent=True)
        self.export_dir = export_dir

        # map corrections/changes here
        self.newfunctions = {849296:'Loca',
                             825329:'Loca',
                             828081:'Cmpl',
                             774349:'Adju',
                             774352:'Adju',
                             775948:'Adju',
                             775985:'Adju',
                             876172:'Adju',}
        
        # features to be exported
        self.meta = {'':metadata,
                     'note':{'valueType':'str',
                             'description':'notes on objects for tracking issues throughout my research'},
                     'function':{'valueType':'str'}}
        
    def execute(self):
        '''
        Executes the remappings.
        '''
        api = self.api
        newfunctions = self.newfunctions
        
        # build new phrase function features
        # phrase2funct is a dict with phrase node to function strip mappings
        # the alteration is made by simply updating phrase2funct with the newfunctions dict
        phrase2funct = dict((ph, api.F.function.v(ph)) for ph in api.F.otype.s('phrase'))
        note = {}
        
        print('Updating phrase functions...')
        phrase2funct.update(newfunctions)
        
        print('\tadding notes feature to new functions...')
        for phrase, new_funct in newfunctions.items():
            note[phrase] = f'function changed from Time to {new_funct}'

        print('\texporting tf...')
        # export the .tf files
        TFexport = Fabric(locations=self.export_dir, silent=True)
        TFexport.save(metaData=self.meta, nodeFeatures={'function':phrase2funct, 'note':note})

        print('\tSUCCESS')