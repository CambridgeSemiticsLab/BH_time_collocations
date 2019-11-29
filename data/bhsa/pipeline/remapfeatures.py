import collections, os
from tf.fabric import Fabric

class RemapFeatures:
    
    '''
    -- Correcting Features On Nodes in BHSA -- 
    
    --Phrase Functions--
    Phrase functions are redrawn based on the manual 
    corrections stored in the dictionary named `newfunctions`
    
    -- Subphrase Relations --
    A number of corrections are made to subphrase relations.
    '''

    def __init__(self, locs, metadata):
    
        '''
        Stages the changes to be made.
        Provides attributes for reporting.
        '''
        
        bhsa_dir, export_dir = locs['bhsa'], locs['output']
        
        # load TF api + functions
        TF = Fabric(locations=bhsa_dir, silent=True)
        self.api = TF.load('function rela', silent=True)
        self.export_dir = export_dir
        
        # features to be exported
        self.meta = {'':metadata,
                     'note':{'valueType':'str',
                             'description':'notes on objects for tracking issues throughout my research'},
                     'function':{'valueType':'str'},
                     'rela':{'valueType':'str'}
                    }
        
        # store features here
        self.nodeFeatures = collections.defaultdict(lambda:collections.defaultdict())
        
    def remap_functions(self):
        '''
        Remaps select phrase functions in the BHSA.
        '''
        
        # shortform TF methods
        F = self.api.F
        
        # map corrections/changes here
        newfunctions = {
            849296:'Loca',
            825329:'Loca',
            828081:'Cmpl',
            774349:'Adju',
            774352:'Adju',
            775948:'Adju',
            775985:'Adju',
            876172:'Adju',
            881665:'Objc', # phrase belongs with previous as adjectival element
        }
        
        # reconstruct default BHSA phrase function mappings
        self.nodeFeatures['function'] = {ph:F.function.v(ph) for ph in F.otype.s('phrase')}
        
        print('Updating phrase functions...')
        self.nodeFeatures['function'].update(newfunctions)
        
        print('\tadding notes feature to new functions...')
        for phrase, new_funct in newfunctions.items():
            self.nodeFeatures['note'][phrase] = f'function changed from Time to {new_funct}'
            
        self.newfunctions = newfunctions # make available for reporting in the pipeline notebook
        
    def remap_relations(self):
        '''
        Remap select subphrase relations in the BHSA.
        '''
        F, N = self.api.F, self.api.N
        
        newrelas = {1308697: 'adj',
                    1341455: 'adj',
                    1351683: 'adj',
                    1351797: 'adj',
                    1353159: 'adj',
                    1353279: 'adj',
                    1353681: 'adj'}
        
        # reconstruct default relations
        self.nodeFeatures['rela'] = {node:F.rela.v(node) for node in N() 
                                         if F.rela.v(node)
                                         and F.rela.v(node) != 'NA'}
        
        print('Updating subphrase relations...')
        self.nodeFeatures['rela'].update(newrelas)
    
        for sp, new_funct in newrelas.items():
            self.nodeFeatures['note'][sp] = f'relation changed from {F.rela.v(sp)} to {new_funct}'
        
        self.newrelas = newrelas
    
    def execute(self):
        '''
        Executes the remappings.
        '''        
        print('\texporting tf...')
        # export the .tf files
        TFexport = Fabric(locations=self.export_dir, silent=True)
        TFexport.save(metaData=self.meta, nodeFeatures=self.nodeFeatures)

        print('\tSUCCESS')
