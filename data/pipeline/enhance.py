'''
This class exports helping features on chunks.
These features rely on pre-processed chunk data
and thus are best built using loaded TF data.
The embedding feature tells whether a chunk is embedded in another chunk
of the same kind. The exported feature is a string on 
chunk nodes with a value of true or false.
Time roles are also added on timephrase chunks with embedded quantifiers.
''' 

from tf.fabric import Fabric
import collections

class Enhance:
    
    def __init__(self, locs, metadata):
        
        print('preparing to enhance chunks...')
        print('\tloading fresh TF data...')
        
        # intialize Text-Fabric methods
        locas = [locs['bhsa'], locs['output']]
        TF = Fabric(locations=locas, silent=True)
        load_features = ['label', 'function', 'role']
        self.bhsa = TF.load(' '.join(load_features), silent=True)    
        self.output_dir = locs['output']
        self.meta = metadata
        self.nodeFeatures = collections.defaultdict(lambda:collections.defaultdict())
        self.edgeFeatures = collections.defaultdict(lambda:collections.defaultdict())
        
        # add existing role edges
        F, E = self.bhsa.F, self.bhsa.E
        for w in F.otype.s('word'):
            if E.role.f(w):
                self.edgeFeatures['role'][w] = dict(E.role.f(w))
        
    def embeddings(self):
        '''
        Identifies chunks contained in another chunk
        of the same kind (i.e. label). If it is, feature
        `embed`=true; else =false. 
        '''
        
        F, L = self.bhsa.F, self.bhsa.L
        
        for chunk in F.otype.s('chunk'):
            # check all embedding chunks for a 
            # label that matches this chunk's label
            # if matched, embed=true
            label = F.label.v(chunk)
            contained = [ch for ch in L.u(chunk, 'chunk')
                            if F.label.v(chunk) == label]
            if contained:
                self.nodeFeatures['embed'][chunk] = 'true'
            else:
                self.nodeFeatures['embed'][chunk] = 'false'
                
        # report new features to be added
        new_features = len(self.nodeFeatures['embed'])
        feature_breakdown = collections.Counter(feat for node, feat 
                                                    in self.nodeFeatures['embed'].items())
        
        print(f'{new_features} new embed features added...')
        for feat, count in feature_breakdown.most_common():
            print(f'\t{count} embed features with value of {feat}...')
            
            
    def quanttimes(self):
        '''
        Adds time roles to timephrase chunks with 
        quantifiers, also adds edge relation from
        the quantifiers to the timephrase chunk.
        '''
        
        print('Adding new quanttime role data...')
        
        F, E, L = self.bhsa.F, self.bhsa.E, self.bhsa.L
        
        new_roles = 0
        
        for chunk in F.label.s('quant_NP'):
            
            # isolate time chunk
            timephrase = [tp for tp in L.u(chunk, 'chunk')
                             if F.label.v(tp) == 'timephrase']
            
            # skip non-times
            if not timephrase:
                continue
                
            # add time roles
            for word, role in E.role.t(chunk):
                if role == 'subs':
                    self.edgeFeatures['role'][word].update({timephrase[0]:'time'})
                    new_roles += 1
                elif role == 'quant':
                    self.edgeFeatures['role'][word].update({timephrase[0]:'quant'})
                    new_roles += 1
                    
        print(f'\t{new_roles} new time roles added...')
            
    def export(self):
        '''
        Exports new TF data.
        '''
    
        # export to TF files
        metadata = {'': self.meta,
                    'embed': {'valueType':'str',
                              'description':'a feature of chunk objects, describing whether they are embedded in another chunk of the same label'
                             },
                    'role':{'edgeValues':True,
                            'valueType':'str', 
                            'description':'role of a word in a chunk object'}
                   }
        
        print('exporting tf...')
        TFexport = Fabric(locations=self.output_dir, silent=True)
        TFexport.save(metaData=metadata, nodeFeatures=self.nodeFeatures, edgeFeatures=self.edgeFeatures)
        print('SUCCESS')