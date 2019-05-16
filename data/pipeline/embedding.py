'''
This class exports a simple helping feature on chunks
that tells whether a chunk is embedded in another chunk
of the same kind. The exported feature is a string on 
chunk nodes with a value of true or false.
'''

from tf.fabric import Fabric
import collections

class EmbedDetector:
    
    def __init__(self, bhsa_dir, output_dir, metadata):
        
        locas = [bhsa_dir, output_dir]
        TF = Fabric(locations=locas, silent=True)
        load_features = ['label',]
        self.bhsa = TF.load(' '.join(load_features), silent=True)    
        self.output_dir = output_dir
        self.meta = metadata
        self.nodeFeatures = collections.defaultdict(lambda:collections.defaultdict())
        
    def execute(self):
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
            
        # export to TF files
        metadata = {'': self.meta,
                    'embed': {'valueType':'str',
                              'description':'a feature of chunk objects, describing whether they are embedded in another chunk of the same label'
                             }}
        
        print('exporting tf...')
        TFexport = Fabric(locations=self.output_dir, silent=True)
        TFexport.save(metaData=metadata, nodeFeatures=self.nodeFeatures)
        print('SUCCESS')