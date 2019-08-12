'''
This class exports enhancements based on the final
underlying dataset.

--Helping Features on Chunks--
• embedding [embeddings] - tells whether a chunk is embedded in another chunk
of the same kind. The feature is a string on chunk nodes 
with a value of "true" or "false".
• role [quanttimes] - Time roles are added on timephrase chunks with embedded quantifiers.

--Edited BHSA Features
• vt [add_weqatal] - the verb tense feature is edited to include weqatal
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
        load_features = ['label', 'function', 'role', 'pdp', 'vt',
                         'lex', 'mother']
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
                            if F.label.v(ch) == label]
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
    
    def add_weqatal(self):
        '''
        Adds weqatal parsings to the vt (verb tense) feature.
        '''
        
        print('Adding weqatals to vt (verb tense) feature...')
        
        F, E, L = self.bhsa.F, self.bhsa.E, self.bhsa.L
        
        def get_grandma(clause_atom):
            '''
            Recursively climbs up a clause's ancestorial tree.
            Stops upon identifying either wayyiqtol
            or a yiqtol|impv grand(mother).
            Returns a string of the ancestor's tense, or nothing.
            '''

            this_verb = next((F.vt.v(w) for w in L.d(clause_atom) if F.pdp.v(w)=='verb'), '')
            mother = next((m for m in E.mother.f(clause_atom)), 0)
            mom_verb = next((F.vt.v(w) for w in L.d(mother) if F.pdp.v(w)=='verb'), '')    

            if mom_verb in {'wayq', 'impf', 'impv'}:
                return mom_verb
            elif not mother:
                return this_verb
            else:
                return get_grandma(mother)
        
        # Begin weqatal calculation
        new_weqtls = []
        print('\tre-calculating verb tenses...')
        for word in F.otype.s('word'):
            
            # skip non-verbs
            if F.vt.v(word) == 'NA':
                continue
            
            # calculate presence of the weqatal
            if F.vt.v(word) == 'perf' and F.lex.v(word-1) == 'W':
                
                # get tense of the ancestor of the verb's clause
                clause = L.u(word, 'clause_atom')[0]
                qatal_ancestor = get_grandma(clause)
                
                # check for whether ancestor triggers weqatal analysis
                if qatal_ancestor in {'impf', 'impv'}:
                    self.nodeFeatures['vt'][word] = 'weqt' # change tense to weqt
                    new_weqtls.append(word)
                else:
                    self.nodeFeatures['vt'][word] = F.vt.v(word) # no change on tense
                
            # re-add tense unchanged
            else:
                self.nodeFeatures['vt'][word] = F.vt.v(word)
                
        print(f'\t{len(new_weqtls)} qatals changed to weqatals...')
                
        
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
                            'description':'role of a word in a chunk object'},
                    'vt':{'valueType':'str',
                          'description':'These are default tense values from the BHSA but edited to include weqatal where it is present.'
                         }
                   }
        
        print('exporting tf...')
        TFexport = Fabric(locations=self.output_dir, silent=True)
        TFexport.save(metaData=metadata, nodeFeatures=self.nodeFeatures, edgeFeatures=self.edgeFeatures)
        print('SUCCESS')