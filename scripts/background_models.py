#
# Copyright John Reid 2007
#

"""
Code to save and access background models trained on the TRANSFAC chip-chip data
"""



import hmm, hmm.pssm, numpy, os.path

from fragments import *
from sequences import *

class BackgroundModelCache(object):
    _dir = 'background-models'
    _cache = {}

    if not os.access(_dir, os.X_OK):
        os.makedirs(_dir)

    def index_as_string(self, index):
        return 'order_%d-num_mosaics_%d-%s' % index

    def index(self, order, num_mosaics, fragment):
        return (order, num_mosaics, fragment)

    def file_for(self, index):
        return 'bg-%s.pickle' % self.index_as_string(index)

    def path_for(self, index):
        return os.path.join(self._dir, self.file_for(index))

    def load(self, index):
        order, num_mosaics, fragment = index
        filename = self.path_for(index)
        if not os.access(filename, os.R_OK):
            raise RuntimeError('Have no model for %s' % self.index_as_string(index))
        builder = hmm.pssm.ModelBuilder(order)
        return hmm.as_model(builder.load_background_mosaic_model(filename))

    def get_model(self, order, num_mosaics, fragment):
        index = self.index(order, num_mosaics, fragment)
        if index in self._cache:
            return self._cache[index]
        return self.load(index)

    def save_model(self, model, order, num_mosaics, fragment):
        index = self.index(order, num_mosaics, fragment)
        self._cache[index] = model
        hmm.pssm.ModelBuilder(order).dump_background_mosaic_model(hmm.as_state_model(model), self.path_for(index))

_global_cache = None
def global_background_model_cache():
    global _global_cache
    if None == _global_cache:
        _global_cache = BackgroundModelCache()
    return _global_cache

def background_model(order, N):
    return hmm.as_state_model(global_background_model_cache().load((order, N, 'T00594')))

if '__main__' == __name__:
    try:
        import cookbook
        cookbook.make_current_process_nice()
    except:
        print 'Could not set process priority'


    def bw_callback(LL): pass; #print 'LL: %.3f' % LL
    order = 0
    num_mosaics = 1
    cache = BackgroundModelCache()

    print '  fragment      order  # mosaics LL/base'
    for fragment in all_fragments:
        seqs = seqs_for_fragment(fragment)
        model_builder = hmm.pssm.ModelBuilder(order)
        training_sequences = [ model_builder.converter.to_order_n(s) for s in seqs ]
        known_bases = sum(model_builder.converter.num_known_bases_order_n(s) for s in training_sequences)
        print '%10s %10d %10d' % (fragment, order, num_mosaics),
        try:
            model = cache.get_model(order, num_mosaics, fragment)
        except:
            model_by_states = model_builder.create_background_mosaic_model(num_mosaics, 0.01, 100.0)
            model = hmm.model_states_2_model(model_by_states)
            tolerance = 1e-4 * known_bases
            model.baum_welch(
                    training_sequences,
                    tolerance = tolerance,
                    callback = bw_callback
            )
            cache.save_model(model, order, num_mosaics, fragment)
        LL = sum(model.forward(s)[0] for s in training_sequences)
        print LL/known_bases
