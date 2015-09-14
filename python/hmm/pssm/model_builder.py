#
# Copyright John Reid 2007, 2008
#


"""
Code to build HMM models of PSSMs of various Markov orders.
"""

import hmm, numpy, pickle


def _load_model_parameters(model, f):
    parameters = pickle.load(f)
    if len(parameters) != len(model.parameters):
        raise RuntimeError('Different number of parameters in file (%d) and in model (%d)' % (len(parameters), len(model.parameters)))
    for i, t in enumerate(parameters):
        model.parameters[i] = t

def _dump_model_parameters(model, f):
    pickle.dump(list(model.parameters), f)



class ModelBuilder(object):
    """
    Helps to build models over order n alphabets
    """

    def __init__(self, order, alphabet_size=4):
        self.order = order
        self.alphabet_size = alphabet_size
        self.M = alphabet_size ** (order + 1)
        self.converter = hmm.MarkovOrderConverter(alphabet_size, order)

    def new_model_by_states(self):
        "@return: A new hmm.ModelByStates."
        return hmm.ModelByStates(self.M, self.order)

    def add_fully_parameterised_state(self, model, pi = 1.0, emission_dist = None):
        """
        Adds a state with a different emission parameter for each of self.M
        possible order-n output characters
        """

        # add the state
        state = model.add_state(pi = model.add_parameter(pi))

        # allocate the parameters to the emissions
        for m in xrange(self.M):
            if None != emission_dist:
                assert len(emission_dist) == self.M
                state.b[m] = model.add_parameter(emission_dist[m])
            else:
                state.b[m] = model.add_parameter( 1.0 / self.M )

    def add_order_0_parameterised_state(self, model, pi = 1.0, emission_dist = None):
        """
        Adds a state with shared emission parameters for each of the self.M
        possible order-n output characters that represent the same order-0 character
        """

        # add the state
        state = model.add_state(pi = model.add_parameter(pi))

#               from IPython.Shell import IPShellEmbed
#               ipshell = IPShellEmbed()
#               ipshell()

        # add the parameters for the emissions
        if None != emission_dist:
            assert len(emission_dist) == self.alphabet_size
            params = [ model.add_parameter( x ) for x in emission_dist ]
        else:
            params = [ model.add_parameter( 1.0 / self.M ) for i in xrange(self.alphabet_size) ]

        # allocate the parameters to the emissions
        for m in xrange(self.M):
            state.b[m] = params[m % self.alphabet_size]

        return state

    def add_order_0_rev_comp_state(self, model, forward_state, pi = 1.0):
        """
        Adds a state with shared emission parameters for each of the self.M
        possible order-n output characters that represent the same order-0 character.
        This state's output is the reverse complement of the given state
        """

        # add the state
        state = model.add_state(pi = model.add_parameter(pi))

        # allocate the parameters to the emissions
        for m in xrange(self.M):
            state.b[m] = forward_state.b[self.alphabet_size - 1 - (m % self.alphabet_size)]

        return state

    def create_uniform_background_model(self):
        """
        @return: A HMM with one mosaic with uniform emission probabilities.
        """
        model = hmm.ModelByStates(self.M, self.order)
        self.add_fully_parameterised_state(
                model,
                emission_dist = numpy.ones(self.M)/4.
        )
        transition_param = model.add_parameter(1.0)
        model.states[0].add_successor(model.states[0], transition_param)
        return model

    def create_background_mosaic_model(self, num_mosaics, p_transition, dirichlet_prior_strength):
        """
        Create a mosaic model
        """
        model = hmm.ModelByStates(self.M, self.order)
        transition_param = model.add_parameter(p_transition)
        no_transition_param = model.add_parameter(1.0 - p_transition)
        for i in xrange(num_mosaics):
            self.add_fully_parameterised_state(
                    model,
                    emission_dist = hmm.dirichlet_draw(numpy.ones(self.M)*dirichlet_prior_strength)
            )
        for state_1 in model.states:
            for state_2 in model.states:
                if state_1 == state_2: state_1.add_successor(state_2, no_transition_param)
                else: state_1.add_successor(state_2, transition_param)
        return model

    def load_background_mosaic_model(self, f):
        '''
        Load a background model from the given file (or filename)
        '''
        if isinstance(f, str):
            f = open(f)

        # how many mosaics?
        num_mosaics = pickle.load(f)

        # create a model of the desired structure
        model = self.create_background_mosaic_model(num_mosaics, 0.0, 1.0)

        # load the parameters into the model
        _load_model_parameters(model, f)

        return model

    def dump_background_mosaic_model(self, model, f):
        '''
        Dump a background model into the given file (or filename)
        '''
        if isinstance(f, str):
            f = open(f, 'w')

        # how many mosaics?
        pickle.dump(model.N, f)

        # load the parameters into the model
        _dump_model_parameters(model, f)

def create_background_model(order, N):
    return hmm.pssm.ModelBuilder(order).create_background_mosaic_model(N, .01, 100.0)
