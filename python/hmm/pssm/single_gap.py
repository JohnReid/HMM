#
# Copyright John Reid 2008
#

"""
Code to build single gap PSSM models using HMMs.
"""

import hmm

class MotifModelPositionMap(object):
    """
    Maps between positions in a motif and positions in a model representing the motif.
    """
    def __init__(self, K):
        """
        @arg K: The length of the motif.
        """
        self.K = K
        "The length of the ungapped motif."

        self.N = self.K*2
        "The number of states in the model."

    def model_idx(self, k, positive):
        "@return: The index into the model's states that represents the k'th base in the given orientation."
        if positive:
            return k
        else:
            return k+self.N/2

    def motif_position(self, n):
        "@return: (k, positive) where k is the index into the motif and positive represents the orientation."
        if n < self.K:
            return n, True
        else:
            return n-self.K, False




class SingleGappedPssmBuilder(object):
    """
    Knows how to build a single gapped pssm.
    """

    def __init__(
      self,
      K,
      gap_index,
      markov_order=0,
      M=4
    ):
        """
        @arg K: The number of positions in the gapped pssm.
        @arg gap_index: The index of the position in the motif which is the gap. I.e. 1 would
        place the gap after the first base.
        @arg markov_order: the Markov order of the model.
        @arg M: the output alphabet size.
        """
        self.map = MotifModelPositionMap(K)
        "Maps positions in the pssm to the model and back again."

        self.K = K
        "The number of positions in the gapped pssm."

        self.gap_index = gap_index
        """
        The index of the position in the motif which is the gap. I.e. 1 would
        place the gap after the first base.
        """
        if self.gap_index < 1 or self.gap_index > self.K-2:
            raise RuntimeError('Gap must be in middle of motif.')

        self.markov_order = markov_order
        "The Markov order of the model."

        self.M = M
        "The output alphabet size."

    def num_states(self):
        "The number of states in a model of the shape defined by this builder."
        return self.map.N

    def _set_emissions(self, model, positive_state, negative_state, emissions):
        """
        Set the emissions of the states to be the complement of each other.
        """
        assert len(emissions) == self.M
        emission_params = [model.add_parameter(p_e) for p_e in emissions]
        for m, p in enumerate(emission_params):
            positive_state.b[m] = p
            negative_state.b[-m-1] = p

    def get_emissions_and_gap_probabilities(self, model, offset=0):
        "@return: emissions, gap_probabilities"
        import numpy
        emissions = model.B[offset:offset+self.K]
        gap_probabilities = numpy.ones(self.K)
        gap_probabilities[self.gap_index] = model.A[offset+self.gap_index-1, offset+self.gap_index]
        return emissions, gap_probabilities

    def create(
      self,
      p_gap,
      emissions
    ):
        """
        @arg p_gap: the probability of a gap.
        @arg emissions: the emission distributions of the bases (including the gap).
        @returns: A tuple (model, positive_start, positive_end, negative_start, negative_end).
        The model is defined by its states and
        includes both the motif and its reverse complement. positive_start indexes the first state in the
        positive motif and negative_start indexes the first state in the negative motif.
        """
        # create the model
        model = hmm.ModelByStates(M=self.M, markov_order=self.markov_order)

        # add enough states to the models
        for k in xrange(self.num_states()):
            model.add_state()

        # link the states
        transition_param_one = model.add_parameter(1.)
        transition_param_gap = model.add_parameter(p_gap)
        transition_param_not_gap = model.add_parameter(1. - p_gap)
        # positive transitions
        for k in xrange(self.K-1):
            if k+1 != self.gap_index:
                # this is not the base before the gap
                # so just connect to next base
                model.states[self.map.model_idx(k, True)].add_successor(
                  model.states[self.map.model_idx(k+1, True)],
                  transition_param_one
                )
            else:
                # this is the base before the gap
                # so connect to the gap and the base after the gap
                model.states[self.map.model_idx(k, True)].add_successor(
                  model.states[self.map.model_idx(k+1, True)],
                  transition_param_gap
                )
                model.states[self.map.model_idx(k, True)].add_successor(
                  model.states[self.map.model_idx(k+2, True)],
                  transition_param_not_gap
                )
        # negative transitions
        for k in xrange(self.K-1):
            if k != self.gap_index:
                # this is not the base before the gap
                # so just connect to next base
                model.states[self.map.model_idx(k+1, False)].add_successor(
                  model.states[self.map.model_idx(k, False)],
                  transition_param_one
                )
            else:
                # this is the base before the gap
                # so connect to the gap and the base after the gap
                model.states[self.map.model_idx(k+1, False)].add_successor(
                  model.states[self.map.model_idx(k, False)],
                  transition_param_gap
                )
                model.states[self.map.model_idx(k+1, False)].add_successor(
                  model.states[self.map.model_idx(k-1, False)],
                  transition_param_not_gap
                )

        # fill in the emission distributions
        assert len(emissions) == self.K
        for k, base_emissions in enumerate(emissions):
            self._set_emissions(
              model,
              model.states[self.map.model_idx(k, True)],
              model.states[self.map.model_idx(k, False)],
              base_emissions
            )

        return (
          model,
          [
            self.map.model_idx(k=0, positive=True),
            self.map.model_idx(k=self.K-1, positive=False),
          ],
          [
            self.map.model_idx(k=self.K-1, positive=True),
            self.map.model_idx(k=0, positive=False),
          ]
        )




def extend_model(model, extension):
    """
    Copies all the states, and emission and transition parameters from the extension model into
    the model.
    @arg model: The model to be extended.
    @arg extension: The model that is the extension.
    """
    assert model.M == extension.M

    # add a parameter to model for each parameter in extension
    extended_parameters = [
      model.add_parameter(p)
      for p in extension.parameters
    ]

    # add a state to model for each state in extension
    extended_states = [
      model.add_state()
      for state in extension.states
    ]
    state_map = dict((state, i) for i, state in enumerate(extension.states))

    # add the transitions to the model
    for state in extension.states:
        extended_state = extended_states[state_map[state]]
        for successor in state.successors:
            extended_state.add_successor(
              extended_states[state_map[successor.state]],
              extended_parameters[successor.a.idx]
            )

    # add the emissions to the model
    for state in extension.states:
        extended_state = extended_states[state_map[state]]
        for m, b in enumerate(state.b):
            extended_state.b[m] = extended_parameters[b.idx]





def simplest_background_model(markov_order=0, M=4):
    model = hmm.ModelByStates(M=M, markov_order=markov_order)
    state = model.add_state()
    state.add_successor(state, model.add_parameter(1.))
    for m in xrange(M):
        state.b[m] = model.add_parameter(1./M)
    return model



def add_to_simple_background_model(model, in_states, out_states, p_binding_site):
    """
    Create a simple background model and extend it with a copy of the given model.

    @arg model: The model to extend the background model with.
    @arg in_states: Indices of those states the background model should transition to.
    @arg out_states: Indices of those states that should transition back to the background model.
    """
    complete_model = simplest_background_model(model.converter.n, model.M)
    extend_model(complete_model, model)
    # link the background model to the positive and negative parts of the single gapped pssm
    binding_site_transition_param = complete_model.add_parameter(p_binding_site/len(in_states))
    back_to_bg_transition_param = complete_model.add_parameter(1.)
    for in_state in in_states:
        complete_model.states[0].add_successor(
          complete_model.states[1+in_state],
          binding_site_transition_param
        )
    for out_state in out_states:
        complete_model.states[1+out_state].add_successor(
          complete_model.states[0],
          back_to_bg_transition_param
        )
    complete_model.states[0].pi = complete_model.add_parameter(1.)
    return complete_model




if '__main__' == __name__:
    import numpy

    # build a single gapped pssm with some random emissions
    builder = SingleGappedPssmBuilder(K=6, gap_index=1, markov_order=0, M=4)
    emissions = numpy.array(
      [
        hmm.dirichlet_draw(numpy.ones(builder.M) * .1)
        for k in xrange(builder.K)
      ]
    )
    emissions[builder.gap_index] = hmm.dirichlet_draw(numpy.ones(builder.M) * .3)
    model_by_states, in_states, out_states = builder.create(
      p_gap=.6,
      emissions=emissions
    )

    # create a background model and add the single gapped pssm to it
    complete_model = add_to_simple_background_model(
      model_by_states,
      in_states,
      out_states,
      p_binding_site=.01)

    # convert to other type of model
    model = hmm.as_model(complete_model)

    # write as a graph
    hmm.graph_as_svg(
      model,
      'single-gapped-hmm',
      graphing_keywords={'include_emissions':False},
      neato_properties={'-Elen':2}
    )

    # get the emissions and gap probabilities and write a logo
    emissions_copy, gap_probs = builder.get_emissions_and_gap_probabilities(model, offset=1)
    assert (emissions_copy - emissions).sum() < 1e-10
    import hmm.pssm.logo as logo
    image = logo.pssm_as_image(emissions, transparencies=gap_probs)
    image.save("single-gapped-pssm-logo.png", "PNG")
