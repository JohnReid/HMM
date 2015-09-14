#
# Copyright John Reid 2007, 2008
#

"""
Code to test model building.
"""


import hmm, hmm.pssm, unittest, infpy

class TestModelBuilder(unittest.TestCase):
    def test_pickling(self):
        model_builder = hmm.pssm.ModelBuilder(3)
        model = model_builder.create_background_mosaic_model(3, .01, 1.0)
        model_builder.dump_background_mosaic_model(model, 'model.pickle')
        model_copy = model_builder.load_background_mosaic_model('model.pickle')

        converted_model = hmm.model_states_2_model(model)
        converted_model_copy = hmm.model_states_2_model(model_copy)

        assert converted_model.A.all() == converted_model_copy.A.all()
        assert converted_model.B.all() == converted_model_copy.B.all()
        assert converted_model.pi.all() == converted_model_copy.pi.all()

    def test_order_0_states(self):
        model_builder = hmm.pssm.ModelBuilder(2)
        model = model_builder.create_background_mosaic_model(3, .01, 1.0)
        positive = model_builder.add_order_0_parameterised_state(model, emission_dist=[.1,.2,.3,.4])
        negative = model_builder.add_order_0_rev_comp_state(model, positive)
        B = hmm.model_states_2_model(model).B[3:]
        assert infpy.check_is_close_2(B[0,0], B[1,3])
        assert infpy.check_is_close_2(B[0,1], B[1,2])
        assert infpy.check_is_close_2(B[0,2], B[1,1])
        assert infpy.check_is_close_2(B[0,3], B[1,0])




if __name__ == '__main__':
    unittest.main()
