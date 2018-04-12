# -*- coding: utf-8 -*-
"""
Created on Tue May  2 17:25:21 2017

:author: Ruyman Azzollini

"""


class Catalog(dict):
    """Catalog Class of vison.


    """

    def burst(self, inputs, indexes=None):
        """ """

    def insert(self, inputs, indexes=None):
        """ """

    def zip(self, f, inputs, outputs, indexes, **kwargs):
        """ """

        self.insert(outputs, indexes) = f(self.burst(inputs, indexes), **kwargs)


def proof_of_concept():
    """ """

    C = Catalog()
    C.add_column('Fluence', dimensions=['CCD', 'Q', 'Spot'])
    C.add_column('Exptime', dimensions=['CCD'])

    def get_flux(fluence, exptime):
        return fluence/exptime

    C.zip(get_flux, inputs=['Fluence', 'Exptime'], outputs=['Flux'])


if __name__ == '__main__':

    proof_of_concept()
