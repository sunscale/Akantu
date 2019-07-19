#!/usr/bin/env python3

import numpy as np


class Mazars:
    def __init__(self, **kwargs):
        self.K0 = kwargs.pop("K0", 1e-4)
        self.At = kwargs.pop("At", 1.0)
        self.Bt = kwargs.pop("Bt", 5e3)
        self.Ac = kwargs.pop("Ac", 0.8)
        self.Bc = kwargs.pop("Bc", 1391.3)

        self.E = kwargs.pop("E", 25e9)
        self.nu = kwargs.pop("nu", 0.2)

        self.dam = 0

        self.Gf = 0
        self.ε_p = 0
        self.σ_p = 0

    def compute_step(self, ε, σ, dam, trace):
        dam_t = 0
        # dam_c = 0

        if trace:
            import pdb
            pdb.set_trace()

        σ = self.E * ε
        if ε > self.K0:
            dam_t = 1 - self.K0*(1 - self.At)/ε - \
                self.At * np.exp(-self.Bt*(ε - self.K0))
            # dam_c = 1 - self.K0*(1 - self.Ac)/ε - \
            #     self.Ac * np.exp(-self.Bc*(ε - self.K0))

        dam = max(dam, dam_t)
        dam = min(dam, 1)

        σ = (1 - dam) * σ
        return σ, dam

    # def compute(self, **kwargs):
    #     epsilons = np.array(kwargs['epsilons'], copy=False)
    #     sigmas = np.array(kwargs['sigmas'], copy=False)
    #     damages = np.array(kwargs['damages'], copy=False)

    #     for t, ε in enumerate(epsilons):
    #         σ = 0.
    #         dam = 0.
    #         self.compute_step(ε, σ, dam)

    #         self.Gf = self.Gf + (σ + σ_p) * (ε - ε_p) / 2.
    #         self.σ_p = σ
    #         self.ε_p = ε

    #         sigmas[t] = σ
    #         damages[t] = self.dam
    #     return self.Gf
