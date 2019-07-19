#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
# -----------------------------------------------------------------------------


class ImageSaver:
    # -------------------------------------------------------------------------
    # Constructors/Destructors
    # -------------------------------------------------------------------------
    def __init__(self, mesh, field, component, Lbar):

        self.mesh = mesh

        self.field_copy = None
        self.field = field

        self.component = component
        self.max_value = 0
        self.Lbar = Lbar

        # compute the number of nodes in one direction
        self.nb_nodes = 0
        epsilon = 1e-8
        nodes = mesh.getNodes()
        for n in range(0, mesh.getNbNodes()):
            if np.abs(nodes[n, 1]) < epsilon:
                self.nb_nodes += 1

    # -------------------------------------------------------------------------
    # Methods
    # -------------------------------------------------------------------------
    def storeStep(self):
        if self.field_copy is None:
            current_size = 0
            self.field_copy = np.zeros(self.nb_nodes)
        else:
            current_size = self.field_copy.shape[0]
            self.field_copy.resize(current_size + self.nb_nodes)

        epsilon = 1e-8
        h = self.Lbar / (self.nb_nodes-1)

        nodes = self.mesh.getNodes()
        for n in range(0, self.mesh.getNbNodes()):
            if np.abs(nodes[n, 1]) < epsilon:
                normed_x = nodes[n, 0]/h + h/10.
                index = int(normed_x)
                self.field_copy[current_size +
                                index] = self.field[n, self.component]
                if self.max_value < self.field[n, self.component]:
                    self.max_value = self.field[n, self.component]

    def getImage(self):
        width = int(self.nb_nodes)
        height = int(self.field_copy.shape[0] / self.nb_nodes)

        if np.abs(self.max_value) > 1e-8:
            for n in range(0, self.field_copy.shape[0]):
                self.field_copy[n] = 1 - self.field_copy[n] / self.max_value

        img = self.field_copy.reshape((height, width))
        return img

    def saveImage(self, filename):
        img = self.getImage()
        plt.imshow(img)
        plt.savefig(filename)
