from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from matplotlib.patches import Circle, Wedge, Polygon
from matplotlib.collections import PatchCollection

from simhom.sequences import *

plt.ion()

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig = plt.figure(1)
ax = plt.subplot(111)
plt.tight_layout()

ax.set_xlim(-1, 1)
ax.set_ylim(-1, 1)

class MakeComplex:
    def __init__(self, fig, ax):
        self.fig, self.ax = fig, ax
        self.P, self.S, self.L = [], [], []
        self.cur = []
        self.shift = False
        self.subcomplex = False
        self.cidclick = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
        self.cidpress = self.fig.canvas.mpl_connect('key_press_event', self.onpress)
        self.cidrelease = self.fig.canvas.mpl_connect('key_release_event', self.onrelease)
    def onclick(self, event):
        p = np.array([event.xdata, event.ydata])
        if not self.shift:
            s = [len(self.P)]
            self.S.append(s)
            self.P.append(p)
            if self.subcomplex:
                self.L.append(s)
            self.ax.scatter(event.xdata, event.ydata, c='black', zorder=3)
            plt.pause(0.1)
        else:
            i = min(range(len(self.P)), key=lambda i: np.linalg.norm(p - self.P[i]))
            self.cur.append(i)
    def onpress(self, event):
        if event.key == 'shift':
            self.shift = True
        if event.key == 'a':
            self.subcomplex = not self.subcomplex
    def onrelease(self, event):
        if event.key == 'shift':
            self.shift = False
            self.S.append(self.cur)
            if self.subcomplex:
                self.L.append(self.cur)
            if len(self.cur) > 1:
                for i,j in combinations(self.cur, 2):
                    ax.plot([self.P[i][0], self.P[j][0]], [self.P[i][1], self.P[j][1]], c='black', zorder=0)
            if len(self.cur) > 2:
                patches = []
                for f in combinations(self.cur, 3):
                    f = np.vstack([self.P[i] for i in f])
                    polygon = Polygon(f, True)
                    patches.append(polygon)
                p = PatchCollection(patches, alpha=0.4, color='blue', zorder=1)
                # p.set_array(np.array([0 for _ in patches]))
                ax.add_collection(p)
            if len(self.cur) > 3:
                patches = []
                for f in combinations(self.cur, 4):
                    f = np.vstack([self.P[i] for i in f])
                    polygon = Polygon(f, True)
                    patches.append(polygon)
                p = PatchCollection(patches, alpha=0.4, color='red', zorder=2)
                # p.set_array(np.array([30 for _ in patches]))
                ax.add_collection(p)
            plt.pause(0.1)
            self.cur = []
    def run(self):
        K = SimplicialComplex(list(map(Simplex, self.S)), 'K_1')
        L = SimplicialComplex(list(map(Simplex, self.L)), 'L_1')
        S = ChainComplexPairSES.simplicial_init(K, L)
        return HomologyPairLES.complex_init(S.cycles() / S.boundaries())
    def clear(self):
        self.P, self.S, self.L = [], [], []
        self.cur = []
        self.shift = False
        self.subcomplex = False
        self.ax.cla()
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        plt.pause(0.1)





if __name__ == '__main__':
    M = MakeComplex(fig, ax)
