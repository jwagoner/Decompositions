#!/usr/bin/env python


import numpy as np
from matplotlib import pyplot as pp
import seaborn as sns


sns.set_style('white')

nx = ny = nz = 100
n_lambda = 40

x = np.linspace(-5, 5, nx)
y = np.linspace(-5, 5, ny)
z = np.linspace(-5, 5, nz)

xv, yv, zv = np.meshgrid(x, y, z)

lambdas = np.linspace(0, 1, n_lambda)
l1v, l2v = np.meshgrid(lambdas, lambdas)

A_grid = np.zeros((n_lambda, n_lambda))

for i, lambda1 in enumerate(lambdas):
    for j, lambda2 in enumerate(lambdas):
        A_grid[i, j] = -np.log(np.sum(np.exp(-zv**2 - 1.0 * lambda1 * xv**2 - 1.0 * lambda2 * yv**2)))


dl1, dl2 = np.gradient(A_grid)

WIDTH = 0.0075
HEADLENGTH = 5
SCALE = 0.7

pp.figure(figsize=(4, 4))
pp.subplot(2, 2, 1, aspect='equal')
pp.imshow(A_grid.T, origin='bottom', extent=(0, 1, 0, 1), cmap=pp.get_cmap('viridis'))
pp.locator_params(nbins=3, tight=True)
pp.xlabel('$\lambda_1$')
pp.ylabel('$\lambda_2$')
pp.text(-0.4, 0.9, 'A', fontsize=16, fontweight='bold', transform=pp.gca().transAxes)
pp.xlim(0, 1)
pp.ylim(0, 1)

pp.subplot(2, 2, 2, aspect='equal')
pp.quiver(l1v[::4, ::4], l2v[::4, ::4], dl1[::4, ::4], dl2[::4, ::4], pivot='mid',
          width=WIDTH, headlength=HEADLENGTH, scale=SCALE)
pp.xlabel('$\lambda_1$')
pp.locator_params(nbins=3, tight=True)
pp.gca().get_yaxis().set_ticklabels([])
pp.text(-0.15, 0.9, 'B', fontsize=16, fontweight='bold', transform=pp.gca().transAxes)
pp.xlim(0, 1)
pp.ylim(0, 1)

pp.subplot(2, 2, 3, aspect='equal')
pp.quiver(l1v[::4, ::4], l2v[::4, ::4], dl1[::4, ::4], 0, pivot='mid',
          width=WIDTH, headlength=HEADLENGTH, scale=SCALE)
pp.xlabel('$\lambda_1$')
pp.locator_params(nbins=3, tight=True)
pp.text(-0.4, 0.9, 'C', fontsize=16, fontweight='bold', transform=pp.gca().transAxes)
pp.xlim(0, 1)
pp.ylim(0, 1)

pp.subplot(2, 2, 4, aspect='equal')
pp.quiver(l1v[::4, ::4], l2v[::4, ::4], 0, dl2[::4, ::4], pivot='mid',
          width=WIDTH, headlength=HEADLENGTH, scale=SCALE)
pp.xlabel('$\lambda_1$')
pp.locator_params(nbins=3, tight=True)
pp.gca().get_yaxis().set_ticklabels([])
pp.text(-0.15, 0.9, 'D', fontsize=16, fontweight='bold', transform=pp.gca().transAxes)
pp.xlim(0, 1)
pp.ylim(0, 1)

pp.tight_layout(pad=0.3)
pp.savefig('vectors.pdf')
