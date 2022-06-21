#!/usr/bin/env python3

################################################################################
#  plot_fcfc_logo.py: this file is for creating the logo of FCFC
#  FCFC: Fast Correlation Function Calculator.
#  Github repository: https://github.com/cheng-zhao/FCFC
#  Copyright (c) 2020 -- 2021 Cheng Zhao <zhaocheng03@gmail.com>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
################################################################################

import numpy as np

class fcfc_logo:
  def __init__(self, figwidth=1100, outradius=150, linewidth=50, angle=54, \
      margin=0.02, font='Arial Black', fontsize=20, textoffx=0, textoffy=0):
    # Check parameters.
    if linewidth >= outradius:
      raise ValueError('linewidth must be smaller than outradius')
    # Compute the size of the logo.
    self.angle = angle * np.pi / 180
    width = outradius * (1 / np.sin(self.angle) + np.cos(self.angle)) - \
            linewidth * 0.5 * (1 / np.sin(self.angle) - np.sin(self.angle))
    width = np.ceil(2 * width * (1 + 2 * margin) / 10) * 10
    height = np.ceil(2 * outradius * (1 + 2 * margin) / 10) * 10

    # Rescale the logo given the figure size.
    self.width = round(figwidth)
    rescale = self.width / width
    self.height = round(height * rescale)

    # Computed specifications.
    self.cenx = int(self.width / 2)
    self.ceny = int(self.height / 2)
    self.outradius = outradius * rescale
    self.linewidth = linewidth * rescale
    self.inradius = self.outradius - self.linewidth
    self.font = font
    self.fontsize = fontsize
    self.textoffx = textoffx * rescale
    self.textoffy = textoffy * rescale

  def arc_ctrl_pts(self, eta1, eta2, radius):
    '''Compute the end and control points of cubic Bezier curves,
       for approximating an arc given the angles of the endpoints.
       ref: http://www.spaceroots.org/documents/ellipse/elliptical-arc.pdf
    '''
    if eta1 >= eta2:
      raise ValueError('the first angle must be smaller than the second one')
    delta = eta2 - eta1
    alpha = np.sin(delta) * (np.sqrt(4 + 3 * np.tan(delta/2)**2) - 1) / 3
    p1x = self.cenx + radius * np.cos(eta1)
    p1y = self.ceny + radius * np.sin(eta1)
    p2x = self.cenx + radius * np.cos(eta2)
    p2y = self.ceny + radius * np.sin(eta2)
    q1x = p1x - alpha * radius * np.sin(eta1)
    q1y = p1y + alpha * radius * np.cos(eta1)
    q2x = p2x + alpha * radius * np.sin(eta2)
    q2y = p2y - alpha * radius * np.cos(eta2)
    return ((p1x,p1y),(q1x,q1y),(q2x,q2y),(p2x,p2y))

  def get_path(self, rotate=0):
    '''Draw the path of half of the logo.'''
    from matplotlib.path import Path
    codes = []
    verts = []

    # First part of the outer circle.
    eta0 = self.angle + np.arcsin(self.linewidth*0.5/self.outradius) + rotate
    eta2 = self.angle + np.pi * 0.5 + rotate
    eta1 = (eta0 + eta2) * 0.5    # split the arc into two parts
    p0x = self.cenx + self.outradius * np.cos(eta0)
    p0y = self.ceny + self.outradius * np.sin(eta0)
    codes.extend((Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4))
    verts.extend(self.arc_ctrl_pts(eta0, eta1, self.outradius))
    codes.extend((Path.CURVE4, Path.CURVE4, Path.CURVE4))
    verts.extend(self.arc_ctrl_pts(eta1, eta2, self.outradius)[1:])

    # Endpoints of the straight line (rectangle).
    l0x0 = self.outradius * (1 / np.sin(self.angle) + np.cos(self.angle)) - \
           self.linewidth * 0.5 * (1 / np.sin(self.angle) - np.sin(self.angle))
    l0y0 = self.outradius * np.sin(self.angle) - \
           self.linewidth * 0.5 * np.cos(self.angle)
    l1x0 = self.outradius * (1 / np.sin(self.angle) + np.cos(self.angle)) - \
           self.linewidth * 0.5 * (1 / np.sin(self.angle) + np.sin(self.angle))
    l1y0 = self.outradius * np.sin(self.angle) + \
           self.linewidth * 0.5 * np.cos(self.angle)
    l0x = self.cenx - l0x0 * np.cos(rotate) + l0y0 * np.sin(rotate)
    l0y = self.ceny - l0x0 * np.sin(rotate) - l0y0 * np.cos(rotate)
    l1x = self.cenx - l1x0 * np.cos(rotate) + l1y0 * np.sin(rotate)
    l1y = self.ceny - l1x0 * np.sin(rotate) - l1y0 * np.cos(rotate)
    codes.extend((Path.LINETO, Path.LINETO))
    verts.extend(((l0x,l0y), (l1x,l1y)))

    # Intersection point of the outer circle and the inner straight line
    delta = self.outradius - self.linewidth
    gamma = np.sqrt(self.outradius**2 - delta**2)
    ix = -delta * np.sin(self.angle) - gamma * np.cos(self.angle)
    iy = delta * np.cos(self.angle) - gamma * np.sin(self.angle)

    # Second part of the outer circle
    eta0 = np.arctan2(iy, ix) + rotate
    eta1 = self.angle + np.pi - np.arcsin(self.linewidth*0.5/self.outradius) + \
           rotate
    codes.extend((Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4))
    verts.extend(self.arc_ctrl_pts(eta0, eta1, self.outradius))

    # Inner circle.
    eta0 = self.angle + np.pi - np.arcsin(self.linewidth*0.5/self.inradius) + \
           rotate
    eta3 = self.angle + np.arcsin(self.linewidth*0.5/self.inradius) + rotate
    eta1 = eta0 * 2 / 3 + eta3 / 3      # split the arc into three parts
    eta2 = eta0 / 3 + eta3 * 2 / 3
    codes.extend((Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4))
    verts.extend(self.arc_ctrl_pts(eta1, eta0, self.inradius)[::-1])
    codes.extend((Path.CURVE4, Path.CURVE4, Path.CURVE4))
    verts.extend(self.arc_ctrl_pts(eta2, eta1, self.inradius)[:-1][::-1])
    codes.extend((Path.CURVE4, Path.CURVE4, Path.CURVE4))
    verts.extend(self.arc_ctrl_pts(eta3, eta2, self.inradius)[:-1][::-1])

    # Close the path.
    codes.append(Path.CLOSEPOLY)
    verts.append((p0x, p0y))
    return Path(verts, codes)

  def text_pos(self, rotate=0):
    '''Compute the location of texts.'''
    tx0 = self.outradius * (1 / np.sin(self.angle) + np.cos(self.angle)) - \
          self.textoffx * np.cos(self.angle) - \
          self.textoffy * np.sin(self.angle) - \
          self.linewidth * 0.5 * (1 / np.sin(self.angle) - np.sin(self.angle))
    ty0 = self.outradius * np.sin(self.angle) - \
          self.textoffx * np.sin(self.angle) - \
          self.textoffy * np.cos(self.angle) - \
          self.linewidth * 0.5 * np.cos(self.angle)
    tx = self.cenx - tx0 * np.cos(rotate) + ty0 * np.sin(rotate)
    ty = self.ceny - tx0 * np.sin(rotate) - ty0 * np.cos(rotate)
    return (tx, ty)

  def export(self, filename, dpi=300):
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    # Create the figure.
    fig = plt.figure(figsize=(float(self.width)/dpi,float(self.height)/dpi), \
          dpi=dpi, frameon=False)
    ax = plt.axes([0,0,1,1])
    ax.set_axis_off()
    ax.set_xlim(0, self.width)
    ax.set_ylim(0, self.height)

    # Draw and fill the patches.
    path = self.get_path(rotate=0)
    patch = mpatches.PathPatch(path, facecolor='k', edgecolor='none')
    ax.add_patch(patch)
    path = self.get_path(rotate=np.pi)
    patch = mpatches.PathPatch(path, facecolor='k', edgecolor='none')
    ax.add_patch(patch)

    # Draw the texts.
    tx, ty = self.text_pos(rotate=0)
    ax.text(tx, ty, 'FCFC', fontfamily=self.font, \
        fontsize=self.fontsize, rotation=(self.angle+np.pi)*180/np.pi, \
        ha='center', va='center', color='w')
    tx, ty = self.text_pos(rotate=np.pi)
    ax.text(tx, ty, 'FCFC', fontfamily=self.font, \
        fontsize=self.fontsize, rotation=(self.angle)*180/np.pi, \
        ha='center', va='center', color='w')

    plt.savefig(filename)


logo = fcfc_logo(textoffx=24, textoffy=38)
logo.export('FCFC_logo.svg')
