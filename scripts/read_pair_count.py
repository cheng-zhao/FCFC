#!/usr/bin/env python3
# read_pair_count.py: this file is part of the FCFC program.
#
# FCFC: Fast Correlation Function Calculator.
#
# Github repository:
#       https://github.com/cheng-zhao/FCFC
#
# Copyright (c) 2020 -- 2022 Cheng Zhao <zhaocheng03@gmail.com>
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

import numpy as np

def read_pair_count(ifile, verbose=False):
  '''
  Read a binary-format pair count file.
  Return:
    bin_type:   type of separation bins (0: s, 1: (s,mu), 2: (s_perp,pi))
    nsbin:      number of s or s_perp bins
    nmbin:      number of mu bins
    npbin:      number of pi bins
    s_edges:    edges of s or s_perp bins
    pi_edges:   edges of pi bins
    num1:       (weighted) number of objects from the first catalog
    num2:       (weighted) number of objects from the second catalog
    norm:       normalization factor for pair counts
    cnt_norm:   the normalized pair counts
    cnt:        the raw pair counts
  '''
  with open(ifile, 'rb') as f:
    # Read the header
    header_keys = f.readline().decode('ascii').split()
    assert len(header_keys) == 10 and \
        header_keys[0:3] == ['#','Created','by'], \
        'unrecognized pair count file'
    code_name = header_keys[3]
    assert code_name == 'FCFC_2PT' or code_name == 'FCFC_2PT_BOX', \
        'the pair count file was not created by FCFC'
    assert header_keys[6:] == [';','format','=','0'], \
        'the pair count file is not in binary format'

    # Read the pair count configurations
    size = np.fromfile(f, dtype=np.uint64, count=1)[0]
    assert size == (4 * 7) or size == (4 * 8), 'invalid pair count file'
    spec = np.fromfile(f, dtype=np.int32, count=int(size/4))
    assert np.fromfile(f, dtype=np.uint64, count=1)[0] == size, \
        'invalid pair count file'

    is_double, with_mu_one, is_box, is_cross, with_wt = spec[:5].astype('bool')
    bin_type = spec[5]                  # 0,1,2 for ISO,SMU,SPI
    assert bin_type >= 0 and bin_type <= 2, 'invalid separation bin type'
    nsbin = spec[6]                     # number of s or s_perp bins
    ntot = nsbin
    if bin_type == 1:
      nmbin = spec[7]                   # number of mu bins
      ntot *= nmbin
    else: nmbin = None
    if bin_type == 2:
      npbin = spec[7]                   # number of pi bins
      ntot *= npbin
    else: npbin = None

    # Read the box sizes
    if is_box:
      size = np.fromfile(f, dtype=np.uint64, count=1)[0]
      if is_double:
        assert size == 3 * 8, 'invalid pair count file'
        bsize = np.fromfile(f, dtype=np.float64, count=3)
      else:
        assert size == 3 * 4, 'invalid pair count file'
        bsize = np.fromfile(f, dtype=np.float32, count=3)
      assert np.fromfile(f, dtype=np.uint64, count=1)[0] == size, \
          'invalid pair count file'

    # Read the separation (or s_perp) bin edges
    size = np.fromfile(f, dtype=np.uint64, count=1)[0]
    if is_double:
      assert size == (nsbin + 1) * 8, 'invalid pair count file'
      s_edges = np.fromfile(f, dtype=np.float64, count=int(nsbin+1))
    else:
      assert size == (nsbin + 1) * 4, 'invalid pair count file'
      s_edges = np.fromfile(f, dtype=np.float32, count=int(nsbin+1))
    assert np.fromfile(f, dtype=np.uint64, count=1)[0] == size, \
        'invalid pair count file'

    # Read the pi bin edges
    if bin_type == 2:
      size = np.fromfile(f, dtype=np.uint64, count=1)[0]
      if is_double:
        assert size == (npbin + 1) * 8, 'invalid pair count file'
        pi_edges = np.fromfile(f, dtype=np.float64, count=int(npbin+1))
      else:
        assert size == (npbin + 1) * 4, 'invalid pair count file'
        pi_edges = np.fromfile(f, dtype=np.float32, count=int(npbin+1))
      assert np.fromfile(f, dtype=np.uint64, count=1)[0] == size, \
          'invalid pair count file'
    else:
      pi_edges = None

    # Read the (weighted) number of tracers
    if with_wt: num_dtype = np.float64
    else: num_dtype = np.uint64
    size = np.fromfile(f, dtype=np.uint64, count=1)[0]
    assert size == (int(is_cross) + 1) * 8, 'invalid pair count file'
    if is_cross:
      num1, num2 = np.fromfile(f, dtype=num_dtype, count=2)
    else:
      num1 = np.fromfile(f, dtype=num_dtype, count=1)[0]
      num2 = num1
    assert np.fromfile(f, dtype=np.uint64, count=1)[0] == size, \
        'invalid pair count file'

    # Read the normalization factor
    size = np.fromfile(f, dtype=np.uint64, count=1)[0]
    assert size == 8, 'invalid pair count file'
    norm = np.fromfile(f, dtype=np.float64, count=1)[0]
    assert np.fromfile(f, dtype=np.uint64, count=1)[0] == size, \
        'invalid pair count file'

    # Read the normalized pair counts
    size = np.fromfile(f, dtype=np.uint64, count=1)[0]
    assert size == ntot * 8, 'invalid pair count file'
    cnt_norm = np.fromfile(f, dtype=np.float64, count=ntot)
    assert np.fromfile(f, dtype=np.uint64, count=1)[0] == size, \
        'invalid pair count file'

    # Read the raw pair counts
    if with_wt: cnt_dtype = np.float64
    else: cnt_dtype = np.int64
    size = np.fromfile(f, dtype=np.uint64, count=1)[0]
    assert size == ntot * 8, 'invalid pair count file'
    cnt = np.fromfile(f, dtype=cnt_dtype, count=ntot)
    assert np.fromfile(f, dtype=np.uint64, count=1)[0] == size, \
        'invalid pair count file'

    if bin_type == 1:
      cnt_norm = cnt_norm.reshape([nmbin,nsbin])
      cnt = cnt.reshape([nmbin,nsbin])
    elif bin_type == 2:
      cnt_norm = cnt_norm.reshape([npbin,nsbin])
      cnt = cnt.reshape([npbin,nsbin])

    if verbose:
      print('Pair count file:', ifile)
      print('Created by:', code_name)
      print('Version:', header_keys[4], header_keys[5])
      if is_double: print('Precision: double')
      else: print('Precision: float')
      if bin_type == 1: print('Include mu = 1:', with_mu_one)
      print('Periodic boundary:', is_box)
      if is_box: print('Box sizes:', bsize)
      if is_cross: print('Pair count type: cross')
      else: print('Pair count type: auto')
      print('With weights:', with_wt)
      print('Separation bin type:',
          ['isotropic','s & mu','s_perp & pi'][bin_type])

  return bin_type, nsbin, nmbin, npbin, s_edges, pi_edges, \
      num1, num2, norm, cnt_norm, cnt

