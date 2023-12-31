#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  test_propersubtract.py
#
#  Copyright 2020 QuatroPe
#
# This file is part of ProperImage (https://github.com/quatrope/ProperImage)
# License: BSD-3-Clause
# Full Text: https://github.com/quatrope/ProperImage/blob/master/LICENSE.txt
#

"""
test_propersubtract module from ProperImage,
for subtracting astronomical images.

Written by Bruno SANCHEZ

PhD of Astromoy - UNC
bruno@oac.unc.edu.ar

Instituto de Astronomia Teorica y Experimental (IATE) UNC
Cordoba - Argentina

Of 301
"""

import os
import shutil
import tempfile
import unittest

from astropy.io import fits

import numpy as np
from numpy.random import default_rng

from properimage import SingleImage, simtools, subtract

random = default_rng(seed=42)


class PropersubtractBase(object):
    def setUp(self):
        print("setting up")
        self.tempdir = tempfile.mkdtemp()
        self.paths = []
        # mock data
        psf = simtools.Psf(11, 2.5, 3.0)

        self.mock_image_data = random.random((256, 256)) * 10.0
        for i in range(50):
            x = random.integers(7, 220)
            y = random.integers(7, 120)
            # print x, y
            self.mock_image_data[x : x + 11, y : y + 11] += (
                psf * float(i + 1) * 2000.0
            )

        psf = simtools.Psf(11, 3.0, 1.9)
        for i in range(50):
            x = random.integers(7, 220)
            y = random.integers(122, 220)
            # print x, y
            self.mock_image_data[x : x + 11, y : y + 11] += (
                psf * float(i + 1) * 2000.0
            )

        # generate 2 images to subtract

        # a reference
        refimage_data = (
            self.mock_image_data * 10.0
            + random.random((256, 256)) * 50.0
            + 350
        )

        refimage_data[123, 123] = np.nan

        # a fits file
        refmockfits_path = os.path.join(self.tempdir, "refmockfits.fits")
        fits.writeto(refmockfits_path, refimage_data, overwrite=True)
        self.paths.append(refmockfits_path)

        # a new image
        rand_pix = random.random((256, 256)) * 50.0 + 750
        newimage_data = self.mock_image_data * 3.0 + rand_pix

        newimage_data[123, 123] = np.nan

        # a fits file
        newmockfits_path = os.path.join(self.tempdir, "newmockfits.fits")
        fits.writeto(newmockfits_path, newimage_data, overwrite=True)
        self.paths.append(newmockfits_path)

    def tearDown(self):
        if os.path.isdir(self.tempdir):
            shutil.rmtree(self.tempdir)
        try:
            self.si_ref._clean()
            self.si_new._clean()
        except OSError:
            pass


class TestSubtract(PropersubtractBase, unittest.TestCase):
    def setUp(self):
        super(TestSubtract, self).setUp()
        self.si_ref = SingleImage(self.paths[0])
        self.si_new = SingleImage(self.paths[1])

    def testSubtractNoBeta(self):
        D, P, S_corr, mix_mask = subtract(self.si_ref, self.si_new, beta=False)
        self.assertIsInstance(D, np.ndarray)
        self.assertIsInstance(P, np.ndarray)
        self.assertIsInstance(S_corr, np.ndarray)
        self.assertIsInstance(mix_mask, np.ndarray)

    def testSubtractBetaIterative(self):
        D, P, S_corr, mix_mask = subtract(
            self.si_ref, self.si_new, beta=True, iterative=True, shift=False
        )
        self.assertIsInstance(D, np.ndarray)
        self.assertIsInstance(P, np.ndarray)
        self.assertIsInstance(S_corr, np.ndarray)
        self.assertIsInstance(mix_mask, np.ndarray)

    def testSubtractBetaShift(self):
        D, P, S_corr, mix_mask = subtract(
            self.si_ref, self.si_new, beta=True, iterative=False, shift=True
        )
        self.assertIsInstance(D, np.ndarray)
        self.assertIsInstance(P, np.ndarray)
        self.assertIsInstance(S_corr, np.ndarray)
        self.assertIsInstance(mix_mask, np.ndarray)

    def testSubtractOnlyBeta(self):
        D, P, S_corr, mix_mask = subtract(
            self.si_ref, self.si_new, beta=True, iterative=False, shift=False
        )
        self.assertIsInstance(D, np.ndarray)
        self.assertIsInstance(P, np.ndarray)
        self.assertIsInstance(S_corr, np.ndarray)
        self.assertIsInstance(mix_mask, np.ndarray)

    def testSubtractOnlyShift(self):
        D, P, S_corr, mix_mask = subtract(
            self.si_ref, self.si_new, beta=False, iterative=False, shift=True
        )
        self.assertIsInstance(D, np.ndarray)
        self.assertIsInstance(P, np.ndarray)
        self.assertIsInstance(S_corr, np.ndarray)
        self.assertIsInstance(mix_mask, np.ndarray)

    def testSubtractNoFitPSF(self):
        D, P, S_corr, mix_mask = subtract(
            self.si_ref,
            self.si_new,
            beta=False,
            iterative=False,
            shift=False,
            fitted_psf=False,
        )
        self.assertIsInstance(D, np.ndarray)
        self.assertIsInstance(P, np.ndarray)
        self.assertIsInstance(S_corr, np.ndarray)
        self.assertIsInstance(mix_mask, np.ndarray)

    def tearDown(self):
        self.si_ref._clean()
        self.si_new._clean()
