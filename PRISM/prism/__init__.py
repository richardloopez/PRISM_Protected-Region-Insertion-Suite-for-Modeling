#!/usr/bin/env python3
# Author: Richard Lopez Corbalan
# GitHub: github.com/richardloopez
#
# PRISM (Protected-Region Insertion Suite for Modeling)
#
# Citation:
# If you use this software in your research, please cite:
# Lopez-Corbalan, R.

"""
PRISM: Protected-Region Insertion Suite for Modeling

A professional MODELLER pipeline for high-fidelity homology modeling and loop
refinement while maintaining experimental core coordinates completely fixed.

This package prevents coordinate drift of experimental regions during optimization
and uses HETATM repulsion shields to model loops in complex environments.
"""

__version__ = '1.0.0'
__author__ = 'Richard Lopez Corbalan'
__all__ = ['config', 'controller', 'utils', 'fixed_region_utils', 'custom_models',
           'homology_modeling', 'loop_refinement']
