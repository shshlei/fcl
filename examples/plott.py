#!/usr/bin/env python

######################################################################
# Software License Agreement (BSD License)
#
#  Copyright (c) 2012, Rice University
#  All rights reserved.
#
#  Redistribution and use in source and binary forms, with or without
#  modification, are permitted provided that the following conditions
#  are met:
#
#   * Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#   * Redistributions in binary form must reproduce the above
#     copyright notice, this list of conditions and the following
#     disclaimer in the documentation and/or other materials provided
#     with the distribution.
#   * Neither the name of the Rice University nor the names of its
#     contributors may be used to endorse or promote products derived
#     from this software without specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
#  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
#  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
#  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
#  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
#  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
#  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
#  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
#  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
#  POSSIBILITY OF SUCH DAMAGE.
######################################################################

# Author: Shi Shenglei

import argparse
import os

import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, sqrt

if __name__ == "__main__":
    # Create an argument parser
    parser = argparse.ArgumentParser(description='Draw time efficiency.')
    parser.add_argument('-t1', '--time1', default=None, help='Filename of time1.')
    parser.add_argument('-t2', '--time2', default=None, help='Filename of time2.')
    args = parser.parse_args()

    time1 = []
    time2 = []
    if args.time1:
        for line in open(args.time1, 'r').readlines():
            l = line.strip()
            if not l:
                continue
            time1.append([float(x) for x in l.split(' ')])
    if args.time2:
        for line in open(args.time2, 'r').readlines():
            l = line.strip()
            if not l:
                continue
            time2.append([float(x) for x in l.split(' ')])

    names = ["Box", "Ellipsoid", "Sphere", "Cylinder", "Capsule", "Cone"]
    labels= []
    percentage = []
    for i in range(len(time1)):
        name1 = names[i]
        for j in range(i, len(time1[0])):
            name2 = names[j]
            labels.append(name1 + "-" + name2)
            percentage.append(100.0 * (time1[i][j] - time2[i][j]) / time1[i][j])

    plt.style.use(['seaborn-deep', 'seaborn-paper'])
    plt.rcParams.update({'axes.grid': False})
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = prop_cycle.by_key()['color']

    fig, axs = plt.subplots()#(figsize=(3.5, 3.5))

    width = .5
    x = range(len(percentage))

    #colors = ['#1072b4', '#ef7f00', '#1f9c3a', '#d6191b', '#9763a6', '#965947', '#db79ae', '#c0c205', '#26b9ce', '#2a3377', '#7f3b71', '#9fa0a0']
    barplot = axs.bar(x, percentage, width)#, color = colors)#, sym='k+', vert=1, whis=1.5, bootstrap=1000)

    axs.set_xticks(x)
    #axs.set_xticklabels(labels)

    #axs.set_xlabel('RRT Step Size')
    axs.set_ylabel('Improvement (%)')

    #axs.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
    #axs.yaxis.set_major_locator(ticker.FixedLocator([0, 0.2, 0.4, 0.6, 0.8]))
    #axs.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.0f} s"))
    #axs.yaxis.set_minor_formatter(ticker.NullFormatter())
    axs.yaxis.set_tick_params(which='major', direction = 'in', width=1.00, length=3)
    axs.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)

    #axs.xaxis.set_major_locator(ticker.NullLocator())
    #axs.xaxis.set_minor_locator(ticker.NullLocator())
    #axs.xaxis.set_major_formatter(ticker.NullFormatter())
    axs.xaxis.set_tick_params(which='major', direction = 'in', width=1.00, length=0)

    plt.savefig('improvement.png')
    plt.show()
