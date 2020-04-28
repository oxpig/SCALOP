'''
 -------------------------------------------------------------------------
 Function:
 [RD,CD,order]=optics(x,k)
 -------------------------------------------------------------------------
 Aim:
 Ordering objects of a data set to obtain the clustering structure
 -------------------------------------------------------------------------
 Input:
 x - data set (m,n); m-objects, n-variables
 k - number of objects in a neighborhood of the selected object
 (minimal number of objects considered as a cluster)
 -------------------------------------------------------------------------
 Output:
 RD - vector with reachability distances (m,1)
 CD - vector with core distances (m,1)
 order - vector specifying the order of objects (1,m)
 -------------------------------------------------------------------------
 Example of use:
 x=[randn(30,2)*.4;randn(40,2)*.5+ones(40,1)*[4 4]];
 [RD,CD,order]=optics(x,4)
 -------------------------------------------------------------------------
 References:
 [1] M. Ankrest, M. Breunig, H. Kriegel, J. Sander,
 OPTICS: Ordering Points To Identify the Clustering Structure,
 available from www.dbs.informatik.uni-muenchen.de/cgi-bin/papers?query=--CO
 [2] M. Daszykowski, B. Walczak, D.L. Massart, Looking for natural
 patterns in analytical data. Part 2. Tracing local density
 with OPTICS, J. Chem. Inf. Comput. Sci. 42 (2002) 500-507
 -------------------------------------------------------------------------
 Written by Michal Daszykowski
 Department of Chemometrics, Institute of Chemistry,
 The University of Silesia
 December 2004
 http://www.chemometria.us.edu.pl


ported to python Jan, 2009 by Brian H. Clowers, Pacific Northwest National Laboratory.
Dependencies include scipy, numpy, and hcluster.
bhclowers at gmail.com
'''


import numpy as N
import pylab as P
import scipy.spatial.distance as H
from numpy import argmin
from numpy import where


def optics(x, k, distMethod = 'euclidean', p=2):
    if len(x.shape)>1:
        m,n = x.shape
    else:
        m = x.shape[0]
        n == 1

    # try:
    D = H.squareform(H.pdist(x, distMethod, p))
        # distOK = True
    # except:
    #     print "squareform or pdist error"
    #     distOK = False

    print("here")

    CD = N.zeros(m)
    RD = N.ones(m)*1E10

    for i in range(m):
        #again you can use the euclid function if you don't want hcluster
#        d = euclid(x[i],x)
#        d.sort()
#        CD[i] = d[k]

        tempInd = D[i].argsort()
        tempD = D[i][tempInd]
#        tempD.sort() #we don't use this function as it changes the reference
        CD[i] = tempD[k]#**2

    print("here")
    order = []
    seeds = N.arange(m, dtype = N.int)

    ind = 0
    while len(seeds) != 1:
#    for seed in seeds:
        ob = seeds[ind]
        seedInd = N.where(seeds != ob)
        seeds = seeds[seedInd]

        order.append(ob)
        tempX = N.ones(len(seeds))*CD[ob]
        tempD = D[ob][seeds]#[seeds]
        #you can use this function if you don't want to use hcluster
        #tempD = euclid(x[ob],x[seeds])

        temp = N.column_stack((tempX, tempD))
        mm = N.max(temp, axis = 1)
        ii = N.where(RD[seeds]>mm)[0]
        RD[seeds[ii]] = mm[ii]
        ind = N.argmin(RD[seeds])

    print("here")
    order.append(seeds[0])
    RD[0] = 0 #we set this point to 0 as it does not get overwritten
    return RD, CD, order

def optics_precomputed(D, k):

    m = D.shape[0]

    CD = N.zeros(m)
    RD = N.ones(m)*1E10

    for i in range(m):
        #again you can use the euclid function if you don't want hcluster
#        d = euclid(x[i],x)
#        d.sort()
#        CD[i] = d[k]

        tempInd = D[i].argsort()
        tempD = D[i][tempInd]
#        tempD.sort() #we don't use this function as it changes the reference
        CD[i] = tempD[k]#**2


    order = []
    seeds = N.arange(m, dtype = N.int)

    ind = 0
    while len(seeds) != 1:
#    for seed in seeds:
        ob = seeds[ind]
        seedInd = N.where(seeds != ob)
        seeds = seeds[seedInd]

        order.append(ob)
        tempX = N.ones(len(seeds))*CD[ob]
        tempD = D[ob][seeds]#[seeds]
        #you can use this function if you don't want to use hcluster
        #tempD = euclid(x[ob],x[seeds])

        temp = N.column_stack((tempX, tempD))
        mm = N.max(temp, axis = 1)
        ii = N.where(RD[seeds]>mm)[0]
        RD[seeds[ii]] = mm[ii]
        ind = N.argmin(RD[seeds])


    order.append(seeds[0])
    RD[0] = 0 #we set this point to 0 as it does not get overwritten
    return RD, CD, order

def euclid(i, x):
    """euclidean(i, x) -> euclidean distance between x and y"""
    y = N.zeros_like(x)
    y += 1
    y *= i
    if len(x) != len(y):
        raise ValueError("vectors must be same length")

    d = (x-y)**2
    return N.sqrt(N.sum(d, axis = 1))

class OpticsClusters:

    def __init__(self, Xi, ordered_RD, minPts):

        self.up_steep_areas = []

        self.down_steep_areas = []

        self.down_steep_areas_mibs = []

        self.valid_clusters = []

        self.counter = 0

        self.up_step = 0

        self.down_step = 0

        self.start_index_up_point = 0

        self.start_index_down_point = 0

        self.end_index_down_point = 0 

        self.end_index_up_point = 0

        self.Xi = Xi

        self.minPts = minPts

        self.ordered_RD = ordered_RD

        self.RD_size = len(self.ordered_RD)

        self.i = 0

        while self.i <= self.RD_size - 2:

            RD_i0 = self.ordered_RD[self.i]

            RD_i1 = self.ordered_RD[self.i+1]

            self.check_for_downsteep_start(RD_i1, RD_i0)

            self.check_for_upsteep_start(RD_i1, RD_i0)

            self.i+=1 

    def up_step_creation(self, RD_i1, RD_i0):

        check = 1

        if RD_i1 >= RD_i0:

            if RD_i1*(1-self.Xi) < RD_i0:

                self.counter += 1

            else:

                self.end_index_up_point = self.i + 1

                self.counter = 0

            if self.counter >= self.minPts:

                check = 0

        else:

            check = 0

        if self.i == self.RD_size - 2:

            check = 0 

            self.i += 1

        if check == 0:

            self.up_steep_areas.append((self.start_index_up_point, self.end_index_up_point))

            self.up_step = 0

            self.i -= 1

        else:

            self.i += 1


    def down_step_creation(self, RD_i1, RD_i0):

        check = 1

        if RD_i1 <= RD_i0:

            if RD_i1 > RD_i0*(1-self.Xi):

                self.counter += 1

            else:

                self.end_index_down_point = self.i + 1

                self.counter = 0

            if self.counter >= self.minPts:

                check = 0

        else:

            check = 0

        if self.i == self.RD_size - 2:

            check = 0

            self.i += 1

        if check == 0:

            self.down_steep_areas.append((self.start_index_down_point, self.end_index_down_point))

            self.down_step = 0

            self.i -= 1

        else:

            self.i += 1


    def check_for_upsteep_start(self, RD_i1, RD_i0):

        if RD_i1*(1-self.Xi) >= RD_i0 and self.up_step == 0:

            self.start_index_up_point = self.i

            self.end_index_up_point = self.i + 1

            self.up_step = 1

            self.counter = 0

            while self.up_step == 1:

                RD_i0 = self.ordered_RD[self.i]

                RD_i1 = self.ordered_RD[self.i+1]

                self.up_step_creation(RD_i1, RD_i0)

            if self.down_steep_areas != []:

                self.find_valid_clusters()


    def find_valid_clusters(self):

        current_up_steep_area = self.up_steep_areas[-1]

        for steep_down_area in self.down_steep_areas:

            if self.ordered_RD[steep_down_area[0]]*(1 - self.Xi) >= self.ordered_RD[current_up_steep_area[1]]:

                # print str(steep_down_area[0]) + ' ' + str(steep_down_area[1])

                # print str(current_up_steep_area[0]) + ' ' + str(current_up_steep_area[1])

                step_down_values = self.ordered_RD[steep_down_area[0]:steep_down_area[1]]

                distance_array = step_down_values - self.ordered_RD[current_up_steep_area[1]]

                # print distance_array

                start_index = steep_down_area[0] + where(distance_array >= 0)[0][-1]

                # print start_index

                #start_index = steep_down_area[0] + argmin(distance_array)

            else:

                start_index = steep_down_area[0]

            if self.ordered_RD[current_up_steep_area[1]]*(1 - self.Xi) >= self.ordered_RD[steep_down_area[0]]:

                step_up_values = self.ordered_RD[current_up_steep_area[0]:current_up_steep_area[1]]

                distance_array = step_up_values - self.ordered_RD[steep_down_area[0]]

                end_index = current_up_steep_area[1] - where(distance_array >= 0)[0][0]

                #end_index = current_up_steep_area[1] - (len(distance_array) - 1 + argmin(distance_array))

            else:

                end_index = current_up_steep_area[1]


            cluster_valid = self.check_cluster_validity(start_index,end_index)

            if cluster_valid:

                self.valid_clusters.append((start_index, end_index))


    def check_cluster_validity(self, start_index, end_index):

        cluster_valid = 1

        if end_index - start_index < self.minPts:

            print("here 1")

            cluster_valid = 0

        elif max(self.ordered_RD[start_index + 1 :end_index - 1]) > self.ordered_RD[start_index] * (1 - self.Xi):

            print("here 2")

            print(str(max(self.ordered_RD[start_index + 1 :end_index - 1])) + ' ' + str(start_index) +  ' ' + str(end_index) + ' ' + str(self.ordered_RD[start_index]))

            cluster_valid = 0

        elif max(self.ordered_RD[start_index + 1:end_index - 1]) > self.ordered_RD[end_index] * (1 - self.Xi):#

            print("here 3")

            cluster_valid = 0

        return cluster_valid




    def check_for_downsteep_start(self, RD_i1, RD_i0):

        if RD_i1 <= RD_i0*(1-self.Xi) and self.down_step == 0:

            self.start_index_down_point = self.i

            self.end_index_down_point = self.i + 1

            self.down_step = 1

            self.counter = 0

            while self.down_step == 1:

                RD_i0 = self.ordered_RD[self.i]

                RD_i1 = self.ordered_RD[self.i+1]

                self.down_step_creation(RD_i1, RD_i0)


    def find_clusters_hierarchy(self):

        contained_in_clusters = {}

        for i in range(0,len(self.valid_clusters)):

            clusteri = self.valid_clusters[i]

            contained_in_clusters[clusteri] = []

            for j in range(0,len(self.valid_clusters)):

                clusterj = self.valid_clusters[j]

                if i!=j and clusterj[0] <= clusteri[0] and clusterj[1] >= clusteri[1]:

                    contained_in_clusters[clusteri].append(clusterj)

            contained_in_clusters[clusteri] = set(contained_in_clusters[clusteri])


        return contained_in_clusters

def find_steep_regions(Xi, ordered_RD, minPts):

    Clusters = OpticsClusters()

    Clusters.Xi = Xi

    Clusters.minPts = minPts

    for i in range(len(ordered_RD) - 1):

        RD_i0 = ordered_RD[i]

        RD_i1 = ordered_RD[i+1]

        if Clusters.up_step == 1:

            Clusters.up_step_creation(RD_i1, RD_i0, i)

            continue

        Clusters.check_for_upsteep_start(self, RD_i1, RD_i0, i)

    return OpticsClusters

        





if __name__ == "__main__":

    pass