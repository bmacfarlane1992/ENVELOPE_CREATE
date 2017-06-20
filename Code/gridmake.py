'''

gridmake.py

Routine to take SPH distribution of particles, and construct grid based on
particle location limits for world size. Method grids based on main.py variable
grid_method, for either regular, or octree grid (dev.).

Code outputs to amr_grid.inp as per required RADMC-3D IO, and returns grid
to main.py to build density grid.

Author: Benjamin MacFarlane
Date: 12/05/2017
Contact: bmacfarlane@uclan.ac.uk

'''
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - VARIABLE DEFINITIONS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
nbins = 64     # For 'reg' grid_make selection, how many bins in each dimension
#
tree_params = ["nodes",5]    # Tree restriction parameters - either ["nodes", max_value], or ["depth", max_value"]
#
leafbranch_test = False
output_test = False         # Choose whether or not to output octree checks (may be large terminal output...)
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MODULE IMPORTS - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
import time
import sys
try:
    import numpy as np
except ImportError:
    np = None
#
#
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
        # # # - - - MAIN PROGRAM - - - # # #
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
#
#
def regular(arch_dir, pos):
#
    # Define world size and bin size based on nbins variable
#
    binlims = [-1.005*np.amax(pos), 1.005*np.amax(pos)]
    bin_it = (binlims[1] - binlims[0]) / float(nbins)
#
    # Write to amr_grid.inp for RADMC-3D the locations of bin edges in all
    # dimensions. As grid is cubic, r_it = x_it = y_it = z_it
#
    f = open(arch_dir+'amr_grid.inp','w')
    f.write('1 \n')
    f.write('0 \n')
    f.write('0 \n')
    f.write('0 \n')
    f.write('1 1 1 \n')
    f.write(str(nbins)+' '+str(nbins)+' '+str(nbins)+'\n')
    for i in range(0,nbins+1):
        x = ( binlims[0] + (i * bin_it) )
        f.write(str(x)+'\n')
    for i in range(0,nbins+1):
        y = ( binlims[0] + (i*bin_it) )
        f.write(str(y)+'\n')
    for i in range(0,nbins+1):
        z = ( binlims[0] + (i*bin_it) )
        f.write(str(z)+'\n')
    f.close()
#
    # Set up a 3D float array for all bins in cubic grid to be used in
    # densgrid.py
#
    grid =  np.array(  \
       [ [ [[0.] for i in range(nbins)] for j in range(nbins) ] for k in range(nbins) ] \
       )
#
    return nbins, bin_it, binlims, grid

### ------------------------------------------------------------------------ ###

def octree(arch_dir, pos, rho, h, r_dust):

    sys.setrecursionlimit(2500)

    class OctNode(object):
        """
        New Octnode Class, can be appended to as well i think
        """
        def __init__(self, position, size, depth, data):
            """
            OctNode Cubes have a position and size
            position is related to, but not the same as the objects the node contains.
            Branches (or children) follow a predictable pattern to make accesses simple.
            Here, - means less than 'origin' in that dimension, + means greater than.
            branch: 0 1 2 3 4 5 6 7
            x:      - - - - + + + +
            y:      - - + + - - + +
            z:      - + - + - + - +
            """
            self.position = position
            self.size = size
            self.depth = depth

        ## All OctNodes will be leaf nodes at first
        ## Then subdivided later as more objects get added
            self.isLeafNode = True

        ## store our object, typically this will be one, but maybe more
            self.data = data

        ## might as well give it some emtpy branches while we are here.
            self.branches = [None, None, None, None, None, None, None, None]

            half = size / 2

        ## The cube's bounding coordinates
            self.lower = (position[0] - half, position[1] - half, position[2] - half)
            self.upper = (position[0] + half, position[1] + half, position[2] + half)

        def __str__(self):
            data_str = u", ".join((str(x) for x in self.data))
            return u"position: {0}, size: {1}, depth: {2} leaf: {3}, data: {4}".format(
                self.position, self.size, self.depth, self.isLeafNode, data_str
            )


    class Octree(object):
        """
        The octree itself, which is capable of adding and searching for nodes.
        """
        def __init__(self, worldSize, origin=(0, 0, 0), max_type="nodes", max_value=10):
            """
            Init the world bounding root cube
            all world geometry is inside this
            it will first be created as a leaf node (ie, without branches)
            this is because it has no objects, which is less than MAX_OBJECTS_PER_CUBE
            if we insert more objects into it than MAX_OBJECTS_PER_CUBE, then it will subdivide itself.
            """
            self.root = OctNode(origin, worldSize, 0, [])
            self.worldSize = worldSize
            self.limit_nodes = (max_type=="nodes")
            self.limit = max_value

        @staticmethod
        def CreateNode(position, size, objects):
            """This creates the actual OctNode itself."""
            return OctNode(position, size, objects)

        def insertNode(self, position, objData=None):
            """
            Add the given object to the octree if possible
            Parameters
            ----------
            position : array_like with 3 elements
                The spatial location for the object
                objData : optional
                The data to store at this position. By default stores the position.
                If the object does not have a position attribute, the object
                itself is assumed to be the position.
            Returns
                -------
                node : OctNode or None
                The node in which the data is stored or None if outside the
                octree's boundary volume.
            """
            if np:
                if np.any(position < self.root.lower):
                    return None
                if np.any(position > self.root.upper):
                    return None
            else:
                if position < self.root.lower:
                    return None
                if position > self.root.upper:
                    return None

            if objData is None:
                objData = position

            return self.__insertNode(self.root, self.root.size, self.root, position, objData)

        def __insertNode(self, root, size, parent, position, objData):
            """Private version of insertNode() that is called recursively"""
            if root is None:
            # we're inserting a single object, so if we reach an empty node, insert it here
            # Our new node will be a leaf with one object, our object
            # More may be added later, or the node maybe subdivided if too many are added
            # Find the Real Geometric centre point of our new node:
            # Found from the position of the parent node supplied in the arguments
                pos = parent.position

            ## offset is halfway across the size allocated for this node
                offset = size / 2

            ## find out which direction we're heading in
                branch = self.__findBranch(parent, position)

            ## new center = parent position + (branch direction * offset)
                newCenter = (0, 0, 0)

                if branch == 0:
                    newCenter = (pos[0] - offset, pos[1] - offset, pos[2] - offset )
                elif branch == 1:
                    newCenter = (pos[0] - offset, pos[1] - offset, pos[2] + offset )
                elif branch == 2:
                    newCenter = (pos[0] - offset, pos[1] + offset, pos[2] - offset )
                elif branch == 3:
                    newCenter = (pos[0] - offset, pos[1] + offset, pos[2] + offset )
                elif branch == 4:
                    newCenter = (pos[0] + offset, pos[1] - offset, pos[2] - offset )
                elif branch == 5:
                    newCenter = (pos[0] + offset, pos[1] - offset, pos[2] + offset )
                elif branch == 6:
                    newCenter = (pos[0] + offset, pos[1] + offset, pos[2] - offset )
                elif branch == 7:
                    newCenter = (pos[0] + offset, pos[1] + offset, pos[2] + offset )

            # Now we know the centre point of the new node
            # we already know the size as supplied by the parent node
            # So create a new node at this position in the tree
            # print "Adding Node of size: " + str(size / 2) + " at " + str(newCenter)
                return OctNode(newCenter, size, parent.depth + 1, [objData])

        #else: are we not at our position, but not at a leaf node either
            elif (
                not root.isLeafNode
                and
                (
                    (np and np.any(root.position != position))
                    or
                    (root.position != position)
                )
            ):

            # we're in an octNode still, we need to traverse further
                branch = self.__findBranch(root, position)
            # Find the new scale we working with
                newSize = root.size / 2
            # Perform the same operation on the appropriate branch recursively
                root.branches[branch] = self.__insertNode(root.branches[branch], newSize, root, position, objData)

        # else, is this node a leaf node with objects already in it?
            elif root.isLeafNode:
            # We've reached a leaf node. This has no branches yet, but does hold
            # some objects, at the moment, this has to be less objects than MAX_OBJECTS_PER_CUBE
            # otherwise this would not be a leafNode (elementary my dear watson).
            # if we add the node to this branch will we be over the limit?
                if (
                    (self.limit_nodes and len(root.data) < self.limit)
                    or
                    (not self.limit_nodes and root.depth >= self.limit)
                ):
                # No? then Add to the Node's list of objects and we're done
                    root.data.append(objData)
                #return root
                else:
                # Adding this object to this leaf takes us over the limit
                # So we have to subdivide the leaf and redistribute the objects
                # on the new children.
                # Add the new object to pre-existing list
                    root.data.append(objData)
                # copy the list
                    objList = root.data
                # Clear this node's data
                    root.data = None
                # It is not a leaf node anymore
                    root.isLeafNode = False
                # Calculate the size of the new children
                    newSize = root.size / 2
                # distribute the objects on the new tree
                # print "Subdividing Node sized at: " + str(root.size) + " at " + str(root.position)
                    for ob in objList:
                    # Use the position attribute of the object if possible
                        if hasattr(ob, "position"):
                            pos = ob.position
                        else:
                            pos = ob
                            branch = self.__findBranch(root, pos)
                            root.branches[branch] = self.__insertNode(root.branches[branch], newSize, root, pos, ob)
            return root

        def findPosition(self, position):
            """
            Basic lookup that finds the leaf node containing the specified position
            Returns the child objects of the leaf, or None if the leaf is empty or none
            """
            if np:
                if np.any(position < self.root.lower):
                    return None
                if np.any(position > self.root.upper):
                    return None
                else:
                    if position < self.root.lower:
                        return None
                    if position > self.root.upper:
                        return None
                return self.__findPosition(self.root, position)

        @staticmethod
        def __findPosition(node, position, count=0, branch=0):
            """Private version of findPosition """
            if node.isLeafNode:
            #print("The position is", position, " data is", node.data)
                return node.data
            branch = Octree.__findBranch(node, position)
            child = node.branches[branch]
            if child is None:
                return None
            return Octree.__findPosition(child, position, count + 1, branch)

        @staticmethod
        def __findBranch(root, position):
            """
            helper function
            returns an index corresponding to a branch
            pointing in the direction we want to go
            """
            index = 0
            if (position[0] >= root.position[0]):
                index |= 4
            if (position[1] >= root.position[1]):
                index |= 2
            if (position[2] >= root.position[2]):
                index |= 1
            return index

        def iterateDepthFirst(self):
            """Iterate through the octree depth-first"""
            gen = self.__iterateDepthFirst(self.root)
            for n in gen:
                yield n

        @staticmethod
        def __iterateDepthFirst(root):
            """Private (static) version of iterateDepthFirst"""

            for branch in root.branches:
                if branch is None:
                    continue
                for n in Octree.__iterateDepthFirst(branch):
                    yield n
                if branch.isLeafNode:
                    yield branch

#
#

        def leaf_or_branch(self):
            """Output list of leaf/branch as per amr_grid.inp requirements"""
            gen = self.__leaf_or_branch(self.root)
            for n in gen:
                yield n

        @staticmethod
        def __leaf_or_branch(root):
            """Private (static) version of leaf_or_branch)"""
            for branch in root.branches:
                for n in Octree.__iterateDepthFirst(branch):
                    yield n.branches, n.size

### ------------------------------------------------------------------------ ###
    ### Create SPH Particle Object, then construct Octree with world size ###
    ### defined by limits of envelope                                     ###
### ------------------------------------------------------------------------ ###

    plims = [-np.amax(pos), np.amax(pos)]

    class SPHpart(object):
        """Object class to contain SPH particles"""
        def __init__(self, name, position, rho, h):
            self.name = name
            self.position = position
            self.rho = rho
            self.h = h

        def __str__(self):
            return "name: {0} position: {1} rho: {2} h: {3}".format(self.name, self.position, self.rho, self.h)

    # Number of objects we intend to add.
    NUM_SPHparts = int(len(pos[0]))

    # Size that the octree covers
    WORLD_SIZE = 1.1*(plims[1]-plims[0])

    #ORIGIN = (WORLD_SIZE, WORLD_SIZE, WORLD_SIZE)
    ORIGIN = (0, 0, 0)

    # create random test objects
    SPHparts = []
    for x in range(NUM_SPHparts):
        the_name = str(x)
#        the_name = "Part__" + str(x)
        the_pos = (
            pos[0][x],
            pos[1][x],
            pos[2][x]
        )
        the_rho = rho[x]
        the_h = h[x]
        SPHparts.append(SPHpart(the_name, the_pos, the_rho, the_h))

        # Create a new octree, size of world
    myTree = Octree(
        WORLD_SIZE,
        ORIGIN,
        max_type=tree_params[0],
        max_value=tree_params[1]
    )

        # Insert some random objects and time it
    Start = time.time()
    for SPHpart in SPHparts:
        myTree.insertNode(SPHpart.position, SPHpart)
    End = time.time() - Start

### ------------------------------------------------------------------------ ###
    ### Print results, outputting tests in length units of r_dust  ###
### ------------------------------------------------------------------------ ###

    if leafbranch_test is True:

        for i, x in enumerate(myTree.leaf_or_branch()):
            print x


    if output_test is True:
#
        print NUM_SPHparts, "Node Tree Generated in", End, "Seconds"
        print "Tree centered at", ORIGIN, "with size", WORLD_SIZE/r_dust
        if myTree.limit_nodes:
            print "Tree Leaves contain a maximum of", myTree.limit, "objects each."
        else:
            print "Tree has a maximum depth of", myTree.limit

        for i, x in enumerate(myTree.iterateDepthFirst()):
            print "Node: ",i,", position: ", x.position[0]/r_dust,", size: ", x.size/r_dust
            for j, y in enumerate(x.data):
                print "Particle", y.name,", position(x):", y.position[0]/r_dust, "density:", y.rho, "h:", y.h
            print '\n'


### ------------------------------------------------------------------------ ###
    ### Use results to ouput grid data into amr_grid.inp for RADMC-3D  ###
### ------------------------------------------------------------------------ ###

    octbase =  max(x.size for x in myTree.iterateDepthFirst())
    n_octbase =  int(WORLD_SIZE / octbase)
    base0 = 0. - 0.5*WORLD_SIZE
    levelmax = int( max(x.depth for x in myTree.iterateDepthFirst()) )
    nleafsmax = 0 ; nbranchmax = 0
    for i, x in enumerate(myTree.iterateDepthFirst()):
        nleafsmax += 1 ; nbranchmax += 7
#
    f = open(arch_dir+'amr_grid.inp','w')
    f.write('1 \n')
    f.write('1 \n')
    f.write('0 \n')
    f.write('0 \n')
    f.write('1 1 1 \n')
    f.write(str(n_octbase)+' '+str(n_octbase)+' '+str(n_octbase)+'\n')
    f.write(str(levelmax)+' '+str(nleafsmax)+' '+str(nbranchmax)+'\n')
    for i in range(0,n_octbase+1):
        x = ( base0 + (i * octbase) )
        f.write(str(x)+'\t')
    f.write('\n')
    for i in range(0,n_octbase+1):
        y = ( base0 + (i * octbase) )
        f.write(str(y)+'\t')
    f.write('\n')
    for i in range(0,n_octbase+1):
        z = ( base0 + (i * octbase) )
        f.write(str(z)+'\t')
    f.write('\n')
#
    # Loop over potential node locations in an ordered manner, outputting whether
    # leaf or branch in RADMC-3D output requirements

    for i in range(0, n_octbase):
        for j in range(0, n_octbase):
            for k in range(0, n_octbase):

                print i, j, k

                pos_find = (
                base0 + 0.5*octbase + (octbase * float(i)),
                base0 + 0.5*octbase + (octbase * float(j)),
                base0 + 0.5*octbase + (octbase * float(k))
                )

                result = myTree.findPosition(pos_find)

                if result is None:
                    print "No result for test at:", pos_find[0]/r_dust, pos_find[1]/r_dust, pos_find[2]/r_dust
                else:
                    print "Results for test at:", pos_find[0]/r_dust, pos_find[1]/r_dust, pos_find[2]/r_dust
                    if result is not None:
                        for ii in result:
                            print "    ", ii.name, ii.position[0]/r_dust,  ii.position[1]/r_dust, ii.position[2]/r_dust
                print '\n'

    f.close()

    print "\nOctree not fully supported, exiting\n"
    exit()
    return
