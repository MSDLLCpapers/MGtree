# Constructors and methods for Node and Tree classes,
#   created from a Newick string.
# Alfredo Gonzalez; John M. Gaspar
# Version 1.0  2024-10-03

import sys
import re
import unittest
import os

class Node:
  '''
  Node class for a hierarchical tree. 
  Contains a unique name, an (optional) genotype,
    data label(s), pointers to parent and children Nodes,
    and basic counters.
  '''
  def __init__(self, id=None, name=None, genotype=None, branch_length=0,
      children=None, parent=None, record=None):
    self.id = id               # number unique to this node
    self.name = name           # string unique to this node
    self.genotype = genotype   # genotype (optional)
    self.branch_length = branch_length
    self.children = children if children is not None else []
    self.parent = parent       # parent node
    self.additional_info = []  # necessary?
    self.lcaTotal = 0          # count of reads aligning here plus the subtree (LCA only)
    self.lcaCount = 0          # count of reads aligning here (LCA only)
    self.readCount = 0.0       # count of reads aligning here (leafs only)
                               #   - may be fractional
    self.record = record if record is not None else []
                               # list of reads assigned here (LCA only)

  def append_info(self, info):
    self.additional_info.append(info)


class Tree:
    '''
    Class for instantiating a hierarchical tree of Nodes from a Newick string.
    Contains pointer to root Node, and dicts of all Nodes and leaf Nodes.
    '''
    def __init__(self, newick, saveGT=False, verbose=False):
      self.newick = newick
      self.nodes = dict()  # dict of all nodes
      self.current_id = 0  # arbitrary id for unnamed nodes
      self.root = self.parse(newick, verbose)
      self.leafs = dict()  # dict of leaf nodes
      self.findLeafs(self.root, saveGT)

    def findLeafs(self, node, saveGT):
      '''
      Identify leaf nodes in a Tree, save to leafs dict.
        If saveGT is False, remove genotype.
      '''
      if not node.children:
        self.leafs[node.name] = node
        if not saveGT:
          node.genotype = None
      for child in node.children:
        self.findLeafs(child, saveGT)

    def _saveName(self, node, name):
      '''
      Save node name. Check for forbidden characters.
        If '_' is present, save first token as genotype.
        Helper method for parse().
      '''
      name = name.strip()
      if name in self.nodes:
        sys.stderr.write('Error! Repeated node name %s\n' % name)
        sys.exit(-1)
      if re.search(r'[;(),:]', name):
        sys.stderr.write('Error! Node name has a forbidden '
          + 'character [;(),:]: %s\n' % name)
        sys.exit(-1)
      spl = name.split('_')
      if len(spl) > 2:
        sys.stderr.write('Error! Too many \'_\' characters in %s\n' % name)
        sys.exit(-1)
      if len(spl) > 1:
        node.genotype = spl[0]
      node.name = name

    def _addNode(self, node, name):
      '''
      Save node id, name, and branch_length.
        Add node to tree dict.
        Helper method for parse().
      '''
      # find unique id (number)
      while str(self.current_id) in self.nodes:
        sys.stderr.write('conflict with %d\n' % self.current_id)
        self.current_id += 1
      node.id = self.current_id

      # save name, or use id if no name given
      spl = name.split(':')
      if spl[0]:
        self._saveName(node, spl[0])
      else:
        self._saveName(node, str(self.current_id))
      # save branch length
      if len(spl) > 1:
        node.branch_length = float(spl[1])

      # save to nodes dict
      self.nodes[node.name] = node
      self.current_id += 1

    def _connectNodes(self, parent, child):
      '''
      Connect parent and child nodes.
        Helper method for parse().
      '''
      parent.children.append(child)
      child.parent = parent

    def _checkID(self, node, count):
      '''
      Check that this node has an id attribute.
        Return count + 1.
        Helper method for parse().
      '''
      if node.id is None:
        sys.stderr.write('Error! Incomplete node found\n')
        sys.exit(-1)
      for child in node.children:
        count = self._checkID(child, count)
      return count + 1

    def parse(self, newick, verbose=False):
      '''
      Parse a Newick string. Construct Nodes of the tree.
      '''
      newick = newick.replace('\n', '')  # remove new line characters
      root = Node()   # dummy node
      current = root  # pointer to current node
      count = 1       # count of nodes created

      token = re.split('([(),;])', newick)
      i = 0
      while i < len(token):

        if token[i] == '':
          # empty token
          if (i+1 < len(token) and token[i+1] == '(') \
              or (i>0 and token[i-1] == ';'):
            # do nothing
            pass
          else:
            # create unnamed node
            node = Node()
            count += 1
            self._connectNodes(current, node)
            self._addNode(node, '')

        elif token[i] == ')':
          # complete internal node, with next token (blank is OK)
          i += 1  # skip ahead
          name = ''
          if i < len(token):
            name = token[i]
          self._addNode(current, name)

          # pop the stack
          if not current.parent:
            sys.stderr.write('Error! Poorly formed Newick string:' \
              + ' %s %s\n' % (token[i], token[i+1]))
            sys.exit(-1)
          current = current.parent

        elif token[i] == '(':
          # create unnamed node
          node = Node()
          count += 1
          self._connectNodes(current, node)
          current = node

        elif token[i] == ',':
          pass

        elif token[i] == ';':
          break

        else:
          # nonempty node name: create new node
          node = Node()
          count += 1
          self._connectNodes(current, node)
          self._addNode(node, token[i])

        i += 1

      # check for failure
      if count == 1:
        sys.stderr.write('Error! No nodes created\n')
        sys.exit(-1)

      # remove root node(s) with only 1 child
      current = root
      while len(current.children) == 1:
        if current.name in self.nodes:
          del self.nodes[current.name]
        current = current.children[0]
        current.parent = None
        count -= 1
      # ensure root node has attributes
      if not current.id:
        self._addNode(current, '')  # could label as 'root'
      current.branch_length = 0

      # verify count of nodes; check for incomplete (no id)
      countCheck = self._checkID(current, 0)
      if countCheck != count:
        sys.stderr.write('Error! Mismatched node counts: ' \
          + '%d, %d\n' % (countCheck, count))
        sys.exit(-1)
      if verbose:
        sys.stderr.write('Nodes created: %d\n' % count)

      return current

    def generateNewick(self, node):
      '''
      Recursively generate a Newick string representation of the tree.
      '''
      newick = ''
      if node.children:
        newick = '(%s)' % ','.join(self.generateNewick(child) \
          for child in node.children)
      newick += '%s:%s' % (node.name, node.branch_length)
      if node == self.root:
        newick += ';'
      return newick

    def print_tree(self, node, level=0):
        """
        Recursively prints the tree structure to stderr starting from
        the given node, including relevant node data and records.
        """
        if node:
            sys.stderr.write('  ' * level)
            sys.stderr.write('%s:%.2f' % (node.name, node.branch_length))
            sys.stderr.write('(ID: %s, ' % node.id)
            sys.stderr.write('LCA count: %d, ' % node.lcaCount)
            sys.stderr.write('Genotype: %s, ' % node.genotype)
            if node.additional_info:
              sys.stderr.write(' [Info: %s]' % ', '.join(node.additional_info))
            if node.record:
              sys.stderr.write(' [Record: %s]' % ', '.join(map(str, node.record)))
            if node.name in self.leafs:
              sys.stderr.write(' LEAF')
            sys.stderr.write(')\n')

            for child in node.children:
                self.print_tree(child, level + 1)


    def get_node_info(self, name):
        node = self.nodes.get(name, None)
        if node:
            return{
                'id': node.id,
                'parent': node.parent.name if node.parent else None,
                'children' : [child.name for child in node.children],
                'additional_info': node.additional_info,
                'record': node.record
            }
        return None

    def append_info_to_node(self, name, info):
        node = self.nodes.get(name)
        if node:
            node.append_info(info)

    def append_info_to_node_and_parents(self, name, info):
        node = self.nodes.get(name)
        while node:
            node.append_info(info)
            node = node.parent

    def rename_node(self, name, new_name):
        node = self.nodes.get(name)
        if node:
            self.nodes[new_name] = self.nodes.pop(name)
            node.name = new_name

    def set_record_for_node(self, name, record):
        node = self.nodes.get(name)
        if node:
            node.record = record

    def get_child_node_records(self, node_name, filename):
        """
        Retrieves and prints the records of child nodes for a given node,
        and writes the results to a CSV file.
        """
        # Get node or output a warning if no node found
        node = self.nodes.get(node_name)
        if not node:
            sys.stderr.write('No node found with name: %s\n' % node_name)
            return

        # Collect child nodes that have records
        record_nodes = [(child.name, len(child.record)) for child in node.children if child.record]

        # Sort the list by record count in descending order
        record_nodes.sort(key=lambda x: x[1], reverse=True)

        # Compute total number of records
        total_records = sum(count for _, count in record_nodes)

        # Print the table
        sys.stderr.write('Total records for children of %s: %d\n' % (node_name, total_records))
        for child_name, count in record_nodes:
            sys.stderr.write('%s: %d\n' % (child_name, count))

        # Write output to csv file
        f = open(filename, 'w')
        f.write('Total Record Count,%d\n' % total_records)
        f.write('Node Name,Record Count\n')
        for n in record_nodes:
          f.write('%s,%d\n' % (n[0], n[1]))
        f.close()


    def updateGenotype(self, name, genotype):
      '''
      Update the genotype of the node with the given name.
      '''
      if name not in self.nodes:
        sys.stderr.write('Warning! Cannot update genotype for '
          + 'node "%s": not found\n' % name)
        return False
      if re.search(r'[;(),:_]', genotype):
        sys.stderr.write('Warning! Genotype has a forbidden '
          + 'character [;(),:_]: %s\n' % genotype)
        return False
      node = self.nodes[name]
      if not node.children:
        sys.stderr.write('Warning! Renaming genotype for leaf node'
          + ' "%s"\n' % name)
      node.genotype = genotype
      # update name to include new genotype
      acc = node.name
      spl = node.name.split('_')
      if len(spl) > 1:
        acc = spl[1]
      node.name = '%s_%s' % (genotype, acc)
      return True

    def increaseReadCount(self, name, val):
      '''
      Increase the readCount of the node with the given name.
      '''
      if name not in self.nodes:
        sys.stderr.write('Warning! Cannot increase readCount for '
          + 'node "%s": not found\n' % name)
        return
      self.nodes[name].readCount += val

    def getLeafReadCounts(self):
      '''
      Return sorted list of leaf node counts.
      '''
      info_list = [(name, node.readCount) for name, node in self.leafs.items()]
      # sort the list by readCount in descending order
      info_list.sort(key=lambda x: (x[1], x[0]), reverse=True)
      return info_list

    def _collectLCA(self, node):
      '''
      Save LCA total as sum of LCA counts of this node
        and all the nodes in its subtree.
        Helper method for getGenotypeInfo().
      '''
      node.lcaTotal = node.lcaCount
      for child in node.children:
        node.lcaTotal += self._collectLCA(child)
      return node.lcaTotal

    def _getGenotypeChildren(self, node, gts):
      '''
      Recursively save nodes that are genotyped in the
        given node's subtree to the gts list.
        Helper method for _collectGenotype().
      '''
      if node.genotype:
        gts.append(node.genotype)
        return gts
      for child in node.children:
        gts = self._getGenotypeChildren(child, gts)
      return gts

    def _collectGenotype(self, info_list, ambigLCA, dupLCA,
        node, parentGen):
      '''
      Recursively append lcaTotals for nodes with genotypes to info_list.
        Warn if an ancestor node has a genotype.
        Counters for ambiguous LCA counts unassigned to genotypes
        (i.e. in nodes without genotypes and with no genotyped ancestor),
        and for LCA counts reported multiple times.
        Helper method for getGenotypeInfo().
      '''
      # return if there are no LCA totals
      if not node.lcaTotal:
        return info_list, ambigLCA, dupLCA

      # save if genotype is given
      if node.genotype:
        if parentGen:
          sys.stderr.write('Warning! Node %s ' % node.name \
            + '(genotype %s) has genotyped ancestor(s): ' % node.genotype \
            + '%s\n' % parentGen)
          parentGen += ', '
          dupLCA += node.lcaTotal
        info_list.append((node.genotype, node.lcaTotal))
        parentGen += 'Node %s (genotype %s)' % (node.name, node.genotype)

      # ambiguous LCA counts: node with no genotype, not printed yet
      elif node.lcaCount and not parentGen:
        gts = []  # list of genotypes in descendant nodes
        gts = self._getGenotypeChildren(node, gts)
        if gts:
          genotype = 'ambig[%s]' % ','.join(gts)
        else:
          genotype = 'ambig[?]'
        info_list.append((genotype, node.lcaCount))
        ambigLCA += node.lcaCount

      for child in node.children:
        info_list, ambigLCA, dupLCA = self._collectGenotype(info_list,
          ambigLCA, dupLCA, child, parentGen)
      return info_list, ambigLCA, dupLCA

    def getGenotypeInfo(self):
      '''
      Return sorted list of counts for genotype nodes.
      '''
      # collect LCA totals for each node
      self._collectLCA(self.root)

      # recursively collect genotype nodes with LCA counts
      info_list = []
      ambigLCA = dupLCA = 0
      info_list, ambigLCA, dupLCA = self._collectGenotype(info_list,
        ambigLCA, dupLCA, self.root, '')

      # sort the list by lcaTotal in descending order
      info_list.sort(key=lambda x: x[1], reverse=True)
      return info_list, ambigLCA, dupLCA, self.root.lcaTotal

    def _returnLCA(self, node, d, count):
      '''
      Traverse tree recursively, return reference to LCA
        as the lowest node whose dict value matches count.
        Helper method for findLCA().
      '''
      for child in node.children:
        if child.name in d and d[child.name] == count:
          return self._returnLCA(child, d, count) 
      return node

    def _countAncestors(self, node, d):
      '''
      Add a count for all ancestor nodes to dict.
        Helper method for findLCA().
      '''
      while node:
        d[node.name] = d.get(node.name, 0) + 1
        node = node.parent

    def findLCA(self, names):
      '''
      Return a reference to the lowest common ancestor node
        of a given list of node names.
      '''
      d = dict()
      count = 0
      for name in names:
        if name not in self.nodes:
          sys.stderr.write('Warning! Cannot find node: %s\n' % name)
          continue
        # add counts to ancestor nodes
        self._countAncestors(self.nodes[name], d)
        count += 1
      # check that root matches count
      if self.root.name not in d or d[self.root.name] != count:
        sys.stderr.write('Error! Cannot find LCA for nodes:\n'
          + '  %s\n' % ','.join(names))
        sys.exit(-1)
      # return ref to LCA (via recursive tree traversal)
      return self._returnLCA(self.root, d, count)

    def increaseLCACount(self, names, val, qname):
      '''
      Find the lowest common ancestor node of given node names.
        Increase its lcaCount; append given qname.
      '''
      lcaNode = self.findLCA(names)
      lcaNode.record.append(qname)
      lcaNode.lcaCount += val
      return lcaNode

    def getRecords(self, node, d):
      '''
      Recursively add leaf node records to given dict.
        Method to assist qname-json dumping.
      '''
      if node.record:
        d[node.name] = node.record
      for child in node.children:
        self.getRecords(child, d)

class TestTrees(unittest.TestCase):
    def setUp(self):
        self.newick_string = "((Sphynx:0.7,(Calico:0.3, Tabby:0.3)OrangeCats:0.5,(Tuxedocat:0.2, Oreocat:0.4)BlackCats:0.3)Cats:1.0);"
        self.tree = Tree(self.newick_string, False, False)

    def test_append_info_to_node(self):
        sys.stderr.write('\nInitial Tree:\n')
        self.tree.print_tree(self.tree.root)
        self.tree.append_info_to_node('Tuxedocat', 'Record-A')
        node_info = self.tree.get_node_info('Tuxedocat')
        self.assertIn('Record-A', node_info['additional_info'])
        sys.stderr.write('\nFinal tree after appending to Tuxedocat:\n')
        self.tree.print_tree(self.tree.root)

    def test_append_info_to_node_and_parents(self):
        self.tree.append_info_to_node_and_parents('Oreocat', 'Record-B')
        node_info = self.tree.get_node_info('Oreocat')
        self.assertIn('Record-B', node_info['additional_info'])

        parent_info = self.tree.get_node_info('BlackCats')
        self.assertIn('Record-B', parent_info['additional_info'])

        parent_level2_info = self.tree.get_node_info('Cats')
        self.assertIn('Record-B', parent_level2_info['additional_info'])

        sys.stderr.write('\nFinal tree after appending to Oreocat:\n')
        self.tree.print_tree(self.tree.root)

    def test_append_record_to_first_common_ancestor(self):
        node_names = ['Sphynx', 'Calico', 'Tabby']
        record = 'Record-C'
        common_ancestor = self.tree.increaseLCACount(node_names, 1, record)
        self.tree.increaseLCACount(['Tabby', 'OrangeCats'], 1, 'Record-D')
        self.tree.increaseLCACount(['Tabby', 'OrangeCats'], 1, 'Record-E')
        self.assertIsNotNone(common_ancestor)
        self.assertIn(record, common_ancestor.record)
        self.tree.increaseReadCount('Tabby', 1)
        self.tree.updateGenotype('Tabby', 'Kitty')
        sys.stderr.write('\nFinal tree after appending record to the first common ancestor:\n')

        # Texts renaming function
        self.tree.rename_node('Tuxedocat', 'Tuxcat')
        self.tree.print_tree(self.tree.root)

        # Get and print the table
        sys.stderr.write('\nNode Info Table:\n')
        table = self.tree.getLeafReadCounts()
        for name, recordtally in table:
            sys.stderr.write('%s: %d\n' % (name, recordtally))

    def test_get_child_node_records(self):

        # Append some records to the nodes to set up for the test
        self.tree.set_record_for_node('Sphynx', ['Record-1', 'Record-2'])
        self.tree.set_record_for_node('Tabby', ['Record-3'])

        sys.stderr.write('\nTree before testing get_child_node_records:\n')
        self.tree.print_tree(self.tree.root)

        # Set a filename for the output CSV
        filename = 'test_get_child_node_records.csv'

        # Call the method to be tested
        self.tree.get_child_node_records('Cats', filename)

        # Check that the file was created
        self.assertTrue(os.path.exists(filename))

        # Clean up by removing the file
        os.remove(filename)

if __name__ == '__main__':
    unittest.main(argv=[''], exit=False)
