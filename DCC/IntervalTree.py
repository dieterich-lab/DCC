"""
Intersects ... faster.  Suports GenomicInterval datatype and multiple chromosomes.
Accept GenomicInterval object, with attributes chrom, start, end, strand (optional).
In case of unstranded data, interval.strand == '.'
"""

import math
import random


class IntervalTree(object):
    def __init__(self):
        self.chroms = {}

    def insert(self, interval, annotation=None):
        # This interval is the interval to construct the tree, e.g. gtf annotations
        chrom = interval.chrom
        start = interval.start
        end = interval.end
        strand = interval.strand
        if interval.chrom in self.chroms:
            self.chroms[chrom] = self.chroms[chrom].insert(start, end, strand,
                                                           annotation)  # self.chroms[chrom] is a IntervalNode object
        else:
            self.chroms[chrom] = IntervalNode(start, end, strand, annotation)

    def intersect(self, interval, report_func):
        # This interval from the query
        chrom = interval.chrom
        start = interval.start
        end = interval.end
        strand = interval.strand  # query interval strand can be '.'
        if chrom in self.chroms:
            self.chroms[chrom].intersect(start, end, strand,
                                         report_func)
            # use the intersect method of IntervalNode class, need make this function aware of strand

    def traverse(self, func):
        for item in self.chroms.values():
            item.traverse(func)


class IntervalNode(object):
    def __init__(self, start, end, strand=None, annotation=None):
        self.priority = math.ceil((-1.0 / math.log(.5)) * math.log(-1.0 / (random.uniform(0, 1) - 1)))
        self.start = start
        self.end = end
        self.maxend = self.end
        self.minend = self.end
        self.left = None
        self.right = None
        self.strand = strand
        self.annotation = annotation

    def insert(self, start, end, strand=None, annotation=None):
        root = self
        if start > self.start:
            # insert to right tree
            if self.right:
                self.right = self.right.insert(start, end, strand, annotation)
            else:
                self.right = IntervalNode(start, end, strand, annotation)
            # rebalance tree
            if self.priority < self.right.priority:
                root = self.rotateleft()
        else:
            # insert to left tree
            if self.left:
                self.left = self.left.insert(start, end, strand, annotation)
            else:
                self.left = IntervalNode(start, end, strand, annotation)
            # rebalance tree
            if self.priority < self.left.priority:
                root = self.rotateright()
        if root.right and root.left:
            root.maxend = max(root.end, root.right.maxend, root.left.maxend)
            root.minend = min(root.end, root.right.minend, root.left.minend)
        elif root.right:
            root.maxend = max(root.end, root.right.maxend)
            root.minend = min(root.end, root.right.minend)
        elif root.left:
            root.maxend = max(root.end, root.left.maxend)
            root.minend = min(root.end, root.left.minend)
        return root

    def rotateright(self):
        root = self.left
        self.left = self.left.right
        root.right = self
        if self.right and self.left:
            self.maxend = max(self.end, self.right.maxend, self.left.maxend)
            self.minend = min(self.end, self.right.minend, self.left.minend)
        elif self.right:
            self.maxend = max(self.end, self.right.maxend)
            self.minend = min(self.end, self.right.minend)
        elif self.left:
            self.maxend = max(self.end, self.left.maxend)
            self.minend = min(self.end, self.left.minend)
        return root

    def rotateleft(self):
        root = self.right
        self.right = self.right.left
        root.left = self
        if self.right and self.left:
            self.maxend = max(self.end, self.right.maxend, self.left.maxend)
            self.minend = min(self.end, self.right.minend, self.left.minend)
        elif self.right:
            self.maxend = max(self.end, self.right.maxend)
            self.minend = min(self.end, self.right.minend)
        elif self.left:
            self.maxend = max(self.end, self.left.maxend)
            self.minend = min(self.end, self.left.minend)
        return root

    ### TODO
    def intersect(self, start, end, strand, report_func):
        if strand == '.':  # unstranded data, not going to compare strand
            if start < self.end and end > self.start:
                report_func(self)
        else:
            if start < self.end and end > self.start and strand == self.strand:
                report_func(self)
        if self.left and start < self.left.maxend:
            self.left.intersect(start, end, strand, report_func)
        if self.right and end > self.start:
            self.right.intersect(start, end, strand, report_func)

    def traverse(self, func):
        if self.left: self.left.traverse(func)
        func(self)
        if self.right: self.right.traverse(func)
