#! /usr/bin/env python

from __future__ import print_function

# Returns intersect(A,B)/max(A,B), intersect(A,B)/min(A,B), intersect(A,B)/union(A,B), max. breakpoint offset, overlap length
def overlapMetrics((s1, e1), (s2, e2)):
    bpOffset = max(abs(s2-s1), abs(e2-e1))
    overlapLen = float(min(e1, e2) - max(s1, s2))
    lenA = float(e1-s1)
    lenB = float(e2-s2)
    if (e1 < s2) or (s1 > e2) or (lenA <= 0) or (lenB <= 0) or (overlapLen <= 0):
        return (0, 0, 0, bpOffset, 0)
    lenUnion = float(max(e1, e2)-min(s1, s2))
    return(overlapLen/max(lenA, lenB), overlapLen/min(lenA, lenB), overlapLen/lenUnion, bpOffset, overlapLen)

# Checks if an overlap is valid with respect to a given required reciprocal overlap and breakpoint offset
def overlapValid((s1, e1), (s2, e2), reciprocalOverlap, maxOffset):
    (maxO, _, _, bpOffset, _) = overlapMetrics((s1, e1), (s2, e2))
    return (maxO >= reciprocalOverlap) and (bpOffset <= maxOffset)


if __name__ == "__main__":
    import sys
    print("Enter input intervals: sv1Start, sv1End, sv2Start, sv2End")
    (svS1, svE1, svS2, svE2) = tuple(int(x.strip()) for x in sys.stdin.readline().split(','))
    print("intersect(A,B)/max(A,B), intersect(A,B)/min(A,B), intersect(A,B)/union(A,B), max. breakpoint offset, overlap length")
    print(overlapMetrics((svS1, svE1), (svS2, svE2)))
