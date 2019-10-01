Etrf is a simple tool to find exact tandem repeats (i.e. without mismatches or
gaps in the repeat unit) in DNA sequences. It only has two parameters: the
maximum repeat unit length and the minimum total repeat length. For each unit
length, etrf scans an input sequence and obtains a list of non-overlapping
regions no less than twice of the unit length. For two overlapping regions
identified with different unit lengths, etrf chooses the longer one, or the one
found with the shorter unit length if the two regions are of equal length.

Unable to find impure tandem repeats, etrf doesn't replace more sophisticated
tools such as [TRF][trf] or [ULTRA][ultra]. Nonetheless, because etrf
implements an exact algorithm, it avoids ambiguity in the definition of repeats
and its behavior is predicable. Etrf is also faster. It can process a human
genome in 15 minutes on a single CPU thread.

[trf]: https://tandem.bu.edu/trf/trf.html
[ultra]: https://github.com/TravisWheelerLab/ULTRA
