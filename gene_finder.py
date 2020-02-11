    # -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@mmadanguit: Marion Madanguit

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):
    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide

    *** Pre-written unit tests felt sufficient because they ran through all four
        nucleotides.
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('G')
    'C'
    >>> get_complement('T')
    'A'
    """

    if nucleotide == 'A':
        return 'T'
    if nucleotide == 'C':
        return 'G'
    if nucleotide == 'G':
        return 'C'
    if nucleotide == 'T':
        return 'A'

def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string

    *** Pre-written unit tests felt sufficient because both test dna strands
        included all four nucleotides.
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """

    strand = ''
    for char in dna[ : :-1]:
        complement = get_complement(char)
        strand += complement
    return strand

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string

    *** I added two  unit tests: one that has the "TAA" stop codon and another
        that has no stop codons (two cases that the pre-written unit tests did
        not account for).
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("ATGAGATAAG")
    'ATGAGA'
    >>> rest_of_ORF("ATGAGATGGG")
    'ATGAGATGGG'
    """

    orf = ''
    for i in range(0, len(dna), 3):
        if dna[i:i+3] == 'TAA' or dna[i:i+3] == 'TAG' or dna[i:i+3] == 'TGA':
            return orf
        orf += dna[i:i+3]
    return orf

def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    *** Pre-written unit test felt sufficient because test dna strands included
        both non-nested and nested ORFs.
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """

    list_orfs = []
    i = 0
    while i < len(dna):
        if dna[i:i+3] == 'ATG':
            orf = rest_of_ORF(dna[i:])
            list_orfs.append(orf)
            i += len(orf)
        else:
            i += 3
    return list_orfs

def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    *** Pre-written unit test felt sufficient because test dna strand included
        ORFs in multiple reading frames.
    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """

    list_orfs = []
    for i in range(0, 3):
        orf = find_all_ORFs_oneframe(dna[i:])
        list_orfs += orf
    return list_orfs

def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    *** Pre-written unit test felt sufficient because test dna strand included
        ORFs on both strands.
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """

    list_forward = find_all_ORFs(dna)
    list_reverse = find_all_ORFs(get_reverse_complement(dna))
    return list_forward + list_reverse

def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string

    *** Pre-written unit test felt sufficient because test dna strand had multiple
        ORFs of different lengths.
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """

    list_orfs = find_all_ORFs_both_strands(dna)
    if list_orfs:
        return max(list_orfs, key = len)
    else:
        return ''

def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF

    *** To test this, I ran the function about 20 times to ensure that each time
        it outputted a practical length and never ran into an error.
    """

    list_orfs = []
    for i in range(num_trials):
        new_string = shuffle_string(dna)
        orf = longest_ORF(new_string)
        list_orfs.append(orf)
    longest = max(list_orfs, key = len)
    return len(longest)

def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

    *** Pre-written unit tests felt sufficient because test dna strands varied
        in length and checked multiple amino acids.
    >>> coding_strand_to_AA("ATGCGA")
    'MR'
    >>> coding_strand_to_AA("ATGCCCGCTTT")
    'MPA'
    """

    protein = ''
    i=0
    while i < len(dna):
        if len(dna[i:i+3]) == 3:
            amino_acid = aa_table[dna[i:i+3]]
            protein += amino_acid
        i += 3
    return protein

def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.

    *** To test this, I used a sequence of DNA from Salmonella Enterica and
        checked the returned amino acid sequences with online resources.
    """

    list_proteins = []
    threshold = longest_ORF_noncoding(dna, 1500)
    all_orfs = find_all_ORFs_both_strands(dna)
    for orf in all_orfs:
        if len(orf) > threshold:
            protein = coding_strand_to_AA(orf)
            list_proteins.append(protein)
    return list_proteins

if __name__ == "__main__":
    import doctest
    # doctest.testmod()
    # doctest.run_docstring_examples(coding_strand_to_AA, globals(), verbose=True)

    # x = longest_ORF_noncoding("ATGCGAATGTAGCATCAAA", 3)
    # print(x)

    from load import load_seq
    dna = load_seq("./data/X73525.fa")
    print(gene_finder(dna))
