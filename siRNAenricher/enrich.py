import textwrap
import os.path
import csv
import numpy as np


class RefSeq(dict):
    """
    Reference sequence/s object - subclass of dict
    Format = {header:DNA(seq)}
    """

    def load_ref_file(self, in_file):
        """
        Load a FASTA-formatted reference file.  RNA is converted to DNA.
        :param in_file: path to FASTA formated reference fi;e
        :param strip_dash: strip dashes (alignment FASTA introduced) from references
        """
        with open(in_file, "r") as f:
            header = ""
            for line in f:
                line = line.strip()
                if line == "":
                    pass
                elif line[0] == ">":
                    header = line[1:]
                    self[header] = []
                else:
                    self[header].append(line)

        f.close()
        for k, v in self.items():
            self[k] = "".join(v)

    def write_ref_file(self, out_file):
        with open(out_file, "w") as f:
            for header, seq in self.items():
                f.write(">{}\n".format(header))
                for line in textwrap.wrap(seq, 120):
                    f.write("{}\n".format(line))
        f.close()


class DNA(object):
    """
    DNA class
    """

    dna_alphabet = set("AGCTN")

    def __init__(self, sequence):
        self.sequence = sequence.upper()

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, key):
        return self.sequence[key]

    def __hash__(self):
        return hash(self.sequence)

    def __repr__(self):
        return self.sequence

    def __eq__(self, other):
        return self.sequence == other.sequence


class SingleAlignment(object):
    """
    Single sRNA read alignment class
    """

    def __init__(self, srna, position, strand, times_aligned, indv_alignments):
        self.srna = srna
        self.position = position
        self.strand = strand
        self.times_aligned = times_aligned
        self.indv_alignments = indv_alignments

    def srna_len(self):
        return len(self.srna)

    def standard_error(self):  # TODO: check - likely don't neeed here
        return np.std(self.indv_alignments, ddof=1) / np.sqrt(
            np.size(self.indv_alignments)
        )

    def mean_alignments(self):
        return np.mean(self.indv_alignments)

    def __str__(self):
        return "{0}\t{1}\t{2}\t{3}\t{4}".format(
            self.srna,
            self.position,
            self.strand,
            self.times_aligned,
            self.indv_alignments,
        )

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        return (
            self.srna == other.srna
            and self.position == other.position
            and self.strand == other.strand
            and self.times_aligned == other.times_aligned
            and np.array_equal(self.indv_alignments, other.indv_alignments)
        )


class SingleRefProfile(object):
    """
    Single reference sequence class
    """

    def __init__(self):
        self.ref_len = 0
        self.all_alignments = []

    def __str__(self):
        return "{0}\t{1}".format(self.ref_len, self.all_alignments)

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        return (
            self.ref_len == other.ref_len
            and self.all_alignments == other.all_alignments
        )


class RefProfiles(object):
    """
    All references in a file class
    """

    def __init__(self):
        self.srna_len = 0
        self.replicates = 0
        self.single_ref_profiles = {}

    def __str__(self):
        return "{0}\t{1}\t{2}".format(
            self.srna_len, self.replicates, self.single_ref_profiles
        )

    def __eq__(self, other):
        return (
            self.replicates == other.replicates
            and self.single_ref_profiles == other.single_ref_profiles
            and self.srna_len == other.srna_len
        )

    def load_single_ref_profiles(self, in_file):
        """Loads a scram2 profile file for a single sRNA length
        Args:
            in_file (string): scram2 profile file path
        """
        with open(in_file, "r") as in_handle:
            reader = csv.reader(in_handle, delimiter=",")
            for row in reader:
                if row[0] == "Header":
                    continue
                header = row[0]
                ref_len = int(row[1])
                srna = DNA(row[2])
                position = int(row[3])
                strand = row[4]
                times_aligned = int(row[5])
                indv_alignments = np.array([float(x) for x in row[6:]])
                sa = SingleAlignment(
                    srna, position, strand, times_aligned, indv_alignments
                )
                if header not in self.single_ref_profiles:
                    self.single_ref_profiles[header] = SingleRefProfile()
                    self.single_ref_profiles[header].ref_len = ref_len
                    self.srna_len = len(srna)
                self.single_ref_profiles[header].all_alignments.append(sa)
            self.replicates = len(sa.indv_alignments)


class Enrichment(object):
    """
    Enrichment class
    """

    def __init__(
        self, window=200, cutoff=30, abund_count=5, strand_ratio=0.2, padding=30
    ):
        self.window = window
        self.cutoff = cutoff
        self.abund_count = abund_count
        self.padding = padding
        self.strand_ratio = strand_ratio
        self.expanded_results = []
        self.collapsed_results = {}

    def find_enriched_regions(self, scram2_alignment):
        win_count = 0
        start_pos = 0
        positive_strand = 0
        negative_strand = 0
        a = RefProfiles()
        a.load_single_ref_profiles(scram2_alignment)
        for header, alignments in a.single_ref_profiles.items():
            for sa in alignments.all_alignments:
                if sum(sa.indv_alignments) / len(sa.indv_alignments) >= self.cutoff:
                    if win_count == 0:
                        win_count = 1
                        start_pos = sa.position
                        if sa.strand == "+":
                            positive_strand += 1
                        else:
                            negative_strand += 1
                    elif sa.position - start_pos <= self.window:
                        win_count += 1
                        if sa.strand == "+":
                            positive_strand += 1
                        else:
                            negative_strand += 1
                if sa.position - start_pos > self.window:
                    if positive_strand == 0 or negative_strand == 0:
                        pass
                    elif (
                        min(positive_strand, negative_strand)
                        / max(positive_strand, negative_strand)
                        < self.strand_ratio
                    ):
                        pass
                    elif win_count >= self.abund_count:
                        self.expanded_results.append(
                            [header, start_pos, start_pos + self.window]
                        )
                    win_count = 0
                    positive_strand = 0
                    negative_strand = 0

    def collapse_enriched_regions(self):
        end = 0
        header = ""
        for i in self.expanded_results:
            if i[0] not in self.collapsed_results:
                self.collapsed_results[i[0]] = [
                    i[1:],
                ]
            elif i[1] - self.collapsed_results[i[0]][-1][-1] < self.window:
                self.collapsed_results[i[0]][-1][-1] = i[2]
            else:
                self.collapsed_results[i[0]].append(i[1:])

    def extract_enriched(self, ref, out_fa):
        if len(self.collapsed_results) == 0:
            print("\nNo siRNA enriched regions identified\n")
        else:
            b = RefSeq()
            b.load_ref_file(ref)
            c = RefSeq()
            for k, v in self.collapsed_results.items():
                count = 0
                for i in v:
                    count += 1
                    header = (
                        k
                        + "_siRNA_enriched_"
                        + str(i[0] - self.padding)
                        + "-"
                        + str(i[1] + self.padding)
                    )

                    seq = b[k][i[0] - self.padding : i[1] + self.padding]
                    c[header] = seq
            c.write_ref_file(out_fa)


def extract_enriched_seqs(
    scram_alignment_file,
    reference_fa,
    output_fa,
    window=200,
    cutoff=30,
    abund_count=5,
    strand_ratio=0.2,
    padding=30,
):
    a = Enrichment(window, cutoff, abund_count, strand_ratio, padding)
    print(
        "Finding enriched regions with window = {0}, cutoff = {1}, abundance count = {2}, minimum strand ratio = {3} and padding = {4}".format(
            window, cutoff, abund_count, strand_ratio, padding
        )
    )
    a.find_enriched_regions(scram_alignment_file)
    a.collapse_enriched_regions()
    print("Loading reference file and generating FASTA file with enriched regions")
    a.extract_enriched(reference_fa, output_fa)
    print("Complete!")
