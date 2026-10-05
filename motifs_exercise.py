# task 1

bases1 = "ACGT"

lecture_dna = [
    "TGACGTATAAGTTGCGATGGACGAGATAGCAGAGAATAGGCAACGAGAGATAAGCAG",
    "GACGGTAGCAGATAGACAGATGAAGAGTATGAATTGCACAGATAGCAGATAGCAGAT",
    "GGAGTGTGACGTAGCAGAGACGAAAGACGTAGAGTAGCAGTAGCAGATAGAGGGAGT",
    "TAGACAGTATAGAGACAGCGAGTCGGATAGCACCCAGTATGACGATAGCAATGACAG",
    "GCAGTAGAGCAGATTAGCATTGACAGATAGACGATTGGAGAGATGTGTGGATGACGA",
    "GGCAGGTAGCACACTGGGTCGATAAAGAGTAGCATAGAGACATAGACATATTTTAGC",
]

def count_matrix(motifs):
    l = len(motifs[0])
    counts = {"A": [0] * l, "C": [0] * l, "G": [0] * l, "T": [0] * l}

    for motif in motifs:
        for i in range(l):
            base1 = motif[i]
            counts[base1][i] += 1

    return counts

def score(motifs):
    counts = count_matrix(motifs)
    l = len(motifs[0])
    total = 0

    for i in range(l):
        biggest = max(counts["A"][i], counts["C"][i], counts["G"][i], counts["T"][i])
        total += biggest

    return total

def consensus(motifs):
    counts = count_matrix(motifs)
    l = len(motifs[0])
    result = ""

    for i in range(l):
        top = "A"
        for base1 in bases1:
            if counts[base1][i] > counts[top][i]:
                top = base1
        result += top

    return result

def hamming_distance(a, b):
    distance = 0

    for i in range(len(a)):
        if a[i] != b[i]:
            distance += 1

    return distance

def total_distance(pattern, sequences):
    l = len(pattern)
    total = 0

    for sequence in sequences:
        smallest = l
        for i in range(len(sequence) - l + 1):
            window = sequence[i:i + l]
            d = hamming_distance(pattern, window)
            if d < smallest:
                smallest = d
        total += smallest

    return total

red = ["TAAGTT", "TGAATT", "GGAGTG", "CGAGTC", "TGTGTG", "TGGGTC"]  # slide 19
best = ["AGATAG", "AGATAG", "AGATAG", "AGACAG", "AGATAG", "AGGTAG"]

print(score(red))                                # 26
print(consensus(best), score(best))              # AGATAG 34
print(hamming_distance("TAAGTT", "TGAATT"))      # 2
print(total_distance("TGCGTT", lecture_dna))     # 13

# task 2

class MotifProfile:
    def __init__(self, motifs, pseudocount=1):
        self.l = len(motifs[0])
        t = len(motifs)
        counts = count_matrix(motifs)
        self.ppm = {"A": [], "C": [], "G": [], "T": []}

        for base1 in bases1:
            for i in range(self.l):
                prob = (counts[base1][i] + pseudocount) / (t + 4 * pseudocount)
                self.ppm[base1].append(prob)

    def lmer_probability(self, lmer):
        prob = 1

        for i in range(self.l):
            base1 = lmer[i]
            prob *= self.ppm[base1][i]

        return prob

    def most_probable_lmer(self, sequence):
        top_lmer = sequence[0:self.l]
        top_prob = self.lmer_probability(top_lmer)

        for i in range(1, len(sequence) - self.l + 1):
            window = sequence[i:i + self.l]
            prob = self.lmer_probability(window)

            if prob > top_prob:
                top_prob = prob
                top_lmer = window

        return top_lmer

    def consensus(self):
        result = ""

        for i in range(self.l):
            top = "A"
            for base1 in bases1:
                if self.ppm[base1][i] > self.ppm[top][i]:
                    top = base1
            result += top

        return result

profile = MotifProfile(["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"])
print(profile.consensus())                            # ATGCGTA
print(round(profile.lmer_probability("ATGCGTA"), 4))  # 0.0122

two = MotifProfile(["GTAC", "TTAA"])
print(two.most_probable_lmer("ACTGGATGACCC"))         # TGAC
print(round(two.lmer_probability("TGAC"), 4))         # 0.0093

from Bio import motifs
from Bio.Seq import Seq

bio = motifs.create([Seq(site) for site in ["ATCCGTA", "GTGCATA", "AAGCGTA", "ATGCGTG"]])
bio.pseudocounts = 1
print(bio.consensus)     # ATGCGTA
print(bio.pwm["A"])      # the same numbers as your profile.ppm["A"]
print(profile.ppm["A"])

# task 3

import random
import itertools
from Bio import SeqIO

sequences = [str(record.seq) for record in SeqIO.parse("planted_motif.fasta", "fasta")]

class MotifFinder:
    def __init__(self, sequences, l, seed=None):
        self.sequences = sequences
        self.l = l
        self.rng = random.Random(seed)

        self.windows = []
        for sequence in sequences:
            sequence_windows = []
            for i in range(len(sequence) - l + 1):
                sequence_windows.append(sequence[i:i + l])
            self.windows.append(sequence_windows)

    def total_distance(self, pattern):
        total = 0

        for sequence_windows in self.windows:
            smallest = self.l
            for window in sequence_windows:
                d = hamming_distance(pattern, window)
                if d < smallest:
                    smallest = d
            total += smallest

        return total

    def median_string(self):
        top_pattern = ""
        top_dist = len(self.sequences) * self.l + 1

        for letters in itertools.product("ACGT", repeat=self.l):
            pattern = "".join(letters)
            d = self.total_distance(pattern)
            if d < top_dist:
                top_dist = d
                top_pattern = pattern

        return top_pattern, top_dist

    def randomized_search(self):
        current = []
        for sequence_windows in self.windows:
            current.append(self.rng.choice(sequence_windows))
        current_score = score(current)

        improving = True
        while improving:
            profile = MotifProfile(current, 1)

            new_motifs = []
            for sequence in self.sequences:
                new_motifs.append(profile.most_probable_lmer(sequence))
            new_score = score(new_motifs)

            if new_score > current_score:
                current = new_motifs
                current_score = new_score
            else:
                improving = False

        return current, current_score

    def best_of(self, runs):
        top_motifs = []
        top_score = -1

        for r in range(runs):
            got_motifs, got_score = self.randomized_search()
            if got_score > top_score:
                top_score = got_score
                top_motifs = got_motifs

        return top_motifs, top_score

# what to do 1

print("whattodo1:")
finder = MotifFinder(lecture_dna, 6, seed=1)
print(finder.median_string())
lec_motifs, lec_score = finder.best_of(100)
print(consensus(lec_motifs), lec_score)

# what to do 2

print("whattodo2:")
planted = MotifFinder(sequences, 7, seed=1)
median_pattern, median_dist = planted.median_string()
print(median_pattern, median_dist)

pl_motifs, pl_score = planted.best_of(100)
for m in pl_motifs:
    print(m)
print(consensus(pl_motifs), pl_score)
print(consensus(pl_motifs) == median_pattern)

# what to do 3

print("whattodo3:")
ok_runs = 0
for r in range(200):
    got_motifs, got_score = planted.randomized_search()
    if consensus(got_motifs) == median_pattern:
        ok_runs += 1
p = ok_runs / 200
print("p =", p)

print("R    measured    predicted")
for R in [5, 10, 20, 50]:
    hits_ok = 0
    for k in range(50):
        got_motifs, got_score = planted.best_of(R)
        if consensus(got_motifs) == median_pattern:
            hits_ok += 1
    measured = hits_ok / 50
    predicted = 1 - (1 - p) ** R
    print(R, "   ", measured, "   ", round(predicted, 2))

# what to do 4

print("whattodo4:")
found = motifs.create([Seq(m) for m in pl_motifs])  # the motifs from best_of(100)
found.pseudocounts = 1
pssm = found.pssm

for fpr in [0.01, 0.001]:
    threshold = pssm.distribution().threshold_fpr(fpr)
    total_hits = 0
    found_sites = 0
    reverse_hits = 0

    for k in range(len(sequences)):
        sequence = sequences[k]
        for position, hit_score in pssm.search(Seq(sequence), threshold=threshold):
            total_hits += 1
            if position < 0:
                reverse_hits += 1
            else:
                if sequence[position:position + 7] == pl_motifs[k]:
                    found_sites += 1

    expected = 2 * 10 * 54 * fpr
    print("fpr", fpr, "threshold", round(threshold, 2))
    print("hits:", total_hits, "found sites:", found_sites, "reverse strand:", reverse_hits, "expected by chance:", round(expected, 1))