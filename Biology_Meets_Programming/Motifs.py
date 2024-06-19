import math
import random


#generate a count matrix from an arbitrary list of strings Motifs
# Input:  A set of kmers Motifs
# Output: Count(Motifs)
def Count(Motifs):
    count = {}
    for nt in "ACGT":
        count[nt] = [0]*len(Motifs[0])
    for motif in Motifs:
        for index, let in enumerate(motif):
            count[let][index] += 1
    return count

# Input:  A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.
def Profile(Motifs):
    profile = {}
    count = Count(Motifs)
    for key in count.keys():
        profile[key] = [i/len(Motifs) for i in count[key]]
    return profile

def Consensus(Motifs):
    count = Count(Motifs)
    consensus_str = ""
    for posit in range(len(count["A"])):
        consensus_str += max({'A':count['A'][posit], 'C':count['C'][posit], 'G':count['G'][posit], 'T':count['T'][posit]}, 
                             key = {'A':count['A'][posit], 'C':count['C'][posit], 'G':count['G'][posit], 'T':count['T'][posit]}.get)
    return consensus_str

def Score(Motifs):
    count = Count(Motifs)
    consensus = Consensus(Motifs)
    score = 0
    for index, value in enumerate(consensus):
        for key in count.keys():
            if (key != value) & (count[key][index] != 0):
                score += count[key][index]
    return score

# Input:  String Text and profile matrix Profile
# Output: Pr(Text, Profile) which is probability of that specific text/motif given the profile matrix
def Pr(Text, Profile):
    p = 1
    for index, value in enumerate(Text):
        p = p* Profile[value][index]
    return p


def ProfileMostProbableKmer(text, k, profile):
    max_probability = -1  # Initialize max probability as a negative value
    most_probable_kmer = ""
    for start in range(len(text) - k + 1):
        kmer = text[start:start+k]
        probability = Pr(Text = kmer, Profile = profile)  # Use the Pr function to compute probability
        if probability > max_probability:
            max_probability = probability
            most_probable_kmer = kmer
    return most_probable_kmer


# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
#this gives an approximate solution, not the accurate/best that we would find with a brute force search (all possibs)
#because this is based on the profile from one sequence, the zeroes in that profile will cause issues finding 
#best motif really in the other seqs
def GreedyMotifSearch(Dna, k, t):
    BestMotifs = []
    #for every seq, set the first kmer as the first best motif
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])

    for i in range(len(Dna[0])-k+1): #loop over positions in seq
        Motifs = []
        Motifs.append(Dna[0][i:i+k]) #motifs is set to each kmer in the first sequence
        #for each kmer in the first sequence, 
        for j in range(1, t):
            #construct a profile based on the kmer in the first sequence
            P = Profile(Motifs[0:j])
            #use the first seq's kmer profile to find the most probable kmer in the rest of the seqs
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
            #Motifs ends up having 1 kmer per sequence. 1st seq determines profile, the kmers chosen from other seqs
            #are based on the profile of of the first seq
        #now we check if the score of the current collection of motifs is better than BestMotifs
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

##################
def EntropyMatrix(Motifs):
    profile = Profile(Motifs)
    entropy = {'A':[], 'C':[], 'G':[], 'T':[]}
    for key in profile.keys():
        entropylist = []
        for index in range(len(profile[key])):
            if profile[key][index] <= 0:
                entropylist.append(0)
            else:
                entropylist.append(profile[key][index]*math.log2(profile[key][index]))
        entropy[key] = entropylist
        #entropy[key] = [i*math.log2(i) for i in profile[key]]
    return entropy

def EntropyScore(Motifs):
    entropymatrix = EntropyMatrix(Motifs)
    score = 0
    for index, value in enumerate(Motifs[0]):
        score += -(sum((entropymatrix['A'][index], entropymatrix['C'][index], entropymatrix['G'][index], entropymatrix['T'][index])))
    return score
######################

# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
def CountWithPseudocounts(Motifs):
    count = {}
    for nt in "ACGT":
        count[nt] = [1]*len(Motifs[0])
    for motif in Motifs:
        for index, let in enumerate(motif):
            count[let][index] += 1
    return count

# Input:  A set of kmers Motifs
# Output: ProfileWithPseudocounts(Motifs)
def ProfileWithPseudocounts(Motifs):
    profile = {} # output variable
    count = CountWithPseudocounts(Motifs)
    for key in count.keys():
        #add 4 for 1 pseudocount per base
        profile[key] = [i/(len(Motifs)+4) for i in count[key]]
    return profile

def ConsensuswithPseudoCounts(Motifs):
    count = CountWithPseudocounts(Motifs)
    consensus_str = ""
    for posit in range(len(count["A"])):
        consensus_str += max({'A':count['A'][posit], 'C':count['C'][posit], 'G':count['G'][posit], 'T':count['T'][posit]}, 
                             key = {'A':count['A'][posit], 'C':count['C'][posit], 'G':count['G'][posit], 'T':count['T'][posit]}.get)
    return consensus_str

def ScorewithPseudoCounts(Motifs):
    count = CountWithPseudocounts(Motifs)
    consensus = ConsensuswithPseudoCounts(Motifs)
    score = 0
    for index, value in enumerate(consensus):
        for key in count.keys():
            if (key != value) & (count[key][index] != 0):
                score += count[key][index]
    return score

# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    BestMotifs = []
    #for every seq, set the first kmer as the first best motif
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])

    for i in range(len(Dna[0])-k+1): #loop over positions in seq
        Motifs = []
        Motifs.append(Dna[0][i:i+k]) #motifs is set to each kmer in the first sequence
        #for each kmer in the first sequence, 
        for j in range(1, t):
            #construct a profile based on the kmer in the first sequence
            P = ProfileWithPseudocounts(Motifs[0:j])
            #use the first seq's kmer profile to find the most probable kmer in the rest of the seqs
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
            #Motifs ends up having 1 kmer per sequence. 1st seq determines profile, the kmers chosen from other seqs
            #are based on the profile of of the first seq
        #now we check if the score of the current collection of motifs is better than BestMotifs
        if ScorewithPseudoCounts(Motifs) < ScorewithPseudoCounts(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs


# Input:  A profile matrix Profile and a list of strings Dna
# Output: returns a list of the Profile-most probable k-mers in each string from Dna
def Motifs(Profile, Dna):
    kmer_len = len(Profile['A'])
    ProbableMotifs = []
    for str in Dna:
        ProbableMotifs.append(ProfileMostProbableKmer(text = str, k = kmer_len, profile = Profile))
    return ProbableMotifs


# Input:  A list of strings Dna, and integers k and t
# Output: uses random.randint to choose a random k-mer from each of t different strings Dna, and returns a list of t strings
def RandomMotifs(Dna, k, t):
    random_motifs = []
    for seq in range(t):
        random_start = random.randint(1, len(Dna[0])-k)
        random_motifs.append(Dna[seq][random_start:random_start+k])
    return random_motifs


# Input:  Positive integers k and t, followed by a list of strings Dna
# Output: RandomizedMotifSearch(Dna, k, t)
def RandomizedMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    #code below stops running as soon as the score of the motifs that we generate stops improving
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs 

# Input: A dictionary Probabilities, where keys are k-mers and values are the probabilities of these
# k-mers (which do not necessarily sum up to 1)
# Output: A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities
def Normalize(Probabilities):
    pr_sum = sum(Probabilities.values())
    normalized_prs = {}
    for key, value in Probabilities.items():
        normalized_prs[key] = value/pr_sum
    return normalized_prs


# Input:  A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers
# Output: A randomly chosen k-mer with respect to the values in Probabilities
def WeightedDie(Probabilities):
    random_num = random.uniform(0, 1)
    for kmer, prob in Probabilities.items():
        random_num -= prob
        if random_num <= 0:
            return kmer
        
# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: return a randomly generated k-mer from Text whose probabilities are generated from profile
def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {} 
    #range over all possible k-mers in Text, computing the probability of each one and placing this probability into a dictionary
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)


#GibbsSampler is a more cautious iterative algorithm that discards a single k-mer from the current set of motifs
# at each iteration and decides to either keep it or replace it with a new one
# Input:  Integers k (kmer len), t (num seqs), and N (iterations), followed by a collection of strings Dna
# Output: GibbsSampler(Dna, k, t, N)
def GibbsSampler(Dna, k, t, N):
    Motifs = RandomMotifs(Dna, k, t)
    BestMotifs = Motifs
    for round in range(N):
        #pick random int
        i = random.randint(0, t-1) 
        #randomly exclude a motif
        Motifs.pop(i)
        #profile of remaining motifs
        profile = ProfileWithPseudocounts(Motifs)
        #use profile to get most probable kmer from original sequence of popped motif
        kmer = ProfileMostProbableKmer(text = Dna[i], k = k, profile = profile)
        #add the most probable kmer back to motifs
        Motifs.insert(i, kmer)
        #see if score with new kmer is better than previous BestMotifs
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs



