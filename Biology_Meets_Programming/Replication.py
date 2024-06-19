### Patterns and k-mers ###
#Returns the number of times a substring occurs in a string and accounts for overlapping occurrences.
def PatternCount(Text, Pattern):
    count = 0
    #range stops at last possible position for length of pattern in text
    for i in range(len(Text)-len(Pattern)+1):
        if Text[i:i+len(Pattern)] == Pattern:
            count += 1
    return count 

#Returns the string(s) with the highest number of occurrences for a k-mer of a specific length
def FrequentWords(Text, k):
    words = []
    freq = FrequencyMap(Text, k)
    m = max(freq.values())
    for key in freq:
        # add each key to words whose corresponding frequency value is equal to m
        if freq[key] == m:
            words.append(key)
    return words
#Generates a frequency map for a k-mer of a specific length
def FrequencyMap(Text, k):
    freq = {}
    n = len(Text)
    for i in range(n-k+1):
        Pattern = Text[i:i+k]
        freq[Pattern] = 0
    for pattern in freq:
        for i in range(len(Text)-k+1):
            if Text[i:i+k] == pattern:
                freq[pattern] +=1
    return freq



def ReverseComplement(Pattern):   
    Pattern = Reverse(Pattern)
    Pattern = Complement(Pattern)
    return Pattern
def Reverse(Pattern):
    ReversedPattern = ""
    for char in Pattern:
        ReversedPattern = char + ReversedPattern
    return ReversedPattern
def Complement(Pattern):
    Pattern = Pattern.replace("T", "a").replace("A", "t").replace("C", "g").replace("G","c").upper()
    return Pattern

#Finds all the starting positions of a pattern in a string
def PatternMatching(Pattern, Genome):
        positions = []
        for i in range(len(Genome) - len(Pattern) + 1):
            if Genome[i:i+len(Pattern)] == Pattern:
                positions.append(i)
        return positions

#THIS WILL NOT WORK FOR AN ACTUAL GENOME (TIME OUT/INEFFICIENT): When the window is half the length of the genome, this will count occurrences of the specified symbol (nt)
#in each window with starting index i
def SymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]
    for i in range(n):
        array[i] = PatternCount(Pattern = symbol, Text = ExtendedGenome[i:i+(n//2)])
    return array

#genome = open("Ecoli_genome.txt", 'r').read()
#running the function below times out
#print(SymbolArray(Genome = genome, symbol = "C"))

def FasterSymbolArray(Genome, symbol):
    array = {}
    n = len(Genome)
    ExtendedGenome = Genome + Genome[0:n//2]

    # look at the first half of Genome to compute first array value
    array[0] = PatternCount(Pattern = symbol, Text = Genome[0:n//2])

    for i in range(1, n):
        # start by setting the current array value equal to the previous array value
        array[i] = array[i-1]

        # the current array value can differ from the previous array value by at most 1. Basically we account for
        # the removal of one chr, and the addition of a new chr for the count in each window. This avoids the need 
        # to search for the symbol in the entire window for every window

        #if the prev chr that is no longer in the window was equal to the symbol, then we subtract 1 from the count
        #this removes from the count the matching symbol that's no longer in the new window
        if ExtendedGenome[i-1] == symbol:
            array[i] = array[i]-1

        #look at whether the new chr at the end of the new window matches the symbol. If it does, increase the count
        if ExtendedGenome[i+(n//2)-1] == symbol:
            array[i] = array[i]+1
    return array

#SkewArray will track C vs G
#I thought SkewArray would be array output, but then made quick switch to list output
def SkewArray(Genome):
    array = {}
    array[0] = 0
    for i in range(0, len(Genome)):
        #i is genome position
        #i+1 is array position
        if (Genome[i] == "A") | (Genome[i] == "T"):
            array[i+1] = array[i]
        if Genome[i] == "C":
            array[i+1] = array[i] - 1
        if Genome[i] == "G":
            array[i+1] = array[i] + 1
        else:
            next
    out = list(array.values())
    return out


#locate ori based on where skew diagram attains a minimum
def MinimumSkew(Genome):
    skew_array = SkewArray(Genome = Genome)
    # range over the length of the skew array and keep all indices equal the min value of the skew array
    return [index for index, val in enumerate(skew_array) if val == min(skew_array)]

#locate replication terminus based on where skew diagram attains a maximum
def MaximumSkew(Genome):
    skew_array = SkewArray(Genome = Genome)
    # range over the length of the skew array and keep all indices equal the min value of the skew array
    return [index for index, val in enumerate(skew_array) if val == max(skew_array)]


def HammingDistance(p,q):
    count = 0
    for i, n in zip(p,q):
        if i != n:
            count += 1
    return count

# Input:  Strings Pattern and Text along with an integer d
# Output: A list containing all starting positions where Pattern appears
# as a substring of Text with at most d mismatches
def ApproximatePatternMatching(Text, Pattern, d):
    positions = []
    #slide a window of pattern length along text until last position of window start
    for i in range(len(Text) - len(Pattern) + 1):
        if HammingDistance(Pattern, Text[i:i + len(Pattern)]) <= d:
            positions.append(i)
    return positions

# Input:  Strings Pattern and Text, and an integer d
# Output: The number of times Pattern appears in Text with at most d mismatches
def ApproximatePatternCount(Pattern, Text, d):
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        if HammingDistance(Pattern, Text[i:i + len(Pattern)]) <= d:
            count += 1
    return count

