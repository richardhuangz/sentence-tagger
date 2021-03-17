#!/usr/bin/python3

import sys
import os
import math

# Compute P(w|t)
def part1a(file):

    global t_freqs
    t_freqs = {}
    global wt_freqs
    wt_freqs = {}
    global wt_probs
    wt_probs = {}

    with open(file, "r") as fp:
        for line in fp:

            parts = line.strip().split("//")    # [word], [tag]

            # tag already in dict
            if parts[1] in t_freqs:             
                t_freqs[parts[1]] += 1
            else:
                t_freqs[parts[1]] = 1

            # pairing of word and tag already in dict
            if line.strip() in wt_freqs:                
                wt_freqs[line.strip()] += 1
            else:
                wt_freqs[line.strip()] = 1        
    
    # calculate probs
    for wt in wt_freqs:
        parts = wt.split("//")
        wt_probs[wt] = wt_freqs[wt] / t_freqs[parts[1]]


# Compute P(t_i|t_{i−1})
def part1b(file):

    global ti_1_ti_freqs
    ti_1_ti_freqs = {}
    global ti_1_freqs
    ti_1_freqs = {}
    global ti_1_ti_probs
    ti_1_ti_probs = {}

    with open(file, "r") as fp:

        prev = "...EOS..."   # t_{i-1} tag

        for line in fp:

            parts = line.strip().split("//")

            # previous tag already in dict
            if prev in ti_1_freqs:
                ti_1_freqs[prev] += 1
            else:
                ti_1_freqs[prev] = 1
            
            ti_1_ti = parts[1] + "//" + prev    # tag//prev_tag

            # this pairing of prev and curr tag in dict
            if ti_1_ti in ti_1_ti_freqs:
                ti_1_ti_freqs[ti_1_ti] += 1
            else:
                ti_1_ti_freqs[ti_1_ti] = 1     

            # update prev
            prev = parts[1] 
    
    # calculate probs
    for ti_1_ti in ti_1_ti_freqs:
        prev = ti_1_ti.split("//")[1]
        ti_1_ti_probs[ti_1_ti] = ti_1_ti_freqs[ti_1_ti] / ti_1_freqs[prev]
    

# This is fake. Tag an input sentence using weights from 1a*1b
def part1c(sentence):
    taggedsentence = ""
    prev = "...EOS..."
    for word in sentence.split():
        taggedsentence += word + "//"
        maxprob = 0
        maxtag = ""
        for key in wt_probs:
            parts = key.split("//")
            
            if parts[0] != word:
                continue

            for key2 in ti_1_ti_probs:
                parts2 = key2.split("//")

                if parts2[0] != parts[1]:
                    continue

                curr_ti_ti1 = parts2[0] + "//" + prev

                if curr_ti_ti1 not in ti_1_ti_probs:
                    continue
                
                if wt_probs[key] * ti_1_ti_probs[key2] > maxprob:
                    maxprob = wt_probs[key] * ti_1_ti_probs[key2]
                    maxtag = parts[1]

        taggedsentence += maxtag + " "
        prev = maxtag
    
    return taggedsentence

# Use 1c on multiple sentences from a file
def part1d(file):
    with open(file, "r") as fp:
        for line in fp:
            print(line, otherpart1c(line.strip()))
            

# this is the real part 1c
def otherpart1c(sentence):
    graph = []
    for word in sentence.split():
        tags = []
        for key in wt_probs:
            parts = key.split("//")
            if parts[0] == word:
                tags.append(key)
        graph.append(tags)
    

    probgraph = {}
    i = 0
    for word in sentence.split():
        
        if i == 0:
            for tags in graph[i]:
                parts = tags.split("//")
                prev = "...EOS..."
                probgraph[parts[1]] = wt_probs[tags] * ti_1_ti_probs[parts[1] + "//" + prev]
        else:
            currkeys = probgraph.keys()
            for key in list(currkeys):
                oldprob = probgraph[key]
                prevparts = key.split("//")
                prev = prevparts[len(prevparts) - 1]
                probgraph.pop(key)

                for tags in graph[i]:
                    parts = tags.split("//")
                    seq = key + "//" + parts[1]
                    if parts[1] + "//" + prev not in ti_1_ti_probs:
                        continue
                    probgraph[seq] = oldprob * wt_probs[tags] * ti_1_ti_probs[parts[1] + "//" + prev]

                

        i += 1

    maxprob = -1
    maxseq = ""
    for key in probgraph:
        if probgraph[key] > maxprob:
            maxprob = probgraph[key]
            maxseq = key
    
    maxseq = "...BOS...//" + maxseq + "//...EOS..."
    return maxseq

    
   


        

# Compute P(w|w_i-1,w_i-2)
def part2a(file):

    global wi_1_2_freqs
    wi_1_2_freqs = {}
    global _1_2_freqs
    _1_2_freqs = {}
    global wi_1_2_probs
    wi_1_2_probs = {}

    with open(file, "r") as fp:

        # w_{i-1} and w_{i-2}
        prev2 = ""
        prev1 = ""
        
        for line in fp:

            words = line.strip().split(" ")
            for word in words:
                altogether = tuple([prev2, prev1, word])
                justprevs = tuple([prev2, prev1])
                if altogether in wi_1_2_freqs:
                    wi_1_2_freqs[altogether] += 1
                else:
                    wi_1_2_freqs[altogether] = 1
                if justprevs in _1_2_freqs:
                    _1_2_freqs[justprevs] += 1
                else:
                    _1_2_freqs[justprevs] = 1
                
                # update prevs
                prev2 = prev1
                prev1 = word
            
    for key in wi_1_2_freqs:
        justprevs = tuple([key[0], key[1]])
        wi_1_2_probs[key] = wi_1_2_freqs[key] / _1_2_freqs[justprevs]


# Witten Bell Smoothing
# gets counts for uni,bi,trigrams
def part2b_1(file):
    
    global unicount
    unicount = {}
    global bicount
    bicount = {}
    global tricount
    tricount = {}

    with open(file, "r") as fp:

        prev2 = ""
        prev1 = ""

        for line in fp:
            words = line.strip().split(" ")
            for word in words:
                if word in unicount:
                    unicount[word] += 1
                else:
                    unicount[word] = 1
                if prev1 != "":
                    bigram = tuple([prev1, word])
                    if bigram in bicount:
                        bicount[bigram] += 1
                    else:
                        bicount[bigram] = 1
                    
                    if prev2 != "":
                        trigram = tuple([prev2, prev1, word])
                        if trigram in tricount:
                            tricount[trigram] += 1
                        else:
                            tricount[trigram] = 1
                
                prev2 = prev1
                prev1 = word
    
    global unilen
    unilen = len(unicount)
                
# does wb smoothing after getting counts in 2b_1
def part2b_2(i_2, i_1, i):
    trigram = tuple([i_2, i_1, i])
    bigram = tuple([i_2, i_1])
    bigram2 = tuple([i_1, i])
    T = 0
    N = 0
    
    # case: seen this whole trigram
    if trigram in tricount:
        for tri in tricount:
            bi = tuple([tri[0], tri[1]])
            if bi == bigram:
                T += 1
        return tricount[trigram] / (bicount[bigram] + T)
    # case: not seen trigram, but have seen i-2 i-1
    elif bigram in bicount:
        for tri in tricount:
            bi = tuple([tri[0], tri[1]])
            if bi == bigram:
                T += 1
        N = bicount[bigram]
        return T / ((unilen - T) * (N + T))
    # case: not seen bigram, have seen i-1 i
    elif bigram2 in bicount:
        for bi in bicount:
            uni = bi[0]
            if uni == i_1:
                T += 1
        return bicount[bigram2] / (unicount[i_1] + T)
    # case: have seen i-1 
    elif i_1 in unicount:
        for bi in bicount:
            uni = bi[0]
            if uni == i_1:
                T += 1
        return T / ((unilen - T) * (unicount[i_1] + T))
    # case: have seen i
    elif i in unicount:
        return unicount[i] / unilen
    # case: have not seen i
    else:
        uniN = 0
        for uni in unicount:
            uniN += unicount[uni]
        return unilen / (unilen * (uniN + unilen))

# compute probability of a sentence
def part2c(sentence):
    words = sentence.split()
    totalprob = 1
    for i in range(2, len(words)):
        totalprob = totalprob * part2b_2(words[i - 2], words[i - 1], words[i])
    return totalprob

# compute perplexities of sentences in a file
def part2d(file):
    with open(file, "r") as fp:
        for line in fp:
            n = len(line.strip().split())
            exp = (-1 / n)
            prob = part2c(line.strip())
            perp = math.pow(prob, exp)
            print(line, perp)


    

                



# For Part 1, run 1a -> 1b -> 1d
# For Part 2, run 2b_1 -> 2d
part1a("corpus-2.1.txt")
part1b("corpus-2.1.txt")
part1d("2-1.txt")
#part2a("corpus-2.2.txt")
part2b_1("corpus-2.2.txt")
part2d("test-2.2.txt")      # calls 2c which calls 2b_2








"""
def part1c(sentence):
    taggedsentence = ""
    prev = "...EOS..."
    for word in sentence.split():
        taggedsentence += word + "//"
        maxprob = 0
        maxtag = ""
        for key in wt_probs:
            parts = key.split("//")
            
            if parts[0] != word:
                continue

            for key2 in ti_1_ti_probs:
                parts2 = key2.split("//")
                curr_ti_ti1 = parts2[0] + "//" + prev

                if curr_ti_ti1 not in ti_1_ti_probs:
                    continue
                
                if wt_probs[key] * ti_1_ti_probs[key2] > maxprob:
                    maxprob = wt_probs[key] * ti_1_ti_probs[key2]
                    maxtag = parts[1]

        taggedsentence += maxtag + " "
        print(maxtag, maxprob)
        prev = maxtag
    
    return taggedsentence
"""




"""
# Compute P(w|w_i-1,w_i-2)
def part2a(file):

    global wi_1_2_freqs
    wi_1_2_freqs = {}
    global _1_2_freqs
    _1_2_freqs = {}
    global wi_1_2_probs
    wi_1_2_probs = {}

    with open(file, "r") as fp:

        # w_{i-1} and w_{i-2}
        prev2 = ""
        prev1 = ""
        
        for line in fp:

            words = line.strip().split(" ")
            for word in words:
                altogether = tuple([prev2, prev1, word])
                justprevs = tuple([prev2, prev1])
                if altogether in wi_1_2_freqs:
                    wi_1_2_freqs[altogether] += 1
                else:
                    wi_1_2_freqs[altogether] = 1
                if justprevs in _1_2_freqs:
                    _1_2_freqs[justprevs] += 1
                else:
                    _1_2_freqs[justprevs] = 1
                
                # update prevs
                prev2 = prev1
                prev1 = word
            
    for key in wi_1_2_freqs:
        justprevs = tuple([key[0], key[1]])
        wi_1_2_probs[key] = wi_1_2_freqs[key] / _1_2_freqs[justprevs]


# Witten Bell Smoothing
# gets counts for uni,bi,trigrams
def part2b_1(file):
    
    global unicount
    unicount = {}
    global bicount
    bicount = {}
    global tricount
    tricount = {}

    with open(file, "r") as fp:

        prev2 = ""
        prev1 = ""

        for line in fp:
            words = line.strip().split(" ")
            for word in words:
                if word in unicount:
                    unicount[word] += 1
                else:
                    unicount[word] = 1
                if prev1 != "":
                    bigram = tuple([prev1, word])
                    if bigram in bicount:
                        bicount[bigram] += 1
                    else:
                        bicount[bigram] = 1
                    
                    if prev2 != "":
                        trigram = tuple([prev2, prev1, word])
                        if trigram in tricount:
                            tricount[trigram] += 1
                        else:
                            tricount[trigram] = 1
                
                prev2 = prev1
                prev1 = word
                
# does wb smoothing after getting counts in 2b_1
def part2b_2(i_2, i_1, i):
    trigram = tuple([i_2, i_1, i])
    bigram = tuple([i_2, i_1])
    bigram2 = tuple([i_1, i])
    T = 0
    N = 0
    unilen = len(unicount)
    # case: seen this whole trigram
    if trigram in tricount:
        for tri in tricount:
            bi = tuple([tri[0], tri[1]])
            if bi == bigram:
                T += 1
        return tricount[trigram] / (bicount[bigram] + T)
    # case: not seen trigram, but have seen i-2 i-1
    elif bigram in bicount:
        for tri in tricount:
            bi = tuple([tri[0], tri[1]])
            if bi == bigram:
                T += 1
        N = bicount[bigram]
        return T / ((unilen - T) * (N + T))
    # case: not seen bigram, have seen i-1 i
    elif bigram2 in bicount:
        for bi in bicount:
            uni = bi[0]
            if uni == i_1:
                T += 1
        return bicount[bigram2] / (unicount[i_1] + T)
    # case: have seen i-1 
    elif i_1 in unicount:
        for bi in bicount:
            uni = bi[0]
            if uni == i_1:
                T += 1
        return T / ((unilen - T) * (unicount[i_1] + T))
    # case: have seen i
    elif i in unicount:
        for uni in unicount:
            N += unicount[uni]
        return unicount[i] / unilen
    # case: have not seen i
    else:
        for uni in unicount:
            N += unicount[uni]
        return unilen / (unilen * (N + unilen))

"""