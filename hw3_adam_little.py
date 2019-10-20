# Python version == 3.7.2

import math

# Define our three sequences
s1 = 'ATTDEWKKQRKDSHKEVERRRRENINTAINVLSDLLPVRESSKAAILACAAEYIQKLKETDEANIEKWTLQKLLSEQNASQLASANEKLQEELGNAYKEIEYMKRVLRK----------'
s2 = 'HGSEEWHRQRRENHKEVERKRRESINTGIRELARLIPTTDTNKAQILQRAVEYIKRLKENENNNIEKWTLEKLLTEQAVSELSASNEKLKHELESAYREIEQLKRGKK-----------'
s3 = 'TGSTAWKQQRKESHKEVERRRRQNINTAIEKLSDLLPVKETSKAAILSRAAEYIQKMKETETANIEKWTLQKLLGEQQVSSLTSANDKLEQELSKAYKNLQELKKKLKEAGIEDPTEEE'

# Blosum62 substitution matrix
blosum62 = {
    'A':{'C':0,  'S':1,  'T':-1, 'P':-1, 'A':4,  'G':0,  'N':-1, 'D':-2, 'E':-1, 'Q':-1, 'H':-2, 'R':-1, 'K':-1, 'M':-1, 'I':-1, 'L':-1, 'V':-2, 'F':-2, 'Y':-2, 'W':-3},
    'G':{'C':-3, 'S':0,  'T':1,  'P':-2, 'A':0,  'G':6,  'N':-2, 'D':-1, 'E':-2, 'Q':-2, 'H':-2, 'R':-2, 'K':-2, 'M':-3, 'I':-4, 'L':-4, 'V':0,  'F':-3, 'Y':-3, 'W':-2},
    'I':{'C':-1, 'S':-2, 'T':-2, 'P':-3, 'A':-1, 'G':-4, 'N':-3, 'D':-3, 'E':-3, 'Q':-3, 'H':-3, 'R':-3, 'K':-3, 'M':1,  'I':4,  'L':2,  'V':1,  'F':0,  'Y':-1, 'W':-3},
    'L':{'C':-1, 'S':-2, 'T':-2, 'P':-3, 'A':-1, 'G':-4, 'N':-3, 'D':-4, 'E':-3, 'Q':-2, 'H':-3, 'R':-2, 'K':-2, 'M':2,  'I':2,  'L':4,  'V':3,  'F':0,  'Y':-1, 'W':-2},
    'P':{'C':-3, 'S':-1, 'T':1,  'P':7,  'A':-1, 'G':-2, 'N':-1, 'D':-1, 'E':-1, 'Q':-1, 'H':-2, 'R':-2, 'K':-1, 'M':-2, 'I':-3, 'L':-3, 'V':-2, 'F':-4, 'Y':-3, 'W':-4},
    'V':{'C':-1, 'S':-2, 'T':-2, 'P':-2, 'A':0,  'G':-3, 'N':-3, 'D':-3, 'E':-2, 'Q':-2, 'H':-3, 'R':-3, 'K':-2, 'M':1,  'I':3,  'L':1,  'V':4,  'F':-1, 'Y':-1, 'W':-3},
    'F':{'C':-2, 'S':-2, 'T':-2, 'P':-4, 'A':-2, 'G':-3, 'N':-3, 'D':-3, 'E':-3, 'Q':-3, 'H':-1, 'R':-3, 'K':-3, 'M':0,  'I':0,  'L':0,  'V':-1, 'F':6,  'Y':3,  'W':1},
    'W':{'C':-2, 'S':-3, 'T':-3, 'P':-4, 'A':-3, 'G':-2, 'N':-4, 'D':-4, 'E':-3, 'Q':-2, 'H':-2, 'R':-3, 'K':-3, 'M':-1, 'I':-3, 'L':-2, 'V':-3, 'F':1,  'Y':2,  'W':11},
    'Y':{'C':-2, 'S':-2, 'T':-2, 'P':-3, 'A':-2, 'G':-3, 'N':-2, 'D':-3, 'E':-2, 'Q':-1, 'H':2,  'R':-2, 'K':-2, 'M':-1, 'I':-1, 'L':-1, 'V':-1, 'F':3,  'Y':7,  'W':2},
    'D':{'C':-3, 'S':0,  'T':1,  'P':-1, 'A':-2, 'G':-1, 'N':1,  'D':6,  'E':2,  'Q':0,  'H':-1, 'R':-2, 'K':-1, 'M':-3, 'I':-3, 'L':-4, 'V':-3, 'F':-3, 'Y':-3, 'W':-4},
    'E':{'C':-4, 'S':0,  'T':0,  'P':-1, 'A':-1, 'G':-2, 'N':0,  'D':2,  'E':5,  'Q':2,  'H':0,  'R':0,  'K':1,  'M':-2, 'I':-3, 'L':-3, 'V':-3, 'F':-3, 'Y':-2, 'W':-3},
    'R':{'C':-3, 'S':-1, 'T':-1, 'P':-2, 'A':-1, 'G':-2, 'N':0,  'D':-2, 'E':0,  'Q':1,  'H':0,  'R':5,  'K':2,  'M':-1, 'I':-3, 'L':-2, 'V':-3, 'F':-3, 'Y':-2, 'W':-3},
    'H':{'C':-3, 'S':-1, 'T':0,  'P':-2, 'A':-2, 'G':-2, 'N':1,  'D':1,  'E':0,  'Q':0,  'H':8,  'R':0,  'K':-1, 'M':-2, 'I':-3, 'L':-3, 'V':-2, 'F':-1, 'Y':2,  'W':-2},
    'K':{'C':-3, 'S':0,  'T':0,  'P':-1, 'A':-1, 'G':-2, 'N':0,  'D':-1, 'E':1,  'Q':1,  'H':-1, 'R':2,  'K':5,  'M':-1, 'I':-3, 'L':-2, 'V':-3, 'F':-3, 'Y':-2, 'W':-3},
    'S':{'C':-1, 'S':4,  'T':1,  'P':-1, 'A':1,  'G':0,  'N':1,  'D':0,  'E':0,  'Q':0,  'H':-1, 'R':-1, 'K':0,  'M':-1, 'I':-2, 'L':-2, 'V':-2, 'F':-2, 'Y':-2, 'W':-3},
    'T':{'C':-1, 'S':1,  'T':4,  'P':1,  'A':-1, 'G':1,  'N':0,  'D':1,  'E':0,  'Q':0,  'H':0,  'R':-1, 'K':0,  'M':-1, 'I':-2, 'L':-2, 'V':-2, 'F':-2, 'Y':-2, 'W':-3},
    'C':{'C':9,  'S':-1, 'T':-1, 'P':-3, 'A':0,  'G':-3, 'N':-3, 'D':-3, 'E':-4, 'Q':-3, 'H':-3, 'R':-3, 'K':-3, 'M':-1, 'I':-1, 'L':-1, 'V':-1, 'F':-2, 'Y':-2, 'W':-2},
    'M':{'C':-1, 'S':-1, 'T':-1, 'P':-2, 'A':-1, 'G':-3, 'N':-2, 'D':-3, 'E':-2, 'Q':0,  'H':-2, 'R':-1, 'K':-1, 'M':5,  'I':1,  'L':2,  'V':-2, 'F':0,  'Y':-1, 'W':-1},
    'N':{'C':-3, 'S':1,  'T':0,  'P':-2, 'A':-2, 'G':0,  'N':6,  'D':1,  'E':0,  'Q':0,  'H':-1, 'R':0,  'K':0,  'M':-2, 'I':-3, 'L':-3, 'V':-3, 'F':-3, 'Y':-2, 'W':-4},
    'Q':{'C':-3, 'S':0,  'T':0,  'P':-1, 'A':-1, 'G':-2, 'N':0,  'D':0,  'E':2,  'Q':5,  'H':0,  'R':1,  'K':1,  'M':0,  'I':-3, 'L':-2, 'V':-2, 'F':-3, 'Y':-1, 'W':-2},
}


# Function for log 2
def log2(base):
    if base != 0:
        log_2 = math.log(base, 2.0)
        return log_2
    # Need to do this to account for possible domain error if base = 0
    else:
        return 0

# Finding the probability distribution value
def p(aa1, aa2, aa3):

    p1 = 0
    p2 = 0
    p3 = 0
    # Compare the three bases in the sequence
    # If a value occurs once, give it a probability of 1/20

    # Skip values with a -
    if aa1 == aa2 and aa1 == aa3:
        p1 = 1
        p2 = 0
        p3 = 0
    elif aa1 == aa2 and aa2 != aa3:
        p1 = 2/3
        p2 = 1/3
        p3 = 0
    elif aa1 != aa2 and aa2 == aa3:
        p1 = 1/3
        p2 = 2/3
        p3 = 0
    elif aa1 == aa3 and aa3 != aa2:
        p1 = 2/3
        p2 = 1/3
        p3 = 0
    elif aa1 == aa3 and aa1 != aa2:
        p1 = 2/3
        p2 = 1/3
        p3 = 0
    elif aa1 != aa2 and aa2 != aa3 and aa1 != aa3:
        p1 = 1/3
        p2 = 1/3
        p3 = 1/3
    else:
        print("Something is wrong!")

    return p1, p2, p3

def entropy(seq1, seq2, seq3):
    entropy_score = 0
    for aa in range(0, len(seq1)):\
        # Break if the sequence has a dash
        if seq1[aa] == '-' or seq2[aa] == '-' or seq3[aa] == '-':
            break
        p1, p2, p3 = p(seq1[aa], seq2[aa], seq3[aa])
        marginal_score = -(p1*log2(p1) \
                         + p2*log2(p2) \
                         + p3*log2(p3))
        # print("Entropy score for {}{}{}: {}".format(seq1[aa], seq2[aa], seq3[aa], marginal_score))
        entropy_score += marginal_score
        # print("Total entropy score: {}". format(entropy_score))

    return entropy_score

entropy_score = entropy(s1, s2, s3)
print("Entropy score: {}".format(entropy_score))

def sum_of_pair(seq1, seq2, seq3):
    total_sum = 0
    for aa in range(0, len(seq1)):
        if seq1[aa] == '-' or seq2[aa] == '-' or seq3[aa] == '-':
            break
        # print("Column contains: {}{}{}".format(seq1[aa], seq2[aa], seq3[aa]))
        # Pairs of first and second sequences
        pair_1_2 = blosum62[seq1[aa]][seq2[aa]]
        # print("First pair: {}".format(pair_1_2))
        # Pair of second and third sequences
        pair_2_3 = blosum62[seq2[aa]][seq3[aa]]
        # print("Second pair: {}".format(pair_2_3))
        # Pair of first and third sequences
        pair_1_3 = blosum62[seq1[aa]][seq3[aa]]
        # print("Third pair: {}".format(pair_1_3))

        sum = pair_1_2 + pair_1_3 + pair_2_3

        total_sum += sum
        # print("Total sum: {}".format(total_sum))

    return total_sum



sp_score = sum_of_pair(s1, s2, s3)
print("SP-Score: {}".format(sp_score))
