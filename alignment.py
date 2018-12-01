from testSequences import occurences
from genesList import *

gap_penality = -1
match_award = 1
mismatch_penality = -1


def createMatrix(seq1, seq2):
    matrix = [[0]* (seq2)]*(seq1)
    return matrix

# def print_matrix(matrix):
#     for r in matrix:
#         print( *r , sep=" ")

#initializes the scoring matrix with all zeros
def zeros(rows, cols):
    matrix = []
    for _ in range(rows):
        matrix.append([])
        for _ in range(cols):
            matrix[-1].append(0)

    return matrix

#the gives a score to the nucleotides that are being compared
def scoring(a, b):
    if a == b:
        return match_award
    elif a =='-' or b == '-':
        return gap_penality
    else:
        return mismatch_penality

# the alignment algorithm
# seq1 is the original sequence
# seq2 is the mutated one
def needleman_wunsch(seq1, seq2):
    #in case either string is uppercase/lowercase. 
    #this makes sure that we are comparing the strings in one case to avoid any unncessary errors while comparing
    seq1 = seq1.upper()
    seq1 = seq1[20:]
    seq2 = seq2.upper()

    #taking the length of the 2 sequences 
    n = len(seq1)
    m = len(seq2)
    #create a scoring matrix
    score = zeros(m + 1, n + 1)

    #filling in th scoring table
    for i in range(0, m + 1):
        score[i][0] = gap_penality * i
    for j in range(0, n + 1):
        score[0][j] = gap_penality * j
    
    #comparing the 2 sequences and updating the scoring table 
    #if the 2 match there is one score given 
    #if the 2 doesnt there is another score given bases on gap penality
    for i in range(1, m +1):
        for j in range( 1, n + 1):
            match = score[i-1][j-1] + scoring(seq1[j-1],seq2[i-1])
            delete = score[i-1][j] + gap_penality   
            insert = score[i][j-1]  + gap_penality 
            score[i][j] = max(match, delete, insert) 

    #with the score table complete
    #it is now used to backtrack from the bottom right back to the top left to get 
    #the final optimal sequence
    
    #create 2 empty string to store the alignments
    align1 = ""
    align2 = ""

    i = m
    j = n
    
    #this while condition will go back up to the top left of the matrix
    while i > 0 and j > 0:
        # storing the differetn scores to know how the current score got there
        # ie which combination it was
        score_current = score[i][j]
        score_diagonal = score[i-1][j-1]
        score_up = score[i][j-1]
        score_left = score[i-1][j]

        if score_current == score_diagonal + scoring(seq1[j-1], seq2[i-1]):
            align1 += seq1[j-1]
            align2 += seq2[i-1]
            j -= 1
            i -= 1
        elif score_current == score_up + gap_penality:
            align1 += seq1[j-1]
            align2 += "-"
            j -= 1
        elif score_current == score_left + gap_penality:
            align1 += "-"
            align2 += seq2[i-1]
            i -= 1

    #tracing up to top left
    while j > 0:
        align1 += seq1[j-1]
        align2 += "-"
        j -= 1

    while i > 0:
        align1 += "-"
        align2 += seq2[i-1]
        i -= 1

    align1 = align1[::-1]
    align2 = align2[::-1]
    
    return(align1, align2)

def positions_of_diff(aligned_str):
    diff_position = []
    for x in range(len(aligned_str)):
        if aligned_str[x] == '-':
            diff_position.append(x+1)
    
    return diff_position

#this check for mutations method go through the array that stores the positions where
#misalignment happens and returns the 20bp previous to the mutation
def get_guiding_strand(aligned_str, diff_position):
    for i in range(len(diff_position)):
        # if(i >= aligned_str.gene_location_start and i < aligned_str.gene_location_end):
        print(aligned_str)
        return (aligned_str[diff_position[i]-20:diff_position[i]] )
        #,occurences(aligned_str[i-20:i], complete_genome))

def taking_input_and_run_task():
    sequences_list = {}
    range_num = int(input("State the number of genes you'd like to align(ex 1, 3, 10): "))
    if range_num <= 0:
        print("ERROR RANGE CAN NOT BE LESS THAN 1")
        
    for _ in range(range_num):
        key = input("Enter the name of sequences(ex inhA): ")
        value = raw_input("Enter the sequence: ")

        print("key: {}".format(key))
        print("value: {} ".format(value))

        sequences_list[key] = value

    aligning(sequences_list)

def aligning(sequences_list):
    for k, v in sequences_list.items():
        _, aligned = needleman_wunsch(k, v)
        print("ALIGNED: {}".format(aligned))
        diff_position = positions_of_diff(aligned)
        get_guiding_strand(k ,diff_position)
        print("aligning them")


def create_complete_genome():
    builder = []
    with open("m_tuberculosis_complete_genome.txt",  'r') as work_data:    
        for line in work_data:
            builder.append(line.replace("\n",""))
        
    return ''.join(builder)

def main():
    taking_input_and_run_task()

if __name__ == '__main__': 
    main()

    # gene = "tgggacgggttaacctcaaagcc"
    # gene_mutated = "tgctacggggaacctctaagcc"
    # _, aligned = needleman_wunsch(gene, gene_mutated)
    # mutation_positions_list = positions_of_diff(aligned)

    # print("There are mutations at positions: ") 
    # print(mutation_positions_list)
    # print("The leading 20bp sequence to fix the mutation is " + get_guiding_strand(aligned, mutation_positions_list))





    #print(needleman_wunsch(gene, gene_mutated))
    # taking_input();  
    
    # seq1 = "asdfasssdfaadadfsasaddafadfsafaasaaaadsafsadfsdfsdfsssd"
    # seq2 = "asdfassfadadfsasaddafadfsafaasaaaadsafsadfsdfsdfsaaassd"
   
    # _, aligned = needleman_wunsch(seq1, seq2)
    # print(seq1.upper())
    # print(aligned)
    # print(positions_of_diff(aligned))
    # createMatrix(2, 3)
    # zeros(2,3)'
