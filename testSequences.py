import random

def occurences(gene, sequence):
    index = 0
    count = 0
    length = len(gene)

    while(index != -1):
        index = sequence.find(gene, index)

        if(index != -1):
            count+=1
            index += length
    
    return count

def generate_test_seqs():
    nucleotides = ["a","c","t","g"]
    range_num = input("enter size of gene: ")
    # gene =""
    seq = ""
    # for x in range(7):
    #    gene += nucleotides[random.randint(0,3)]
    
    for x in range(range_num):
       seq += nucleotides[random.randint(0,3)]
    

    return seq
    
# def return substring(str):
#     for i in range(26):
        
        
if __name__ == "__main__":

    print(generate_test_seqs())
    g = 'GATATATGCATATACTT'
    q = 'ATAT'
    i = 23
    print(occurences(q,g))
    print(g[i-20:i])
    