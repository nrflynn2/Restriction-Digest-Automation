# Cutsultant is a quick script designed to generate the most optimal cuts for a restriction enzyme given a list of parameters
# Note: argument_handler() will overwrite the .fasta file - be aware of this.
# Note: BioPython libraries (and subsequently this code) only produce accurate results for non-methylated DNA. 

import sys
import time
import itertools
from Bio.Restriction import *
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

def argument_handler(arg):
    with open(arg, 'r') as fin:
        data = fin.read().splitlines(True)
    with open(arg, 'w') as fout:
        fout.writelines(data[1:])
    with open(arg, "r") as myfile:
        ret_string = myfile.read().replace('\n', '')
    return ret_string    

gene_sequence_one = Seq(argument_handler(sys.argv[2]))
if(len(sys.argv) == 4): #If there are 2 sequences given
    gene_sequence_two = Seq(argument_handler(sys.argv[3]))    

enzyme_list = [line.strip() for line in open(sys.argv[1])] #slick
IUPAC = IUPACAmbiguousDNA()    
smallestSize = input('What is the smallest band you can measure? ') 
largestSize = input('What is the large band you can measure? ')
resolution = input('What is the ratio of band lengths you require? ')
max_enzymes = input('What is the maximum number of enzymes you wish you use? ')
start_time = time.time() #Start timer for bragging after all inputs
restrictionbatch = RestrictionBatch(enzyme_list)

def enzymeEliminator(dict1, dict2, rb):

    """Takes in two dictionaries which map restriction enzymes to number of times it cuts and determines whether they cut the same number or not"""
    #Key: Enzyme
    #Value: List containing cut sites
    #This function iterates through a dictionary and removes enzymes that cut the same in both inserts
    for key in dict2: 
        if ( dict1.has_key(key) and dict2.has_key(key) and (dict1[key] == dict2[key]) ): 
            enzyme_list.remove(str(key)) 
            restrictionbatch.remove(key) 

def number_of_cuts(comb_tup, sequence):
    """This function sees how many cuts a given combination of enzymes has, takes in a combination tuple and a sequence"""

    #The comb_tup is what is generated from itertools.combinations - it is an element from the list of tuples
    #The tuples are the various enzyme combinations possible - this is getting iterated in the main function
    
    comb_list = list(comb_tup) #Cast to a list because RestrictionBatch requires a list as an argument
    temp_rb = RestrictionBatch(comb_list) #Make a RestrictionBatch out of our tuple
    temp_dict = temp_rb.search(sequence)  #Create a dictionary that maps enzymes->cut sites
    count = 0 #Set a counter to 0. This will count how many times our plasmid is getting cut 
    
    for key in temp_dict: #For every key inside the dictionary
        for x in temp_dict[key]: #And for every element in every list the key maps to
            count += 1 #Increase the count by 1. 
    return count

def cut_to_size_converter(comb_tup, sequence):
    """Converts cut sites into band lengths using combination tuples and a sequence as an input"""

    comb_list = list(comb_tup) #Cast to a list because RestrictionBatch requires a list as an argument
    temp_rb = RestrictionBatch(comb_list) #Same thing as in number_of_cuts()
    temp_dict = temp_rb.search(sequence) #Same thing as in number_of_cuts()
    cut_site_list = [] #We are starting this one as empty, the circular nature of the plasmid will be incorporated later
    band_list = [] #Make an empty list with which you can throw in band lengths in later

    for key in temp_dict: #For a key in the dictionary
        for x in temp_dict[key]: #And within that list
            cut_site_list.append(x) #Append all the cut sites to it

    cut_site_list = sorted(cut_site_list) #Then sort it, so it goes [0, cutsite1, cutsite2, ..., len(sequence)]
    for x in range(0, len(cut_site_list)-1): #Should this be -1? Check later and how for loop bounds work in python
        band_list.append( cut_site_list[x+1] - cut_site_list[x] ) #Append the difference - i.e., the band length

    if(len(cut_site_list) != 0):
        band_list.append(len(sequence) - cut_site_list[len(cut_site_list) - 1] + cut_site_list[0]) #This makes it circular    

    return sorted(band_list) #Return it for other functions to use

def differentiability_determiner(comb_tup, sequence): #The bigger band should be 20% bigger - instead of semilog
    """Returns a boolean value saying if the inputted tuple demonstrates sufficient differentiability on the inputted sequence"""
    
    #TODO: For non-plebbyness, check if the line 2 lines below this one is needed
    bands_list = cut_to_size_converter(comb_tup, sequence) #Uses the band_list that cut_to_size_converter() returned

    #These two print lines can be used to test if the band lengths are what they ought to be
    #print comb_tup
    #print bands_list

    for y in range(0, len(bands_list)-2): #The -2 is done here

        if(bands_list[y] == 0): #Does this fix the n=31 problem?
            continue
        elif((float(bands_list[y+1])/bands_list[y]) <= resolution): #Division is an integer calculation normally, you need to cast to a float......
            return False    
    return True             

def size_band_eliminator(comb_tup, sequence):
    """Tells us if there are bands that are < min band size or > max band size"""

    bands_list = cut_to_size_converter(comb_tup, sequence)

    for x in range(0, len(bands_list)): #WE CHANGED THE BOUNDS
        if ((bands_list[x] < smallestSize) or (bands_list[x] > largestSize)):
            return False
    return True    

def twoPlasmidAnalysis(restrictionbatch, gene_sequence_two, gene_sequence_one):
    """Finds the ideal restriciton enzymes when examining two plasmids"""

    plasmid_WithoutInsert_analysis = Analysis(restrictionbatch, gene_sequence_one)
    plasmid_WithInsert_analysis = Analysis(restrictionbatch, gene_sequence_two)

    #This eliminates all the enzymes that cut the same in both plasmids
    temp_dict1 = plasmid_WithoutInsert_analysis.full()
    temp_dict2 = plasmid_WithInsert_analysis.full()
    enzymeEliminator(temp_dict1, temp_dict2, restrictionbatch)

    #TODO: See if we can get away with deleting the bottom line
    restrictionbatch2 = RestrictionBatch(enzyme_list) #Updates the restrictionbatch just in case our code didn't before

    final_enzymes = []
    for n in range(1, max_enzymes+1):
        combination_list = list(itertools.combinations(enzyme_list, n))

        for x in combination_list:
            if((differentiability_determiner(x, gene_sequence_one) and differentiability_determiner(x, gene_sequence_two)) and \
            (size_band_eliminator(x, gene_sequence_two) and size_band_eliminator(x, gene_sequence_one)) ):
            
                final_enzymes.append(x)

    #Outputs the three sets of enzymes that produce the most cuts    
    first_place = 0    
    second_place = 0
    third_place = 0

    first_pick = None
    second_pick = None
    third_pick = None

    for x in final_enzymes:
        total_num_sites = number_of_cuts(x, gene_sequence_one) + number_of_cuts(x, gene_sequence_two)

        if(total_num_sites > first_place):
            third_place = second_place
            third_pick = second_pick
            second_place = first_place
            second_pick = first_pick
            first_place = total_num_sites
            first_pick = x
        elif(total_num_sites > second_place):
            third_place = second_place
            third_pick = second_pick
            second_place = total_num_sites
            second_pick = x
        elif(total_num_sites > third_place):
            third_place = total_num_sites
            third_pick = x

    count = 0
    min_enzyme_list = []
    for x in final_enzymes:
        if(count < 3):
            min_enzyme_list.append(x)
            count += 1

    print("\n")
    print("The ideal set of restriction enzymes to digest this plasmid to produce maximum cuts are: {}, {}, and {}" .format(first_pick, second_pick, third_pick))
    print("\n")
    print("The ideal set of restriction enzymes to digest this plasmid in the minimum number of enzymes are:")
    for x in range(0, len(min_enzyme_list)):
        print(min_enzyme_list[x],)
    print("\n")
    print("The entire set of restriction enzymes to digest this plasmid based on the parameters given are:") 
    print(final_enzymes)
    print("\n")
    
def onePlasmidAnalysis(restrictionbatch, gene_sequence_one): #Make sure everything in this function is gene_sequence_one
    """Finds the ideal restriciton enzymes when digesting two plasmids"""

    #Note: For testing, this should all be done on the plasmid with the GFP insert

    plasmid_WithInsert_analysis = Analysis(restrictionbatch, gene_sequence_one) #Same set up as above
    temp_dict = plasmid_WithInsert_analysis.full()

    temp_dict2 = plasmid_WithInsert_analysis.with_N_sites(0, dct = None) 
    for key in temp_dict2: 
        enzyme_list.remove(str(key)) 
        restrictionbatch.remove(key) 
#    for key in temp_dict: #Removes everything that doesn't cut
#        if(len(temp_dict[key]) == 0):
#            enzyme_list.remove(str(key))

    restrictionbatch2 = RestrictionBatch(enzyme_list) #Updates the restrictionbatch just in case our code didn't before

    final_enzymes = []
    for n in range(1, max_enzymes+1): #TODO: Fix the range thing, make it user input perhaps?
        combination_list = list(itertools.combinations(enzyme_list, n))
        for x in combination_list:
            if(differentiability_determiner(x, gene_sequence_one)) and \
            (size_band_eliminator(x, gene_sequence_one)):
                
                final_enzymes.append(x)

    #TODO: Wrap the above and the bottom in a helper function
    
    first_place = 0    
    second_place = 0
    third_place = 0

    first_pick = None
    second_pick = None
    third_pick = None

    #Outputs the three sets of enzymes that produce the most cuts
    for x in final_enzymes:
        total_num_sites = number_of_cuts(x, gene_sequence_one)

        if(total_num_sites > first_place):
            third_place = second_place
            third_pick = second_pick
            second_place = first_place
            second_pick = first_pick
            first_place = total_num_sites
            first_pick = x
        elif(total_num_sites > second_place):
            third_place = second_place
            third_pick = second_pick
            second_place = total_num_sites
            second_pick = x
        elif(total_num_sites > third_place):
            third_place = total_num_sites
            third_pick = x            

    count = 0
    min_enzyme_list = []
    for x in final_enzymes:
        if(count < 3):
            min_enzyme_list.append(x)
            count += 1

    print("\n")
    print("The ideal set of restriction enzymes to digest this plasmid to produce maximum cuts are: {}, {}, and {}" .format(first_pick, second_pick, third_pick))
    print("\n")
    print("The ideal set of restriction enzymes to digest this plasmid in the minimum number of enzymes are: ")
    for x in range(0, len(min_enzyme_list)): print(min_enzyme_list[x],)
    print("\n")
    print("The entire set of restriction enzymes to digest this plasmid based on the parameters given are:")
    print(final_enzymes)
    print("\n")

if(len(sys.argv) == 3):
    onePlasmidAnalysis(restrictionbatch, gene_sequence_one)
elif(len(sys.argv) == 4):
    twoPlasmidAnalysis(restrictionbatch, gene_sequence_two, gene_sequence_one)
else:
    print("The command line input was incorrect, please make it >>python Cutsultant.py enzymes.txt seq1.fasta and if applicable add seq2.fasta after that")

print("Time to execute:", time.time() - start_time, "seconds")
