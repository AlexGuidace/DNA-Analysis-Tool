# DNA Analysis Tool

# A program that compares a CSV table full of peoples' DNA STRs (Short Tandem Repeats), to a file 
# containing a single sequence of DNA, to determine who in the table that sequence belongs to.

# See https://cs50.harvard.edu/x/2020/psets/6/dna/ for a background on DNA profiling.

import sys

def main():
    # Checking to make sure the user puts in the correct number of inputs in the terminal.
    if len(sys.argv) != 3:
        print("Please input two files to compare.")
        return 0

    # Labeling the files to compare. 
    # csv_arg is the database of people and their STRs. DNA_arg is the DNA sequence file.
    csv_arg = sys.argv[1]
    DNA_arg = sys.argv[2]

    # A count of how many times every STR repeats within this DNA sequence.
    # I.e., STR = repetitions of the STR "AGAT": AGATAGATAGATAGATAGAT
    analyzed_sequence = read_and_analyze_DNA(DNA_arg, csv_arg)
    # Comparison used in determining who the DNA belongs to, if any.
    compare(csv_arg, analyzed_sequence)


# Read in and count repetitions of each STR in the DNA sequence.
def read_and_analyze_DNA(DNA_arg, csv_arg):
    with open(DNA_arg, 'r') as DNA:
        sequence = DNA.read() 

    analyzed_sequence = analyze_DNA(sequence, csv_arg)

    # Return results for comparison.
    return analyzed_sequence


# Returns an ordered list of STRs with their highest repetitions (AKA, their longest consecutive sequences).
# We want the highest repetitions of each STR because these are the repetitions we will be comparing our STR 
# csv database file to. If a person's values matches these values, we will have a DNA match.
def analyze_DNA(sequence, csv_arg):
    # Initialize a list to store our STRs.
    STR_list = []
    # Initialize a list to hold the longest repetition chains per STR analyzed.
    analyzed_sequence = []

    # Open our csv dictionary, or "database" to read it.
    dictionary = open(csv_arg, "r")

    # Populate a row-list of STR values from csv database file.
    for row in dictionary:
        # For the first row only, concerning the headers of each STR type.
        if row[0][0] == "n":
            # Split the header row of the csv_arg by comma values.
            STR_list = row.split(",")
            # Remove "name" from the row since we only want STR values.
            STR_list.remove("name")
            # Remove last element from list so you can remove the \n that won't go away any other way.
            last_element = STR_list[len(STR_list) - 1]
            STR_list.remove(last_element)
            # Initialize a new value to put what will be our truncated last_element in.
            new_last_element = ""
            # Iterate through each letter in the string and add each letter to our new string up until
            # (not including) the second to last character.
            for c in last_element:
                if len(new_last_element) < len(last_element) - 1:
                    new_last_element += c
            # Finally, add the new element that does not include \n to the STR_list.
            # These STR labels will be used to identify an STR in a DNA sequence file, and then
            # find its highest repetition.
            STR_list.append(new_last_element)

    # Analyze the DNA file (sequence) to find each STR's highest repetition.
    next_letter = True
    while next_letter:
        # We need to analyze every STR in the DNA file, so we loop through each STR in the STR_list.
        for STR in STR_list:
            # Get the length of the STR to identify and narrow down which STR we will be analyzing.
            STR_length = len(STR)
            # Initialize counters for each STR.
            current_count = 0
            largest_count = 0
            # Summary of for-loops below:

            # Using the len() function, slice format s[i:j], variables letter, letter_chunk, and letter_b, we can detect 
            # if we have a correct STR chunk or not and we can count repetitions of that STR if we do.

            # A concrete example:
            
            # Take the STR "AGATC". The length (len) of this STR is 5. But, since strings are 0-indexed, we actually only 
            # have 4 positions in the string, instead of 5.

            # Now, consider the sequence "AGATCAGATCAGATC": this is the beginning of a sequence in a hypothetical DNA file. 
            # "A" is going to be the variable "letter" in the first for-loop below. letter_chunk starts at "A" (index 0) and 
            # progresses through (index 0 + the length (len) of our STR (5)): that is, sequence[0: 0 + 5]. Since slices are 
            # non-inclusive regarding the last index, the "fifth element" is not going to be included (in 0 + 5). Instead, 
            # indexes 0 to 4 will be the elements included. In this way, we discover that "AGATC" is the STR we are going 
            # to be addressing.

            # So we now have an STR chunk, "AGATC", and we can look for other instances of this STR chunk in the DNA file. 
            # How do we do this? In the second for-loop below, we have the variable "letter_b". On this second iteration, 
            # we need to look at a range of indexes equivalent to the length of the STR "AGATC" to see if the next letters 
            # in our DNA sequence follow the same pattern of letters. To do this, we choose the (previous letter's (the 
            # variable "letter") index + the STR_length) as our starting point. 
            # 
            # This effectively makes it so the starting point of our next chunk analysis is the letter after our first chunk's 
            # letters. So if our first chunk analysis started at index 0 and ended at index 4, our next chunk analysis would 
            # start at index 5 and end at index 9; The chunk after that would start at index 10 and end at index 14. And so on 
            # until we reach a chunk that isn't this STR. "AGATC" is analyzed in chunks of 5 letters every loop. 
            
            # We can use this same algorithm for all STRs we want to look at in a DNA file.

            # Step through each letter in our DNA file.
            for letter in range(0, len(sequence)):
                # Look at a chunk of letters according to the length of our STR, starting at this letter.
                letter_chunk = sequence[letter:letter+STR_length]
                # If the letter_chunk does equal our STR, that means we've found at least one STR sequence 
                # for this STR in the DNA file and we can now start looking for this specific STR throughout
                # the rest of the file.
                if (letter_chunk == STR):
                    # We count this STR sequence to keep track of how many times this STR repeats.
                    current_count += 1
                    # Now we need to determine if this letter_chunk, AKA STR, repeats, and if so, we need to 
                    # count those repeats. We will execute this for-loop every chunk, for the entire DNA sequence,
                    # until we hit a chunk that no longer equals our STR.
                    for letter_b in range((letter+STR_length), len(sequence), STR_length):
                        # Look at a chunk of letters according to the length of our STR.
                        STR_chunk = sequence[letter_b:letter_b+STR_length]
                        # If the chunk we are looking at no longer equals our STR, we break out of this for-loop
                        # block and resume our immediate outer loop to check the rest of the letters of the DNA 
                        # sequence for this STR.
                        if STR_chunk != STR:
                            # Before we break to check the rest of the sequence for more of the same STR
                            # repetitions (starting from the letter index we left off of in our immediate 
                            # outer loop), we need to see if our current_count is our largest STR repetition
                            # for this STR thus far. If it is, we set it to our largest_count.
                            if current_count > largest_count:
                                largest_count = current_count
                            # Reset current_count to zero in preparation for our next repetition count for this 
                            # STR.
                            current_count = 0
                            break
                        # Otherwise, we will have another STR repetition, so we add it to current_count.
                        else:
                            current_count += 1

            # Append every largest repetition for each STR to a list of largest STR repetitions.
            # Make it a string for easy comparison later.
            analyzed_sequence.append(str(largest_count))

        # Break out of while loop if the size of our analyzed list is the same as the number of STRs we've
        # had to analyze, since we will know the job is finished if both sizes are equal.
        if len(analyzed_sequence) == len(STR_list):
            break
    
    # Return the results for comparison. 
    return analyzed_sequence


# Allows us to format our rows from the csv file into a usable format, a list of lists.
def format_rows(csv_arg):
    # A list to collectively hold every person and their values.
    database_list = []
    # A list to hold a single person and their values.
    person_info = []

    # Open the dictionary to read.
    dictionary = open(csv_arg, "r")

    # Populate a list of names and their STR values per each row from the dictionary.
    for row in dictionary:
        # Ignore the first row because it contains headers.
        if row[0][0] != "n":
            # Split the row of the csv_arg by comma values and add to list.
            person_info = row.split(",")
            # Remove last element from list so you can remove the \n that won't go away any other way.
            last_element = person_info[len(person_info) - 1]
            person_info.remove(last_element)
            # Initialize a new value to put what will be our truncated last_element in.
            new_last_element = ""
            # Iterate through each letter in last_element and add each letter to new_last_element up until
            # (not including) the second to last character.
            for c in last_element:
                if len(new_last_element) < len(last_element) - 1:
                    new_last_element += c
            # Finally, add the new element that does not include \n to the person_info list.
            person_info.append(new_last_element)
        # Add our newly formatted item to our database_list.
        database_list.append(person_info)

    # Remove the first empty header row from the list since retaining it will mess with our comparison later.
    database_list.remove(database_list[0])

    return database_list


# Compares the STR database dictionary with our analyzed DNA sequence, and returns the name of the match if there is
# a match.
def compare(csv_arg, analyzed_sequence):
    # Format our CSV file of STR values.
    database_list = format_rows(csv_arg)
    # A list to store the names of people whose DNA we are analyzing.
    name_list = []

    # Transfer database_list names to our name_list. The name indexes in name_list will coincide with the database_list
    # indexes, since we are loading the names one by one from the database_list into our name_list. We can then use
    # these indexes (since they coincide with each other) to print out the appropriate name when we find a match
    # between a row in the database_list and our analyzed_sequence.
    for row in database_list:
        # Add the name at this row, which will always be at index 0, to the name list.
        name_list.append(row[0])
        # Remove the name from the database_list since we no longer need it (as each database_list index will now be matched 
        # to the appropriate name_list index). In database_list, we'll have just the values of the STRs without names, and we 
        # can then compare those pure values to our analyzed_sequence list.
        del row[0]

    # Compare our list of STR repetitions to our analyzed_sequence repetitions.
    row_count = 0
    for row in database_list:
        # If our row in the database has the same repetitions as our analyzed_sequence, print the name of the
        # person that is associated with that row in our original csv file.
        if row == analyzed_sequence:
            print(f"{name_list[row_count]}")
            return 0
        # Count which row we're on so we can locate the name associated
        # with that row.
        row_count += 1

    # Otherwise, if there is no DNA match, we return "No match".
    print("No match")


main()