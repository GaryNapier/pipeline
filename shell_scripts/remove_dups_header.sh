#! /bin/bash

# FILENAME
# remove_dups_header.sh

set -e
set -u
set -o pipefail

# ------------------------------------------------------------------------------

# INFORMATION AND DESCRIPTION
# Remove duplicates in a file with a header e.g. a csv file
# Highly recommended to run remove_hidden.sh first so that hidden characters are removed
# see https://stackoverflow.com/questions/14562423/is-there-a-way-to-ignore-header-lines-in-a-unix-sort

# Arguments to script:
# file

# Input:
# Any file with header and duplicates
# Example:
# $ cat -A <file>
# drug,gene,mutation,samples$
# isoniazid,katG,c.2223A>G,ERR2864254$
# isoniazid,katG,c.2223A>G,ERR2864254$
# isoniazid,katG,c.2223A>G,ERR760749$
# isoniazid,katG,c.2223A>G,ERR760749$
# isoniazid,katG,c.2223A>G,SRR5073823$
# isoniazid,katG,c.2223A>G,SRR5073823$
# isoniazid,katG,p.Arg484His,ERR211998$
# isoniazid,katG,p.Arg484His,ERR211998$

# Main steps:
# Separate the header from the rest of the file with the head and tail commands, using the parentheses to keep commands separate (don't know how this works)
# sort | uniq on the tail 
# store in tmp file
# move back to original file name

# Output:
# de-duped file
# Example:
# $ cat -A <file>
# drug,gene,mutation,samples$
# isoniazid,katG,c.2223A>G,ERR2864254$
# isoniazid,katG,c.2223A>G,ERR760749$
# isoniazid,katG,c.2223A>G,SRR5073823$
# isoniazid,katG,p.Arg484His,ERR211998$

# RUN EXAMPLE
# remove_dups_header.sh <file>

# ------------------------------------------------------------------------------

# Setup

# Variables


# Directories


# Files and suffixes/prefixes


# Parameters/variables


# Print all args and params/variables
printf "\n"
echo "Parameters:"
printf "\n"
echo ${file}
printf "\n"

echo "Removing dups from ${file}"
(head -n 1 ${file} && tail -n +2 ${file} | sort | uniq) > tmp.csv && mv tmp.csv ${file}