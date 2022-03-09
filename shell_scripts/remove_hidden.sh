#! /bin/bash

# FILENAME
# remove_hidden.sh

set -e
set -u
set -o pipefail

# ------------------------------------------------------------------------------

# INFORMATION AND DESCRIPTION
# Remove stupid "^M" ('carriage return') character at the end of each line - often stupidly put there by excel/windows, or, as I found out today, python's writer.writerow()
# https://unix.stackexchange.com/questions/32001/what-is-m-and-how-do-i-get-rid-of-it

# Arguments to script:
# file

# Input:
# Any file
# Example:
# $ cat -A <file>
# drug,gene,mutation,samples^M$
# isoniazid,katG,c.2223A>G,ERR2864254^M$
# isoniazid,katG,c.2223A>G,ERR760749^M$
# isoniazid,katG,c.2223A>G,SRR5073823^M$
# isoniazid,katG,p.Arg484His,ERR211998^M$

# Main steps:
# Run the sed command on the file, removing the stupid characters

# Output:
# Clean file
# Example:
# $ cat -A <file>
# drug,gene,mutation,samples$
# isoniazid,katG,c.2223A>G,ERR2864254$
# isoniazid,katG,c.2223A>G,ERR760749$
# isoniazid,katG,c.2223A>G,SRR5073823$
# isoniazid,katG,p.Arg484His,ERR211998$


# RUN EXAMPLE
# remove_hidden <file>

# ------------------------------------------------------------------------------

# Setup

# Variables


# Directories


# Files and suffixes/prefixes
file=${1}

# Parameters/variables


# Print all args and params/variables
printf "\n"
echo "Parameters:"
printf "\n"
echo ${file}
printf "\n"

echo "Removing hidden characters from ${file}"
sed -i -e 's/\r//g' ${file}