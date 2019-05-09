MASH_DIR=/home/gandhijs/isolates_flagged_mcr-1/mash_output

if [ ! -f $MASH_DIR/master_mash.csv ]
    then touch $MASH_DIR/master_mash.csv
    echo "Reference-ID,Query-ID,MASH-distance,P-value,Matching-hashes" >> $MASH_DIR/master_mash.csv 
fi

for tab_file in $MASH_DIR/*.tab
do
    # sort the top hits by the p-value
    sort -o $tab_file -gk3 $tab_file
    # get the file name
    filename=$(basename -- "$tab_file")
    # grab the sample name
    SAMPLE_NAME="${filename%%_*}"
    echo writing $SAMPLE_NAME to master_mash.csv
    # appends the top mash hit to the master file
    cat $tab_file | awk 'NR==1' | tr "\\t" "," >> $MASH_DIR/master_mash.csv
done
